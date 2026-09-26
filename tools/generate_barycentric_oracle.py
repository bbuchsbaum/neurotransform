"""Generate ordinary Workbench BARYCENTRIC fixtures, without neurotransform.

python3 tools/generate_barycentric_oracle.py /path/to/wb_command output-directory
Requires numpy. The output directory must not exist. GIFTI inputs, outputs,
commands, binary/input hashes and version are retained with the JSON fixture.
"""
import argparse
import base64
import hashlib
import json
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET
import zlib

import numpy as np


def write_gifti(path, arrays):
    root = ET.Element("GIFTI", Version="1.0", NumberOfDataArrays=str(len(arrays)))
    meta = ET.SubElement(root, "MetaData")
    md = ET.SubElement(meta, "MD")
    ET.SubElement(md, "Name").text = "AnatomicalStructurePrimary"
    ET.SubElement(md, "Value").text = "CortexLeft"
    ET.SubElement(root, "LabelTable")
    for intent, values in arrays:
        values = np.asarray(values)
        integer = intent == "NIFTI_INTENT_TRIANGLE"
        attrs = dict(Intent=intent, DataType="NIFTI_TYPE_INT32" if integer else "NIFTI_TYPE_FLOAT32",
                     ArrayIndexingOrder="RowMajorOrder", Dimensionality=str(values.ndim),
                     Encoding="ASCII", Endian="LittleEndian", ExternalFileName="", ExternalFileOffset="")
        attrs.update({"Dim" + str(i): str(n) for i, n in enumerate(values.shape)})
        array = ET.SubElement(root, "DataArray", attrs)
        ET.SubElement(array, "MetaData")
        ET.SubElement(array, "Data").text = " ".join(
            str(int(v)) if integer else format(float(v), ".9g") for v in values.ravel())
    ET.ElementTree(root).write(path, encoding="UTF-8", xml_declaration=True)


def read_gifti(path):
    result = []
    for array in ET.parse(path).getroot().findall("DataArray"):
        dims = tuple(int(array.get("Dim" + str(i))) for i in range(int(array.get("Dimensionality"))))
        dtype = {"NIFTI_TYPE_FLOAT32": "f4", "NIFTI_TYPE_INT32": "i4"}[array.get("DataType")]
        dtype = ("<" if array.get("Endian") == "LittleEndian" else ">") + dtype
        data = array.findtext("Data", "")
        if array.get("Encoding") == "ASCII":
            values = np.fromstring(data, sep=" ", dtype=dtype)
        else:
            raw = base64.b64decode(data)
            if array.get("Encoding") == "GZipBase64Binary":
                raw = zlib.decompress(raw, wbits=47)
            elif array.get("Encoding") != "Base64Binary":
                raise ValueError("Unsupported GIFTI encoding")
            values = np.frombuffer(raw, dtype=dtype)
        order = "C" if array.get("ArrayIndexingOrder") == "RowMajorOrder" else "F"
        result.append(values.reshape(dims, order=order))
    return result


def sphere(level):
    vertices = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
    faces = [[0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4],
             [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5]]
    for _ in range(level):
        edges, refined = {}, []
        def midpoint(a, b):
            key = tuple(sorted((a, b)))
            if key not in edges:
                v = np.asarray(vertices[a]) + vertices[b]
                edges[key] = len(vertices)
                vertices.append((v / np.linalg.norm(v)).tolist())
            return edges[key]
        for a, b, c in faces:
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            refined.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]])
        faces = refined
    return np.asarray(vertices, dtype=float), np.asarray(faces, dtype=np.int32)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbench", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    wb, out = args.workbench.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    version = subprocess.check_output([str(wb), "-version"], text=True)
    commands, cases = [], []
    for name, level, angle in [("octahedron-boundaries", 0, 0),
                                ("octahedron-asymmetric", 0, 0.37),
                                ("irregular-sphere", 1, -0.23)]:
        vertices, faces = sphere(level)
        query, query_faces = sphere(2)
        if level:
            vertices = vertices @ np.array([[1, .17, .09], [.04, 1.1, .12], [.02, -.07, .9]])
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                             [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        query = query @ rotation
        vertices = (100 * vertices / np.linalg.norm(vertices, axis=1)[:, None]).astype(np.float32)
        query = (100 * query / np.linalg.norm(query, axis=1)[:, None]).astype(np.float32)
        folder = out / name
        folder.mkdir()
        write_gifti(folder / "source.surf.gii", [("NIFTI_INTENT_POINTSET", vertices), ("NIFTI_INTENT_TRIANGLE", faces)])
        write_gifti(folder / "query.surf.gii", [("NIFTI_INTENT_POINTSET", query), ("NIFTI_INTENT_TRIANGLE", query_faces)])
        # A basis of impulses identifies every interpolation weight, so a
        # constant-field false positive or a coincidental scalar match cannot pass.
        write_gifti(folder / "input.func.gii", [("NIFTI_INTENT_SHAPE", row) for row in np.eye(len(vertices))])
        command = [str(wb), "-metric-resample", str(folder / "input.func.gii"),
                   str(folder / "source.surf.gii"), str(folder / "query.surf.gii"),
                   "BARYCENTRIC", str(folder / "output.func.gii")]
        run = subprocess.run(command, text=True, capture_output=True, check=True)
        commands.append(dict(argv=command, exit_code=run.returncode, stdout=run.stdout, stderr=run.stderr))
        expected = np.column_stack(read_gifti(folder / "output.func.gii"))
        cases.append(dict(name=name, vertices=vertices.tolist(), faces=faces.tolist(),
                          query=query.tolist(), weights=expected.tolist()))
    fixture = dict(method="BARYCENTRIC", tolerance=2e-6, version=version,
                   binary_sha256=sha(wb), generator_sha256=sha(__file__), cases=cases,
                   commands=commands, files={str(p.relative_to(out)): sha(p) for p in sorted(out.glob("*/*.gii"))})
    (out / "oracle.json").write_text(json.dumps(fixture, indent=2) + "\n")
    print(out / "oracle.json")


if __name__ == "__main__":
    main()
