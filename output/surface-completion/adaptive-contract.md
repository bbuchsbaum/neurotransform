# Native adaptive area contract

Pinned comparator: Workbench 2.2.1, revision
`01164ffa47f2778088bd6ff472ec9cd9a57f5b42`,
`src/Files/SurfaceResamplingHelper.cxx::computeWeightsAdapBaryArea`.
The support rule and equations were independently reviewed before implementation.

Let F be ordinary forward weights (target by source), R the independently built
reverse weights, and G=transpose(R). No ROI enters these constructions.
For target j choose A[j,]=F[j,] if support(G[j,]) is a subset of support(F[j,]);
otherwise choose G[j,] in its entirety. This is not a blend or a comparison of
support cardinalities. Support means precisely positive stored weight.

For positive anatomical source areas a and target areas b, compute
c[i]=sum_j b[j] A[j,i], then H[j,i]=b[j] A[j,i] a[i]/c[i] on supported entries.
Apply source ROI q[i]=1[roi[i]>0] only now; h[j]=sum_i H[j,i] q[i], and
W[j,i]=H[j,i] q[i]/h[j] for h[j]>0. Apply target ROI after this normalization.
An implementation may store unmasked W0=H/rowsum(H), then apply source ROI and
renormalize rows: this gives the same W while retaining a reusable operator.

With no target exclusions and finite data, the exact identity is
sum_j h[j] (W x)[j] = sum_{represented i} a[i] q[i] x[i].
The effective target measure h is generally NOT the supplied b. Final W is
constant preserving; supplied-area integral conservation is not claimed.
Missing-data omission changes the operator and invalidates this fixed identity.
Column normalization is not this algorithm.

Area objects carry exact ordered geometry identity, units and provenance.
All areas must be finite and strictly positive; zero areas are rejected.
An anatomical mesh must have the same ordered faces and vertex count as its
registered geometry. Spherical areas are never silently substituted for anatomy.

ADAP_BARY_AREA qualification remains separate from ordinary interpolation.
Near-zero forward/reverse support changes can switch entire adaptive rows.
The existing ordinary template parity gate failed; adaptive implementation is
therefore exposed only with explicit experimental opt-in until its own gates
pass. No silent fallback or tolerance-based support pruning is allowed.

For explicit omission, a separate data-dependent identity uses the per-column
measure h0*finite_weight_mass and includes only finite represented source
values. It does not supply an adjoint of the original fixed operator. Filter
unavailable rows before forming the integral, and do not apply target exclusion
when asserting full represented-source conservation. Alternative normalization
modes do not expose an effective_target_area value.
