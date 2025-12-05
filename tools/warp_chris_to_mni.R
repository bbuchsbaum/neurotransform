suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  }
  library(neurotransform)
  library(rpyANTs)
  library(reticulate)
})

chris_path <- "inst/extdata/chris/chris_t1.nii.gz"
template_path <- "inst/extdata/template/MNI152NLin2009cAsym_T1w_2mm_brain.nii.gz"
out_dir <- "inst/extdata/chris/ants"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(chris_path), file.exists(template_path))
stopifnot(ants_available())

cat("Running ANTs SyN registration (Chris -> MNI)...\n")
ants <- load_ants()
prefix <- file.path(out_dir, "reg_")
res <- ants_registration(
  fixed = template_path,
  moving = chris_path,
  type_of_transform = "SyN",
  outprefix = prefix,
  write_composite_transform = TRUE,
  verbose = TRUE
)

# Inverse composite H5 (template -> chris)
inv_list <- reticulate::py_to_r(res$invtransforms)
inv_h5 <- inv_list[1]
stopifnot(!is.na(inv_h5), file.exists(inv_h5))
cat("Inverse composite transform:", inv_h5, "\n")

# High-level load + resample
chris_vol <- read_image(chris_path)
template_vol <- read_image(template_path)
morph <- read_transform(inv_h5)  # type inferred as ants_h5

cat("Resampling Chris T1 into MNI space...\n")
out_vol <- resample_to(chris_vol, target = template_vol, transform = morph,
                       method = "linear", modulate = "none")

out_path <- file.path(out_dir, "chris_in_mni.nii.gz")
write_image(out_vol, out_path)
cat("Wrote ", out_path, "\n", sep = "")
