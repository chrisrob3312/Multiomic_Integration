## LUCIDus quasi-mediation: exposure -> latent omics clusters -> relapse.
## Three exposure configurations, each fit independently:
##   (a) germline risk-variant burden
##   (b) ADI quartile (high deprivation = q >= 3 vs q <= 2)
##   (c) Amerindigenous-ancestry indicator (AMI vs other)
##
## Omics layers are passed as raw selected features per modality (LUCIDus'
## native input) — the script trims each modality to the top-variance features
## first so the EM converges on a clinical-scale cohort.

source("R/utils.R")
if (!requireNamespace("LUCIDus", quietly = TRUE))
  stop("install.packages('LUCIDus')")

redial <- load_redial()
cl <- redial$clinical
ids <- cl$sample_id

## Trim each modality to top-variance features (cheap stand-in for the real
## per-modality biology-aware pre-screen you'll do upstream of LUCIDus).
top_var <- function(mat, k) {
  v <- apply(mat, 2, var, na.rm = TRUE)
  mat[, order(-v)[seq_len(min(k, ncol(mat)))], drop = FALSE]
}
Z_list <- list(
  methylation = top_var(redial$omics$methylation[ids, ], 200),
  rna         = top_var(redial$omics$rna[ids, ],         200),
  metabolome  = top_var(redial$omics$metabolome[ids, ],   80),
  cnv         = top_var(redial$omics$cnv[ids, ],         100)
)

Y <- cl$relapsed
COV <- cl[ids, c("subtype", "age", "sex", "blast")]

## NOTE: for lucid_model = "parallel", K must be a vector of length(Z) (one K
## per omics block). For "early" or "serial", K is a single integer.
fit_lucid <- function(exposure_vec, tag, model = "parallel", K = NULL) {
  if (is.null(K)) K <- if (model == "parallel") rep(2, length(Z_list)) else 2
  fit <- LUCIDus::estimate_lucid(
    G = matrix(exposure_vec, ncol = 1, dimnames = list(NULL, tag)),
    Z = Z_list,
    Y = Y,
    CoY = model.matrix(~ ., COV)[, -1, drop = FALSE],
    family = "binary",
    lucid_model = model,
    K = K,
    seed = 1
  )
  saveRDS(fit, file.path(PATHS$results, sprintf("lucid_%s.rds", tag)))
  fit
}

fit_germ <- fit_lucid(redial$germline$burden,            "germline_burden")
fit_adi  <- fit_lucid(as.integer(cl$adi_q >= 3),         "adi_high")
fit_ami  <- fit_lucid(as.integer(cl$ancestry == "AMI"),  "ancestry_AMI")

## Bootstrap inference on the germline arm (the most expensive — comment out
## the others if you're iterating).
if (requireNamespace("LUCIDus", quietly = TRUE)) {
  boot_germ <- LUCIDus::boot_lucid(fit_germ, R = 200, seed = 1)
  saveRDS(boot_germ, file.path(PATHS$results, "lucid_germline_boot.rds"))
}

cat("LUCIDus runs written to", PATHS$results, "\n")
