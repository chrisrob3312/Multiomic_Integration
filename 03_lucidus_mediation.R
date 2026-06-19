## LUCIDus quasi-mediation: exposure(s) -> latent omics clusters -> relapse.
##
## NOTE on LUCIDus G (exposures): G can be CONTINUOUS or CATEGORICAL or a MIX
## of multiple columns. The "categorical" thing in LUCIDus is the latent
## omics cluster (K classes), NOT the exposure. So for germline you can pass:
##   - mode = "snps"   : multi-column matrix of 0/1/2 dosages over the curated
##                       B-ALL risk-variant panel (recommended primary).
##   - mode = "prs"    : single weighted polygenic risk score column.
##   - mode = "burden" : single count of risk-panel carrier variants.
## Use the same `fit_lucid()` wrapper for ADI and ancestry exposures too.
##
## Omics layers are passed as raw selected features per modality. Methylation
## / RNA / metabolome / CNV-from-RNAseq are kept as four SEPARATE Z-blocks so
## LUCIDus' parallel-mode K = c(K_m, K_r, K_b, K_c) per-block clustering is
## interpretable per modality.

source("R/utils.R")
if (!requireNamespace("LUCIDus", quietly = TRUE)) stop("install.packages('LUCIDus')")

redial <- load_redial()
cl <- redial$clinical
ids <- cl$sample_id

top_var <- function(mat, k) {
  v <- apply(mat, 2, var, na.rm = TRUE)
  mat[, order(-v)[seq_len(min(k, ncol(mat)))], drop = FALSE]
}
Z_list <- list(
  methylation = top_var(redial$omics$methylation[ids, ], 200),
  rna         = top_var(redial$omics$rna[ids, ],         200),
  metabolome  = top_var(redial$omics$metabolome[ids, ],   80),
  cnv         = top_var(redial$omics$cnv[ids, ],         100)  # CNV-from-RNAseq layer
)

Y <- cl$relapsed
COV <- cl[ids, c("subtype", "age", "sex", "blast")]

## --- Build the G matrix for a given exposure mode --------------------------
build_G <- function(mode = c("snps", "prs", "burden")) {
  mode <- match.arg(mode)
  switch(mode,
    snps   = redial$germline$matrix[ids, redial$germline$risk_panel, drop = FALSE],
    prs    = matrix(redial$germline$prs[match(ids, rownames(redial$germline$matrix))],
                    ncol = 1, dimnames = list(NULL, "germline_PRS")),
    burden = matrix(redial$germline$burden[match(ids, rownames(redial$germline$matrix))],
                    ncol = 1, dimnames = list(NULL, "germline_burden"))
  )
}

fit_lucid <- function(G, tag, model = "parallel", K = NULL) {
  if (is.null(K)) K <- if (model == "parallel") rep(2, length(Z_list)) else 2
  fit <- LUCIDus::estimate_lucid(
    G = G, Z = Z_list, Y = Y,
    CoY = model.matrix(~ ., COV)[, -1, drop = FALSE],
    family = "binary",
    lucid_model = model,
    K = K,
    seed = 1
  )
  saveRDS(fit, file.path(PATHS$results, sprintf("lucid_%s.rds", tag)))
  fit
}

## --- Germline arm: three encodings, primary is multi-SNP -------------------
fit_snps   <- fit_lucid(build_G("snps"),   "germline_snps")
fit_prs    <- fit_lucid(build_G("prs"),    "germline_prs")
fit_burden <- fit_lucid(build_G("burden"), "germline_burden")

## --- ADI and ancestry exposure arms ----------------------------------------
fit_adi <- fit_lucid(
  matrix(as.integer(cl$adi_q >= 3), ncol = 1, dimnames = list(NULL, "adi_high")),
  "adi_high"
)
fit_ami <- fit_lucid(
  matrix(as.integer(cl$ancestry == "AMI"), ncol = 1, dimnames = list(NULL, "ancestry_AMI")),
  "ancestry_AMI"
)

## --- Bootstrap inference (start with the multi-SNP germline arm) -----------
if (requireNamespace("LUCIDus", quietly = TRUE)) {
  boot_snps <- LUCIDus::boot_lucid(fit_snps, R = 200, seed = 1)
  saveRDS(boot_snps, file.path(PATHS$results, "lucid_germline_snps_boot.rds"))
}

cat("LUCIDus runs written to", PATHS$results, "\n")
