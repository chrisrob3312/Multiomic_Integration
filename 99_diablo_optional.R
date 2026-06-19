## Optional supervised confirmatory: mixOmics DIABLO.
## Run this ONLY if the MOFA-derived high-risk vs low-risk split (scripts
## 01/02) is weak or you want a directly outcome-supervised view.
##
## Two configurations:
##   A. Y = high_risk vs low_risk label from script 01 (MOFA-derived).
##   B. Y = relapse_category (4-class) — direct supervision on the timing.
##
## Caveats: DIABLO does not natively handle block-missing samples; we
## complete-case on the four omics views. At n ~ 20-40 per class do
## leave-one-out CV instead of 5-fold.

source("R/utils.R")
if (!requireNamespace("mixOmics", quietly = TRUE)) stop("install.packages('mixOmics')")

redial <- load_redial()
cl <- redial$clinical

risk_label_path <- file.path(PATHS$results, "01_risk_by_stratum.csv")
if (!file.exists(file.path(PATHS$results, "mofa_overall.rds")))
  stop("Run 01_mofa_overall.R first.")
mofa_overall <- readRDS(file.path(PATHS$results, "mofa_overall.rds"))
Z <- get_factor_scores(mofa_overall)[cl$sample_id, , drop = FALSE]
risk_label <- risk_clusters(Z, cl$os_time, cl$os_event)

## Complete-case across the four omics blocks.
have <- Reduce(`&`, lapply(redial$omics, function(m) complete.cases(m)))
ids  <- cl$sample_id[have]
X <- list(
  methylation = redial$omics$methylation[ids, ],
  rna         = redial$omics$rna[ids, ],
  metabolome  = redial$omics$metabolome[ids, ],
  cnv         = redial$omics$cnv[ids, ]
)

design <- matrix(0.1, nrow = length(X), ncol = length(X),
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

run_diablo <- function(Y, tag, ncomp = 3) {
  fit <- mixOmics::block.splsda(
    X = X, Y = Y, ncomp = ncomp, design = design,
    keepX = lapply(X, function(m) rep(min(50, ncol(m)), ncomp))
  )
  perf <- mixOmics::perf(
    fit, validation = "loo", folds = NULL,
    dist = "centroids.dist", progressBar = FALSE
  )
  saveRDS(list(fit = fit, perf = perf),
          file.path(PATHS$results, sprintf("diablo_%s.rds", tag)))
  cat(sprintf("DIABLO %s: error (centroids) =\n", tag))
  print(perf$error.rate)
  invisible(fit)
}

run_diablo(factor(risk_label[have]),                 "risk_label")
run_diablo(factor(cl[ids, "relapse_category"]),      "relapse_category")
