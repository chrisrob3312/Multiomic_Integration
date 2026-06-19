## Shared helpers for the REDIAL pipeline.
## Sourced by every 0x_ script.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
})

PATHS <- list(
  data    = "data/synthetic_redial.rds",
  results = "results"
)

load_redial <- function(path = PATHS$data) {
  if (!file.exists(path)) stop("Run data/synthetic_redial.R first to generate ", path)
  readRDS(path)
}

## Residualize a feature matrix (samples x features) against covariates.
## Use BEFORE handing matrices to MOFA2 so subtype/age/sex/blast% don't dominate
## the latent space. Missing entries are preserved as NA.
residualize <- function(mat, covariates) {
  stopifnot(nrow(mat) == nrow(covariates))
  out <- mat
  X <- model.matrix(~ ., data = covariates)
  for (j in seq_len(ncol(mat))) {
    y <- mat[, j]
    keep <- !is.na(y)
    if (sum(keep) < ncol(X) + 2) next
    fit <- lm.fit(X[keep, , drop = FALSE], y[keep])
    out[keep, j] <- fit$residuals
  }
  out
}

## Build the MOFA-ready views list. Each element is features x samples.
build_views <- function(redial, sample_ids, residualize_against = NULL) {
  views <- list(
    methylation = t(redial$omics$methylation[sample_ids, , drop = FALSE]),
    rna         = t(redial$omics$rna[sample_ids, , drop = FALSE]),
    metabolome  = t(redial$omics$metabolome[sample_ids, , drop = FALSE]),
    cnv         = t(redial$omics$cnv[sample_ids, , drop = FALSE])
  )
  if (!is.null(residualize_against)) {
    cov_df <- redial$clinical[sample_ids, residualize_against, drop = FALSE]
    views <- lapply(views, function(v) t(residualize(t(v), cov_df)))
  }
  views
}

## Fit MOFA on a list of views with a grouping vector.
## Falls back to ungrouped if `groups` is NULL.
fit_mofa <- function(views, groups = NULL, num_factors = 10, seed = 1) {
  if (!requireNamespace("MOFA2", quietly = TRUE))
    stop("MOFA2 not installed. BiocManager::install('MOFA2')")
  if (is.null(groups)) {
    obj <- MOFA2::create_mofa(views)
  } else {
    samples_per_group <- split(colnames(views[[1]]), groups)
    grouped <- lapply(views, function(mat) {
      lapply(samples_per_group, function(s) mat[, s, drop = FALSE])
    })
    obj <- MOFA2::create_mofa(grouped)
  }
  data_opts  <- MOFA2::get_default_data_options(obj)
  model_opts <- MOFA2::get_default_model_options(obj)
  train_opts <- MOFA2::get_default_training_options(obj)
  model_opts$num_factors <- num_factors
  train_opts$seed <- seed
  train_opts$convergence_mode <- "medium"
  obj <- MOFA2::prepare_mofa(obj,
                             data_options = data_opts,
                             model_options = model_opts,
                             training_options = train_opts)
  tmp <- tempfile(fileext = ".hdf5")
  MOFA2::run_mofa(obj, outfile = tmp, use_basilisk = TRUE)
}

## Pull factor scores (samples x factors) from a fitted MOFA object.
get_factor_scores <- function(mofa) {
  Z <- MOFA2::get_factors(mofa, factors = "all", as.data.frame = FALSE)
  if (is.list(Z)) do.call(rbind, Z) else Z
}

## Cox PH on factor scores for time-to-event.
## Returns a tidy data.frame: factor, HR, CI, p, FDR.
cox_on_factors <- function(factor_scores, time, event, adjust_df = NULL) {
  res <- lapply(colnames(factor_scores), function(fac) {
    df <- data.frame(time = time, event = event, f = factor_scores[, fac])
    if (!is.null(adjust_df)) df <- cbind(df, adjust_df)
    fit <- survival::coxph(survival::Surv(time, event) ~ ., data = df)
    s <- summary(fit)$coefficients["f", , drop = FALSE]
    ci <- summary(fit)$conf.int["f", c("lower .95", "upper .95")]
    data.frame(factor = fac,
               HR = s[, "exp(coef)"],
               lo95 = ci[1], hi95 = ci[2],
               p = s[, "Pr(>|z|)"])
  })
  out <- do.call(rbind, res)
  out$FDR <- p.adjust(out$p, "BH")
  out[order(out$p), ]
}

## Regress each factor on ancestry x ADI (+ adjustments).
## Returns tidy effect estimates.
factor_assoc <- function(factor_scores, clinical, exposure, adjust = NULL) {
  rhs <- if (is.null(adjust)) exposure else c(exposure, adjust)
  fml <- as.formula(paste("f ~", paste(rhs, collapse = " + ")))
  do.call(rbind, lapply(colnames(factor_scores), function(fac) {
    df <- cbind(f = factor_scores[, fac], clinical)
    fit <- lm(fml, data = df)
    s <- summary(fit)$coefficients
    s <- s[grep(exposure, rownames(s)), , drop = FALSE]
    data.frame(factor = fac, term = rownames(s),
               estimate = s[, "Estimate"], se = s[, "Std. Error"],
               p = s[, "Pr(>|t|)"], row.names = NULL)
  }))
}

## Cluster samples on factor scores -> high-risk / low-risk labels.
## Two-class k-means, oriented so "high-risk" has worse OS in the supplied
## time/event vectors.
risk_clusters <- function(factor_scores, time, event, seed = 1) {
  set.seed(seed)
  km <- kmeans(scale(factor_scores), centers = 2, nstart = 25)
  cl <- km$cluster
  fit1 <- survival::survfit(survival::Surv(time, event) ~ cl)
  med <- summary(fit1)$table[, "median"]
  hi <- if (which.min(med) == 1) 1 else 2
  ifelse(cl == hi, "high_risk", "low_risk")
}

dir.create(PATHS$results, showWarnings = FALSE, recursive = TRUE)
