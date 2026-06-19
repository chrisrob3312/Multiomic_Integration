## MRD-stratified MOFA runs: one for MRD-negative, one for MRD-positive.
##
## Grouping: tumor_subtype (MOFA absorbs subtype heterogeneity instead of
## forcing per-subtype refits at small n).
## Residualized: technical noise only — age, sex, blast%.
##
## NOTE on clustering: MOFA itself does NOT cluster samples. It produces
## factor scores (samples x factors). We then run k-means (k=2) on the
## factor scores to derive a high-risk vs low-risk LABEL, oriented by the
## clinical outcome you care about.
##
## MRD-neg arm:
##   - Outcome of interest = relapse y/n  (which factors / pathways drive
##     relapse despite MRD-negativity).
##   - Risk label oriented by relapse y/n.
## MRD-pos arm:
##   - Outcome of interest = OS time-to-event (more events at this baseline
##     risk).
##   - Risk label oriented by OS.

source("R/utils.R")
redial <- load_redial()
cl <- redial$clinical

adj_cov <- c("age", "sex", "blast")

run_one <- function(label, ids, orient_by = c("relapse", "os")) {
  orient_by <- match.arg(orient_by)
  views  <- build_views(redial, ids, residualize_against = adj_cov)
  groups <- as.character(cl[ids, "subtype"])
  mofa   <- fit_mofa(views, groups = groups, num_factors = 8, seed = 1)

  Z <- get_factor_scores(mofa)[ids, , drop = FALSE]
  sub_cl <- cl[ids, , drop = FALSE]

  ## OS Cox per factor (always reported)
  cox_tab <- cox_on_factors(
    Z, time = sub_cl$os_time, event = sub_cl$os_event,
    adjust_df = sub_cl[, c("age", "sex", "subtype")]
  )

  ## Relapse y/n per-factor logistic (the key story for MRD-neg)
  relapse_tab <- do.call(rbind, lapply(colnames(Z), function(fac) {
    df <- data.frame(y = sub_cl$relapsed, f = Z[, fac],
                     sub_cl[, c("age", "sex", "subtype")])
    fit <- glm(y ~ ., data = df, family = binomial())
    s <- summary(fit)$coefficients["f", , drop = FALSE]
    data.frame(factor = fac, OR = exp(s[, "Estimate"]),
               p = s[, "Pr(>|z|)"], row.names = NULL)
  }))
  relapse_tab$FDR <- p.adjust(relapse_tab$p, "BH")

  ## Ancestry / ADI exposure effects on each factor
  anc_assoc <- factor_assoc(
    Z, sub_cl, exposure = "ancestry",
    adjust = c("adi_q", "subtype", "age", "sex", "blast")
  )
  adi_assoc <- factor_assoc(
    Z, sub_cl, exposure = "adi_q",
    adjust = c("ancestry", "subtype", "age", "sex", "blast")
  )

  ## Risk clustering oriented by the chosen outcome.
  risk_label <- if (orient_by == "relapse") {
    set.seed(1)
    km <- kmeans(scale(Z), centers = 2, nstart = 25)$cluster
    rate <- tapply(sub_cl$relapsed, km, mean)
    hi <- as.integer(names(which.max(rate)))
    ifelse(km == hi, "high_risk", "low_risk")
  } else {
    risk_clusters(Z, sub_cl$os_time, sub_cl$os_event)
  }

  risk_by_stratum <- sub_cl |>
    mutate(risk_label = risk_label) |>
    count(ancestry, adi_q, risk_label) |>
    group_by(ancestry, adi_q) |>
    mutate(frac_high = n / sum(n)) |>
    filter(risk_label == "high_risk") |>
    ungroup()

  saveRDS(mofa,            file.path(PATHS$results, sprintf("mofa_%s.rds", label)))
  write.csv(cox_tab,         file.path(PATHS$results, sprintf("02_cox_%s.csv", label)),            row.names = FALSE)
  write.csv(relapse_tab,     file.path(PATHS$results, sprintf("02_relapse_logit_%s.csv", label)), row.names = FALSE)
  write.csv(anc_assoc,       file.path(PATHS$results, sprintf("02_ancestry_%s.csv", label)),       row.names = FALSE)
  write.csv(adi_assoc,       file.path(PATHS$results, sprintf("02_adi_%s.csv", label)),            row.names = FALSE)
  write.csv(risk_by_stratum, file.path(PATHS$results, sprintf("02_risk_by_stratum_%s.csv", label)), row.names = FALSE)
  write.csv(data.frame(sample_id = ids, risk_label = risk_label, orient_by = orient_by),
            file.path(PATHS$results, sprintf("02_risk_labels_%s.csv", label)), row.names = FALSE)

  cat(sprintf("MRD %s (orient by %s): n=%d | factors=%d | Cox-sig %d | relapse-sig %d\n",
              label, orient_by, length(ids), ncol(Z),
              sum(cox_tab$FDR < 0.1), sum(relapse_tab$FDR < 0.1)))
  invisible(list(mofa = mofa, Z = Z, risk = risk_label,
                 cox = cox_tab, relapse = relapse_tab))
}

mrd_neg <- run_one("mrdneg", cl$sample_id[cl$mrd_pos == 0], orient_by = "relapse")
mrd_pos <- run_one("mrdpos", cl$sample_id[cl$mrd_pos == 1], orient_by = "os")
