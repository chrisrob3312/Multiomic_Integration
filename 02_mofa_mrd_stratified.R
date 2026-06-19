## MRD-stratified MOFA runs: one for MRD-negative, one for MRD-positive.
## Grouping: tumor_subtype (so subtype heterogeneity is absorbed by the
## multi-group framework instead of forcing nested per-subtype refits).
## Residualized covariates: only technical/clinical noise (age, sex, blast%).
## Relapse_category, ancestry, ADI are kept free for downstream association.
##
## This is the standalone MRD-neg high-risk-vs-low-risk subtyping arm: fit MOFA
## on the whole MRD-neg pool (group=subtype), derive factor scores, cluster
## into high/low risk oriented by OS, then test ancestry/ADI enrichment in the
## high-risk cluster within each stratum.

source("R/utils.R")
redial <- load_redial()
cl <- redial$clinical

adj_cov <- c("age", "sex", "blast")

run_one <- function(label, ids) {
  views  <- build_views(redial, ids, residualize_against = adj_cov)
  groups <- as.character(cl[ids, "subtype"])
  mofa <- fit_mofa(views, groups = groups, num_factors = 8, seed = 1)

  Z <- get_factor_scores(mofa)[ids, , drop = FALSE]
  sub_cl <- cl[ids, , drop = FALSE]

  cox_tab <- cox_on_factors(
    Z, time = sub_cl$os_time, event = sub_cl$os_event,
    adjust_df = sub_cl[, c("age", "sex")]
  )
  anc_assoc <- factor_assoc(
    Z, sub_cl, exposure = "ancestry",
    adjust = c("adi_q", "subtype", "age", "sex", "blast")
  )
  adi_assoc <- factor_assoc(
    Z, sub_cl, exposure = "adi_q",
    adjust = c("ancestry", "subtype", "age", "sex", "blast")
  )

  risk_label <- risk_clusters(Z, sub_cl$os_time, sub_cl$os_event)
  risk_by_stratum <- sub_cl |>
    mutate(risk_label = risk_label) |>
    count(ancestry, adi_q, risk_label) |>
    group_by(ancestry, adi_q) |>
    mutate(frac_high = n / sum(n)) |>
    filter(risk_label == "high_risk") |>
    ungroup()

  saveRDS(mofa,         file.path(PATHS$results, sprintf("mofa_%s.rds", label)))
  write.csv(cox_tab,         file.path(PATHS$results, sprintf("02_cox_%s.csv", label)),         row.names = FALSE)
  write.csv(anc_assoc,       file.path(PATHS$results, sprintf("02_ancestry_%s.csv", label)),    row.names = FALSE)
  write.csv(adi_assoc,       file.path(PATHS$results, sprintf("02_adi_%s.csv", label)),         row.names = FALSE)
  write.csv(risk_by_stratum, file.path(PATHS$results, sprintf("02_risk_by_stratum_%s.csv", label)), row.names = FALSE)

  cat(sprintf("MRD %s: n=%d, factors=%d, Cox-sig (FDR<0.1)=%d\n",
              label, length(ids), ncol(Z), sum(cox_tab$FDR < 0.1)))
  invisible(list(mofa = mofa, Z = Z, risk = risk_label, cox = cox_tab))
}

mrd_neg <- run_one("mrdneg", cl$sample_id[cl$mrd_pos == 0])
mrd_pos <- run_one("mrdpos", cl$sample_id[cl$mrd_pos == 1])
