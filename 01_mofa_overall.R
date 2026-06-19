## Overall multi-group MOFA on the full REDIAL cohort.
## Grouping: relapse_category (no / early / intermediate / late).
## Residualized covariates: only the technical/clinical noise we want OFF the
## factors (age, sex, blast%). Subtype, MRD, ancestry, ADI are kept free to
## associate with factors downstream — that's the whole point.
##
## Downstream:
##   - Cox PH on factor scores for OS (time-to-event integration).
##   - Linear association of each factor with ancestry, ADI, subtype, MRD.
##   - High-risk vs low-risk clustering on factor scores (orient by OS).

source("R/utils.R")
redial <- load_redial()
cl <- redial$clinical
ids <- cl$sample_id

adj_cov <- c("age", "sex", "blast")

views <- build_views(redial, ids, residualize_against = adj_cov)
groups <- as.character(cl$relapse_category)

mofa_overall <- fit_mofa(views, groups = groups, num_factors = 10, seed = 1)

Z <- get_factor_scores(mofa_overall)
Z <- Z[ids, , drop = FALSE]  # align row order

## ---- OS time-to-event integration via Cox PH on factors --------------------
cox_tab <- cox_on_factors(
  Z, time = cl$os_time, event = cl$os_event,
  adjust_df = cl[, c("age", "sex")]
)

## Multi-factor risk score via penalized Cox.
if (requireNamespace("glmnet", quietly = TRUE)) {
  cv <- glmnet::cv.glmnet(
    x = Z, y = survival::Surv(cl$os_time, cl$os_event),
    family = "cox", alpha = 0.5
  )
  os_risk_score <- as.numeric(predict(cv, newx = Z, s = "lambda.min"))
} else {
  os_risk_score <- rep(NA_real_, nrow(Z))
}

## ---- Ancestry / ADI association on each factor -----------------------------
anc_assoc <- factor_assoc(
  Z, cl, exposure = "ancestry",
  adjust = c("adi_q", "subtype", "mrd_pos", "age", "sex", "blast")
)
adi_assoc <- factor_assoc(
  Z, cl, exposure = "adi_q",
  adjust = c("ancestry", "subtype", "mrd_pos", "age", "sex", "blast")
)
subtype_assoc <- factor_assoc(
  Z, cl, exposure = "subtype",
  adjust = c("ancestry", "adi_q", "mrd_pos", "age", "sex", "blast")
)
mrd_assoc <- factor_assoc(
  Z, cl, exposure = "mrd_pos",
  adjust = c("ancestry", "adi_q", "subtype", "age", "sex", "blast")
)

## ---- High-risk vs low-risk clustering on factor scores ---------------------
risk_label <- risk_clusters(Z, cl$os_time, cl$os_event)

## Test high-risk enrichment in ancestry x ADI strata.
risk_by_stratum <- cl |>
  mutate(risk_label = risk_label) |>
  count(ancestry, adi_q, risk_label) |>
  group_by(ancestry, adi_q) |>
  mutate(frac_high = n / sum(n)) |>
  filter(risk_label == "high_risk") |>
  ungroup()

## ---- Persist ---------------------------------------------------------------
saveRDS(mofa_overall, file.path(PATHS$results, "mofa_overall.rds"))
write.csv(cox_tab,         file.path(PATHS$results, "01_cox_factors_OS.csv"),  row.names = FALSE)
write.csv(anc_assoc,       file.path(PATHS$results, "01_ancestry_assoc.csv"), row.names = FALSE)
write.csv(adi_assoc,       file.path(PATHS$results, "01_adi_assoc.csv"),      row.names = FALSE)
write.csv(subtype_assoc,   file.path(PATHS$results, "01_subtype_assoc.csv"),  row.names = FALSE)
write.csv(mrd_assoc,       file.path(PATHS$results, "01_mrd_assoc.csv"),      row.names = FALSE)
write.csv(risk_by_stratum, file.path(PATHS$results, "01_risk_by_stratum.csv"), row.names = FALSE)
write.csv(data.frame(sample_id = ids, os_risk_score = os_risk_score, risk_label = risk_label),
          file.path(PATHS$results, "01_per_sample_scores.csv"), row.names = FALSE)

cat("MOFA overall: factors =", ncol(Z),
    "| Cox-significant factors (FDR<0.1):", sum(cox_tab$FDR < 0.1), "\n")
