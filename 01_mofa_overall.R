## Overall multi-group MOFA on the full REDIAL cohort.
##
## Role assignments (important):
##   - GROUPING variable (MOFA structure): tumor_subtype
##   - EXPOSURES (tested as predictors of factors): ancestry, ADI quartile
##   - OUTCOMES (tested as functions of factors): MRD status, relapse_category
##     (4-level: no / early / intermediate / late), and OS time-to-event
##   - RESIDUALIZED OUT (technical noise): age, sex, blast%
##
## Downstream:
##   - Cox PH on factor scores for OS (time-to-event integration).
##   - Multinomial logit of relapse_category on factors (4-level outcome).
##   - Nested 2-stage outcome: relapsed y/n, then timing (early/int/late) within
##     relapsers.
##   - Logistic of MRD status on factors.
##   - Ancestry & ADI EXPOSURE effects on each factor.
##   - High-risk vs low-risk clustering on factor scores, oriented by OS.

source("R/utils.R")
redial <- load_redial()
cl <- redial$clinical
ids <- cl$sample_id

adj_cov <- c("age", "sex", "blast")

views  <- build_views(redial, ids, residualize_against = adj_cov)
groups <- as.character(cl$subtype)

mofa_overall <- fit_mofa(views, groups = groups, num_factors = 10, seed = 1)

Z <- get_factor_scores(mofa_overall)[ids, , drop = FALSE]

## ---- OUTCOME: OS time-to-event (Cox PH on each factor) ---------------------
cox_tab <- cox_on_factors(
  Z, time = cl$os_time, event = cl$os_event,
  adjust_df = cl[, c("age", "sex", "subtype")]
)

## Multi-factor OS risk score via penalized Cox.
if (requireNamespace("glmnet", quietly = TRUE)) {
  cv <- glmnet::cv.glmnet(
    x = Z, y = survival::Surv(cl$os_time, cl$os_event),
    family = "cox", alpha = 0.5
  )
  os_risk_score <- as.numeric(predict(cv, newx = Z, s = "lambda.min"))
} else {
  os_risk_score <- rep(NA_real_, nrow(Z))
}

## ---- OUTCOME: relapse_category, 4-level multinomial -----------------------
relapse_multinom <- NULL
if (requireNamespace("nnet", quietly = TRUE)) {
  df_rel <- data.frame(
    y = relevel(cl$relapse_category, ref = "no_relapse"),
    Z, cl[, c("age", "sex", "subtype")]
  )
  fit_mn <- nnet::multinom(y ~ ., data = df_rel, trace = FALSE)
  s <- summary(fit_mn)
  z_stat <- s$coefficients / s$standard.errors
  p <- 2 * (1 - pnorm(abs(z_stat)))
  relapse_multinom <- data.frame(
    timing = rep(rownames(s$coefficients), each = ncol(s$coefficients)),
    term   = rep(colnames(s$coefficients), times = nrow(s$coefficients)),
    estimate = as.vector(t(s$coefficients)),
    se       = as.vector(t(s$standard.errors)),
    p        = as.vector(t(p))
  )
  relapse_multinom <- relapse_multinom[grepl("^Factor", relapse_multinom$term), ]
  relapse_multinom$FDR <- p.adjust(relapse_multinom$p, "BH")
}

## ---- OUTCOME: nested 2-stage relapse --------------------------------------
df_yn <- data.frame(y = cl$relapsed, Z, cl[, c("age", "sex", "subtype")])
relapse_yn_tab <- do.call(rbind, lapply(colnames(Z), function(fac) {
  fit <- glm(as.formula(sprintf("y ~ %s + age + sex + subtype", fac)),
             data = df_yn, family = binomial())
  s <- summary(fit)$coefficients[fac, , drop = FALSE]
  data.frame(factor = fac, OR = exp(s[, "Estimate"]),
             p = s[, "Pr(>|z|)"], row.names = NULL)
}))
relapse_yn_tab$FDR <- p.adjust(relapse_yn_tab$p, "BH")

## Within relapsers, timing as ordinal (early < intermediate < late).
relapsers <- cl$relapsed == 1
relapse_timing_tab <- NULL
if (requireNamespace("MASS", quietly = TRUE) && sum(relapsers) > 30) {
  df_t <- data.frame(
    y = ordered(droplevels(cl$relapse_category[relapsers]),
                levels = c("early", "intermediate", "late")),
    Z[relapsers, ], cl[relapsers, c("age", "sex", "subtype")]
  )
  fit_ord <- MASS::polr(y ~ ., data = df_t, Hess = TRUE)
  s <- summary(fit_ord)$coefficients
  s <- s[grepl("^Factor", rownames(s)), , drop = FALSE]
  relapse_timing_tab <- data.frame(
    factor = rownames(s),
    log_OR = s[, "Value"], se = s[, "Std. Error"],
    p      = 2 * (1 - pnorm(abs(s[, "t value"])))
  )
  relapse_timing_tab$FDR <- p.adjust(relapse_timing_tab$p, "BH")
}

## ---- OUTCOME: MRD status (logistic on factors) -----------------------------
mrd_tab <- do.call(rbind, lapply(colnames(Z), function(fac) {
  df <- data.frame(y = cl$mrd_pos, f = Z[, fac], cl[, c("age", "sex", "subtype")])
  fit <- glm(y ~ ., data = df, family = binomial())
  s <- summary(fit)$coefficients["f", , drop = FALSE]
  data.frame(factor = fac, OR = exp(s[, "Estimate"]),
             p = s[, "Pr(>|z|)"], row.names = NULL)
}))
mrd_tab$FDR <- p.adjust(mrd_tab$p, "BH")

## ---- EXPOSURE: ancestry, ADI effects on each factor ------------------------
anc_assoc <- factor_assoc(
  Z, cl, exposure = "ancestry",
  adjust = c("adi_q", "subtype", "age", "sex", "blast")
)
adi_assoc <- factor_assoc(
  Z, cl, exposure = "adi_q",
  adjust = c("ancestry", "subtype", "age", "sex", "blast")
)

## ---- High-risk vs low-risk clustering on factor scores ---------------------
risk_label <- risk_clusters(Z, cl$os_time, cl$os_event)
risk_by_stratum <- cl |>
  mutate(risk_label = risk_label) |>
  count(ancestry, adi_q, risk_label) |>
  group_by(ancestry, adi_q) |>
  mutate(frac_high = n / sum(n)) |>
  filter(risk_label == "high_risk") |>
  ungroup()

## ---- Persist ---------------------------------------------------------------
saveRDS(mofa_overall, file.path(PATHS$results, "mofa_overall.rds"))
write.csv(cox_tab,           file.path(PATHS$results, "01_cox_factors_OS.csv"),     row.names = FALSE)
write.csv(relapse_yn_tab,    file.path(PATHS$results, "01_relapse_yn.csv"),         row.names = FALSE)
if (!is.null(relapse_multinom))
  write.csv(relapse_multinom, file.path(PATHS$results, "01_relapse_multinom.csv"),  row.names = FALSE)
if (!is.null(relapse_timing_tab))
  write.csv(relapse_timing_tab, file.path(PATHS$results, "01_relapse_timing_ordinal.csv"), row.names = FALSE)
write.csv(mrd_tab,           file.path(PATHS$results, "01_mrd_logit.csv"),          row.names = FALSE)
write.csv(anc_assoc,         file.path(PATHS$results, "01_ancestry_assoc.csv"),     row.names = FALSE)
write.csv(adi_assoc,         file.path(PATHS$results, "01_adi_assoc.csv"),          row.names = FALSE)
write.csv(risk_by_stratum,   file.path(PATHS$results, "01_risk_by_stratum.csv"),    row.names = FALSE)
write.csv(data.frame(sample_id = ids, os_risk_score = os_risk_score, risk_label = risk_label),
          file.path(PATHS$results, "01_per_sample_scores.csv"), row.names = FALSE)

cat("MOFA overall: factors =", ncol(Z),
    "| Cox-sig (FDR<0.1):", sum(cox_tab$FDR < 0.1),
    "| relapse-sig (FDR<0.1):", sum(relapse_yn_tab$FDR < 0.1),
    "| MRD-sig (FDR<0.1):", sum(mrd_tab$FDR < 0.1), "\n")
