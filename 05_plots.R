## Visualization layer for the REDIAL pipeline.
## Reads results/*.rds + results/*.csv from scripts 01-04 and writes PDFs to
## results/plots/. Self-contained — re-run after any upstream rerun.

source("R/utils.R")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(survival)
})

PLOTDIR <- file.path(PATHS$results, "plots")
dir.create(PLOTDIR, showWarnings = FALSE, recursive = TRUE)
ggsave_pdf <- function(name, plot, w = 7, h = 5)
  ggsave(file.path(PLOTDIR, paste0(name, ".pdf")), plot, width = w, height = h)

redial <- load_redial()
cl <- redial$clinical
mofa <- readRDS(file.path(PATHS$results, "mofa_overall.rds"))
Z    <- get_factor_scores(mofa)[cl$sample_id, , drop = FALSE]
scores <- read.csv(file.path(PATHS$results, "01_per_sample_scores.csv"))

## ---- (1) MOFA variance explained per modality x factor --------------------
if (!is.null(MOFA2::get_variance_explained(mofa)$r2_per_factor)) {
  ve <- MOFA2::get_variance_explained(mofa)$r2_per_factor
  ve_df <- do.call(rbind, lapply(names(ve), function(g) {
    m <- ve[[g]]
    data.frame(group = g, factor = rownames(m),
               modality = rep(colnames(m), each = nrow(m)),
               r2 = as.vector(m))
  }))
  p <- ggplot(ve_df, aes(factor, modality, fill = r2)) +
    geom_tile() + facet_wrap(~ group) +
    scale_fill_viridis_c() +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "MOFA variance explained per factor x modality", fill = "R²")
  ggsave_pdf("01_mofa_variance_explained", p, w = 10, h = 6)
}

## ---- (2) Forest plot: Cox HRs per factor for OS ---------------------------
cox_tab <- read.csv(file.path(PATHS$results, "01_cox_factors_OS.csv"))
cox_tab$sig <- cox_tab$FDR < 0.1
p <- ggplot(cox_tab, aes(reorder(factor, HR), HR, ymin = lo95, ymax = hi95, color = sig)) +
  geom_pointrange() + geom_hline(yintercept = 1, lty = 2) +
  coord_flip() + scale_y_log10() +
  scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) +
  theme_minimal() +
  labs(x = NULL, y = "HR (OS) per SD of factor", color = "FDR < 0.1",
       title = "Per-factor Cox PH on OS (script 01)")
ggsave_pdf("02_cox_forest_OS", p)

## ---- (3) KM by MOFA-derived risk_label, overall and by ancestry -----------
df_km <- merge(cl, scores[, c("sample_id", "risk_label")], by = "sample_id")
fit_km <- survfit(Surv(os_time, os_event) ~ risk_label, data = df_km)
pdf(file.path(PLOTDIR, "03_km_risk_label.pdf"), width = 7, height = 5)
plot(fit_km, col = c("steelblue", "firebrick"), lwd = 2,
     xlab = "Days", ylab = "OS", main = "MOFA risk label — overall")
legend("topright", c("high_risk", "low_risk"),
       col = c("firebrick", "steelblue"), lty = 1, lwd = 2)
dev.off()

pdf(file.path(PLOTDIR, "04_km_risk_label_by_ancestry.pdf"), width = 10, height = 8)
op <- par(mfrow = c(2, 2))
for (a in sort(unique(df_km$ancestry))) {
  d <- df_km[df_km$ancestry == a, ]
  if (sum(d$os_event) < 2) next
  f <- survfit(Surv(os_time, os_event) ~ risk_label, data = d)
  plot(f, col = c("steelblue", "firebrick"), lwd = 2,
       xlab = "Days", ylab = "OS", main = sprintf("Ancestry = %s (n=%d)", a, nrow(d)))
}
par(op)
dev.off()

## ---- (4) Factor scatter colored by subtype / ancestry / ADI ---------------
Z_df <- data.frame(sample_id = rownames(Z), Z[, 1:2, drop = FALSE]) |>
  merge(cl, by = "sample_id")
fcols <- colnames(Z)[1:2]
for (cv in c("subtype", "ancestry", "adi_q")) {
  p <- ggplot(Z_df, aes(.data[[fcols[1]]], .data[[fcols[2]]], color = factor(.data[[cv]]))) +
    geom_point(alpha = 0.7) + theme_minimal() +
    labs(title = sprintf("Factor 1 vs Factor 2 by %s", cv), color = cv)
  ggsave_pdf(sprintf("05_factor_scatter_%s", cv), p)
}

## ---- (5) Ancestry & ADI factor-association forest -------------------------
for (tag in c("ancestry", "adi")) {
  f <- file.path(PATHS$results, sprintf("01_%s_assoc.csv", tag))
  if (!file.exists(f)) next
  d <- read.csv(f)
  d$sig <- d$p < 0.05
  p <- ggplot(d, aes(reorder(factor, estimate), estimate,
                     ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se,
                     color = sig)) +
    geom_pointrange() + geom_hline(yintercept = 0, lty = 2) +
    coord_flip() + facet_wrap(~ term, scales = "free_x") +
    scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) +
    theme_minimal() +
    labs(x = NULL, y = "Effect on factor score",
         title = sprintf("Factor association with %s", tag))
  ggsave_pdf(sprintf("06_factor_assoc_%s", tag), p, w = 9, h = 6)
}

## ---- (6) High-risk fraction across ancestry x ADI strata ------------------
rb <- read.csv(file.path(PATHS$results, "01_risk_by_stratum.csv"))
p <- ggplot(rb, aes(factor(adi_q), ancestry, fill = frac_high)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.0f%%\n(n=%d)", 100 * frac_high, n)), size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick", limits = c(0, 1)) +
  theme_minimal() +
  labs(x = "ADI quartile", y = "Ancestry",
       title = "High-risk fraction per ancestry x ADI stratum",
       fill = "Frac high-risk")
ggsave_pdf("07_risk_by_stratum_heatmap", p)

## ---- (7) HIMA mediator effects (if 04 ran) --------------------------------
hima_file <- file.path(PATHS$results, "04_hima_survival_all.csv")
if (file.exists(hima_file)) {
  hima <- read.csv(hima_file)
  ec_col <- grep("Mediation|MED|alpha\\*beta", names(hima), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(ec_col)) ec_col <- intersect(c("alpha*beta", "M2Y", "effect"), names(hima))[1]
  if (!is.na(ec_col)) {
    hima$effect <- hima[[ec_col]]
    top <- hima |> group_by(exposure, modality) |>
      slice_max(abs(effect), n = 10) |> ungroup()
    p <- ggplot(top, aes(reorder(ID, effect), effect, fill = modality)) +
      geom_col() + coord_flip() +
      facet_grid(exposure ~ modality, scales = "free", space = "free_y") +
      theme_minimal(base_size = 9) +
      labs(x = NULL, y = "Mediation effect (alpha*beta)",
           title = "Top HIMA-survival mediators by exposure x modality")
    ggsave_pdf("08_hima_top_mediators", p, w = 12, h = 8)
  }
}

cat("Plots written to", PLOTDIR, "\n")
