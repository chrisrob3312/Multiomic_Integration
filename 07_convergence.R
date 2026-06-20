## Convergence between MOFA factors and LUCIDus latent clusters.
##
## Three questions answered:
##
## Q1. Do the MOFA factors driving OS and relapse y/n overlap with the
##     LUCIDus latent omics clusters that mediate the germline / ADI / AMI
##     exposure -> relapse effect? Convergence = same biology, two angles.
##
## Q2. Within the high-risk MOFA cluster (k-means on factor scores, oriented
##     by OS or relapse), which top-loading features are KNOWN B-ALL drivers
##     (sanity check) vs CANDIDATE NOVEL drivers (discovery)?
##
## Q3. Do those candidate features replicate as LUCIDus Z-coefficients on the
##     omics cluster with worst relapse rate? (Cross-tool corroboration.)
##
## Outcomes inspected: OS Cox HR (script 01) AND relapse y/n logistic
## (scripts 01 + 02 MRD-neg).

source("R/utils.R")
suppressPackageStartupMessages({ library(dplyr) })

source_file <- function(label) {
  p <- file.path(PATHS$results, paste0(label, ".rds"))
  if (file.exists(p)) readRDS(p) else NULL
}

mofa_overall <- source_file("mofa_overall")
mofa_mrdneg  <- source_file("mofa_mrdneg")
mofa_mrdpos  <- source_file("mofa_mrdpos")
lucid_snps   <- source_file("lucid_germline_snps")
lucid_adi    <- source_file("lucid_adi_high")
lucid_ami    <- source_file("lucid_ancestry_AMI")
stopifnot("script 01 must run first" = !is.null(mofa_overall))

cl <- load_redial()$clinical

read_csv_safe <- function(name) {
  p <- file.path(PATHS$results, name)
  if (file.exists(p)) read.csv(p) else NULL
}

## ---- Q1: which MOFA factors drive each outcome ----------------------------
cox_overall      <- read_csv_safe("01_cox_factors_OS.csv")
relapse_overall  <- read_csv_safe("01_relapse_yn.csv")
relapse_mrdneg   <- read_csv_safe("02_relapse_logit_mrdneg.csv")

top_factors_per_outcome <- list(
  OS_overall          = if (!is.null(cox_overall))     head(cox_overall[order(cox_overall$p), "factor"], 5),
  relapse_overall     = if (!is.null(relapse_overall)) head(relapse_overall[order(relapse_overall$p), "factor"], 5),
  relapse_mrdneg      = if (!is.null(relapse_mrdneg))  head(relapse_mrdneg[order(relapse_mrdneg$p), "factor"], 5)
)
print(top_factors_per_outcome)

## ---- Q1 continued: MOFA factor scores vs LUCIDus cluster posteriors -------
## LUCIDus exposes posterior cluster probabilities per sample. In "parallel"
## mode there is one cluster vector per omics block; we collapse to one
## "any-block worst" assignment (per modality, then aggregate).
extract_lucid_posterior <- function(fit) {
  if (is.null(fit)) return(NULL)
  pp <- tryCatch(fit$post.p, error = function(e) NULL)
  if (is.null(pp)) pp <- tryCatch(fit$res_Estimates$post.p, error = function(e) NULL)
  if (is.null(pp)) {
    message("could not locate posterior probabilities on fit object — skipping")
    return(NULL)
  }
  ## pp is typically a list per modality of (samples x K).
  if (is.list(pp)) {
    do.call(cbind, lapply(seq_along(pp), function(i) {
      m <- pp[[i]]
      colnames(m) <- sprintf("blk%d_K%d", i, seq_len(ncol(m)))
      m
    }))
  } else pp
}

mofa_Z <- get_factor_scores(mofa_overall)
mofa_Z <- mofa_Z[intersect(rownames(mofa_Z), cl$sample_id), , drop = FALSE]

convergence_rows <- list()
for (tag in c("snps", "adi", "ami")) {
  fit <- get(paste0("lucid_", tag))
  pp  <- extract_lucid_posterior(fit)
  if (is.null(pp)) next
  pp  <- pp[intersect(rownames(pp), rownames(mofa_Z)), , drop = FALSE]
  Z2  <- mofa_Z[rownames(pp), , drop = FALSE]
  cors <- cor(Z2, pp, use = "pairwise.complete.obs")
  for (i in seq_len(nrow(cors))) for (j in seq_len(ncol(cors)))
    convergence_rows[[length(convergence_rows) + 1]] <- data.frame(
      lucid_exposure = tag,
      mofa_factor    = rownames(cors)[i],
      lucid_cluster  = colnames(cors)[j],
      pearson_r      = cors[i, j]
    )
}
convergence_tab <- do.call(rbind, convergence_rows)
if (!is.null(convergence_tab)) {
  convergence_tab$abs_r <- abs(convergence_tab$pearson_r)
  write.csv(convergence_tab[order(-convergence_tab$abs_r), ],
            file.path(PATHS$results, "07_mofa_lucid_convergence.csv"),
            row.names = FALSE)
}

## ---- Q2: known vs candidate drivers in high-risk MOFA factors -------------
BALL_DRIVERS <- c(
  "IKZF1", "PAX5", "EBF1", "ETV6", "RUNX1", "TCF3", "PBX1", "HLF",
  "MEF2D", "ZNF384", "NUTM1", "DUX4", "ERG",
  "CDKN2A", "CDKN2B", "RB1", "TP53", "BTG1",
  "BCR", "ABL1", "ABL2", "JAK1", "JAK2", "JAK3", "CRLF2", "EPOR",
  "CSF1R", "PDGFRB", "PDGFRA", "FLT3", "NTRK3", "IL7R", "SH2B3",
  "KRAS", "NRAS", "PTPN11", "NF1", "BRAF", "PIK3CA", "PIK3R1",
  "KMT2A", "CREBBP", "SETD2", "EZH2", "WHSC1", "NSD2", "ARID1A",
  "EED", "SUZ12", "SF3B1", "SRSF2", "U2AF1", "DDX3X",
  "TBL1XR1", "FBXW7", "WT1", "PTEN"
)

loadings_file <- file.path(PATHS$results, "06_top_loadings.csv")
if (file.exists(loadings_file)) {
  loadings <- read.csv(loadings_file)
  driver_tag <- loadings |>
    mutate(known_BALL_driver = feature %in% BALL_DRIVERS,
           classification = ifelse(known_BALL_driver,
                                   "known_BALL_driver",
                                   "candidate_novel")) |>
    arrange(factor, modality, desc(abs_loading))

  ## Restrict to top factors driving each outcome.
  high_factors <- unique(unlist(top_factors_per_outcome))
  high_factors <- high_factors[!is.na(high_factors)]
  driver_in_high <- driver_tag[driver_tag$factor %in% high_factors, ]

  write.csv(driver_in_high,
            file.path(PATHS$results, "07_drivers_in_high_risk_factors.csv"),
            row.names = FALSE)

  ## Summary per (factor, modality, classification).
  summ <- driver_in_high |>
    group_by(factor, modality, classification) |>
    summarise(n = dplyr::n(),
              top5 = paste(head(feature[order(-abs_loading)], 5), collapse = ", "),
              .groups = "drop")
  write.csv(summ,
            file.path(PATHS$results, "07_driver_summary_per_factor.csv"),
            row.names = FALSE)

  cat(sprintf("Driver tag: %d known, %d candidate-novel features across %d high-risk-driving factors\n",
              sum(driver_in_high$known_BALL_driver),
              sum(!driver_in_high$known_BALL_driver),
              length(unique(driver_in_high$factor))))
} else {
  message("06_top_loadings.csv not found — run 06_annotation.R first to enable Q2/Q3")
}

## ---- Q3: do the candidate features also load on the high-risk LUCIDus
##         cluster? (cross-tool corroboration) -------------------------------
## LUCIDus exposes per-omics-block Z-coefficients per cluster. Pull the
## features with the largest |coefficient| in the cluster that has the
## strongest relapse association, and intersect with the MOFA candidate list.
extract_lucid_Z_coef <- function(fit, top_n = 50) {
  if (is.null(fit)) return(NULL)
  coefs <- tryCatch(fit$res_Estimates$Z, error = function(e) NULL)
  if (is.null(coefs)) coefs <- tryCatch(fit$pars$mu, error = function(e) NULL)
  if (is.null(coefs)) return(NULL)
  do.call(rbind, lapply(seq_along(coefs), function(i) {
    m <- coefs[[i]]; if (is.null(m)) return(NULL)
    df <- data.frame(modality_index = i, feature = rownames(m),
                     max_abs_coef = apply(m, 1, function(x) max(abs(x), na.rm = TRUE)))
    df <- df[order(-df$max_abs_coef), ]
    head(df, top_n)
  }))
}

lucid_top <- extract_lucid_Z_coef(lucid_snps)
if (!is.null(lucid_top) && file.exists(loadings_file)) {
  loadings <- read.csv(loadings_file)
  cross <- merge(loadings, lucid_top, by = "feature")
  cross$known_BALL_driver <- cross$feature %in% BALL_DRIVERS
  cross <- cross[order(-cross$abs_loading), ]
  write.csv(cross,
            file.path(PATHS$results, "07_mofa_lucid_feature_overlap.csv"),
            row.names = FALSE)
  cat(sprintf("Feature overlap between top MOFA loadings and top LUCIDus Z-coefs: %d features (%d known drivers)\n",
              nrow(cross), sum(cross$known_BALL_driver)))
}
