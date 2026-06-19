## Generate a small synthetic REDIAL-shaped dataset for pipeline dry-runs.
## NOT a substitute for real data — feature counts and noise scales are toy.
## Run once: source("data/synthetic_redial.R")

set.seed(42)

N <- 240
P <- list(methylation = 500, rna = 800, metabolome = 150, cnv = 300, germline = 40)

ancestry <- sample(c("EUR", "AMR", "AFR", "AMI"),
                   N, replace = TRUE, prob = c(0.35, 0.40, 0.10, 0.15))
adi_q    <- sample(1:4, N, replace = TRUE)
subtype  <- sample(c("ETV6_RUNX1", "HighHyperdiploid", "BCR_ABL1_like",
                     "DUX4", "Bother"),
                  N, replace = TRUE, prob = c(0.25, 0.25, 0.15, 0.15, 0.20))
mrd_pos  <- rbinom(N, 1, 0.30)
age      <- round(runif(N, 1, 18), 1)
sex      <- sample(c("F", "M"), N, replace = TRUE)
blast    <- round(runif(N, 60, 99), 1)

## Relapse category — biased by MRD and ancestry to make the toy run nontrivial.
risk_lp <- 0.6 * mrd_pos +
           0.4 * (ancestry == "AMI") +
           0.3 * (ancestry == "AMR") +
           0.3 * (adi_q >= 3) +
           0.2 * (subtype == "BCR_ABL1_like")
p_relapse <- plogis(risk_lp - 0.5)
relapsed  <- rbinom(N, 1, p_relapse)
relapse_category <- ifelse(relapsed == 0, "no_relapse",
                    sample(c("early", "intermediate", "late"),
                           N, replace = TRUE, prob = c(0.4, 0.35, 0.25)))
relapse_category[relapsed == 0] <- "no_relapse"

## Overall survival — exponential with hazard tied to risk_lp.
hazard <- 0.05 * exp(0.8 * risk_lp)
os_time  <- round(pmin(rexp(N, hazard), 8) * 365)  # cap at 8y, in days
os_event <- as.integer(rexp(N, hazard) < 8)

clinical <- data.frame(
  sample_id = sprintf("S%04d", seq_len(N)),
  ancestry, adi_q, subtype, mrd_pos, age, sex, blast,
  relapse_category = factor(relapse_category,
                            levels = c("no_relapse", "early", "intermediate", "late")),
  relapsed, os_time, os_event,
  stringsAsFactors = FALSE
)
rownames(clinical) <- clinical$sample_id

## Omics — each modality gets a baseline + a risk-correlated signal so MOFA
## actually recovers something on the toy data.
gen_view <- function(n_feat, lp, sd_noise = 1, signal_frac = 0.05) {
  signal_n <- max(1, round(signal_frac * n_feat))
  base <- matrix(rnorm(N * n_feat, 0, sd_noise), N, n_feat)
  base[, seq_len(signal_n)] <- base[, seq_len(signal_n)] + outer(lp, rnorm(signal_n, 0.6, 0.2))
  colnames(base) <- paste0("f", seq_len(n_feat))
  rownames(base) <- clinical$sample_id
  base
}

omics <- list(
  methylation = gen_view(P$methylation, risk_lp, sd_noise = 1.0, signal_frac = 0.05),
  rna         = gen_view(P$rna,         risk_lp, sd_noise = 1.0, signal_frac = 0.05),
  metabolome  = gen_view(P$metabolome,  risk_lp, sd_noise = 1.0, signal_frac = 0.08),
  cnv         = gen_view(P$cnv,         risk_lp, sd_noise = 1.0, signal_frac = 0.03)
)

## Inject block missingness: 15% of samples missing metabolome, 5% missing CNV.
miss_metab <- sample(seq_len(N), round(0.15 * N))
omics$metabolome[miss_metab, ] <- NA
miss_cnv <- sample(seq_len(N), round(0.05 * N))
omics$cnv[miss_cnv, ] <- NA

## Germline carrier matrix — used as exposure in LUCIDus and HIMA, not as MOFA view.
## REAL DATA: replace with your imputed/called SNP genotypes (0/1/2 dosage or
## 0/1 carrier) over your curated B-ALL risk-variant list (IKZF1, ARID5B,
## GATA3, CEBPE, PIP4K2A, BMI1-PIP4K2A locus, etc.).
germline <- matrix(rbinom(N * P$germline, 1, 0.10), N, P$germline)
colnames(germline) <- paste0("rs", seq_len(P$germline))
rownames(germline) <- clinical$sample_id

## Make first 10 SNPs the "risk-variant panel" — boost their carrier frequency
## in high-lp samples so the toy LUCIDus / HIMA runs find signal.
risk_panel <- paste0("rs", 1:10)
for (s in risk_panel) {
  bump <- rbinom(N, 1, plogis(risk_lp - 1.0))
  germline[, s] <- pmax(germline[, s], bump)
}

germline_burden <- rowSums(germline[, risk_panel])  # count of risk-panel variants
## Polygenic risk score: linear combo of all 40 SNPs with random effect sizes
## biased so the first 10 carry positive risk-aligned weights.
beta_prs <- rnorm(P$germline, 0, 0.2)
beta_prs[1:10] <- abs(beta_prs[1:10]) + 0.3
germline_prs <- as.numeric(germline %*% beta_prs)

redial <- list(
  clinical = clinical,
  omics    = omics,  # methylation, rna, metabolome, cnv (CNV is from RNA-seq, separate MOFA view)
  germline = list(
    matrix     = germline,        # N x 40 carrier matrix (all SNPs)
    risk_panel = risk_panel,      # column names of curated risk-variant SNPs (top 10)
    burden     = germline_burden, # sum over risk_panel
    prs        = germline_prs     # weighted PRS over all SNPs
  )
)

dir.create("data", showWarnings = FALSE)
saveRDS(redial, "data/synthetic_redial.rds")
cat("Wrote data/synthetic_redial.rds:",
    sprintf("N=%d, modalities=%s\n", N, paste(names(omics), collapse = ",")))
