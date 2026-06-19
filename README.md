# REDIAL multi-omics integration

Two-tool stack for the REDIAL pediatric B-ALL cohort: **MOFA2** for unsupervised
subtyping and survival-aware factor discovery, **LUCIDus** for exposure → omics
→ outcome quasi-mediation. **mixOmics DIABLO** is included as an optional
supervised confirmatory script.

## Pipeline

```
                      ┌───────────────────────────────────────────┐
                      │  data/synthetic_redial.R (toy data)       │
                      │  → R/utils.R (load, QC, residualize)      │
                      └───────────────────┬───────────────────────┘
                                          │
       ┌──────────────────────────────────┼──────────────────────────────────┐
       ▼                                  ▼                                  ▼
 01_mofa_overall.R               02_mofa_mrd_stratified.R          03_lucidus_mediation.R
 group = relapse_category        MRD-neg / MRD-pos runs            germline burden / ADI /
 covariates: subtype, MRD,       group = tumor_subtype             ancestry → omics → relapse
   ancestry, ADI, age, sex,      covariates: relapse, ancestry,
   blast%                          ADI, age, sex, blast%
       │                                  │                                  │
       └────► Cox PH on factors for OS ◄──┘                                  │
       └────► Factor regression on ancestry × ADI                            │
       └────► High-risk vs low-risk clustering on factor scores              │
                                                                             │
                                  99_diablo_optional.R                       │
                                  supervised confirmatory                    │
                                  (run only if MOFA result is weak)          │
```

## Modalities

- methylation array (M-values, gaussian)
- bulk RNA-seq (vst-transformed counts, gaussian)
- tumor metabolome (log-scaled, gaussian)
- RNA-derived CNV (log2 ratios, gaussian; binary calls optional)
- germline risk-variant carrier status (bernoulli, used as **exposure** in LUCIDus, not as MOFA view)

## Why this stack

- **MOFA2** handles block-missing modalities natively, mixed likelihoods, and
  small-n strata via a probabilistic Bayesian factor model with ARD priors.
  The multi-group framework absorbs relapse-category / subtype structure so
  nested stratification is minimized.
- **LUCIDus 3.x** is the only mature R package for the quasi-mediation framing
  we need (latent omics clusters between exposure and outcome) with bootstrap
  inference and g-computation.
- **DIABLO** is supervised — kept as an add-on, not the primary engine, because
  it overfits below n ≈ 50/class.

## Run order

```r
source("data/synthetic_redial.R")   # one-time toy data generation
source("R/utils.R")                  # auto-sourced by all 0x_ scripts
source("01_mofa_overall.R")
source("02_mofa_mrd_stratified.R")
source("03_lucidus_mediation.R")
# optional:
source("99_diablo_optional.R")
```

## Package install (one-time)

```r
install.packages(c("BiocManager", "survival", "glmnet", "ggplot2",
                   "dplyr", "tidyr", "LUCIDus", "mixOmics"))
BiocManager::install(c("MOFA2"))
```
