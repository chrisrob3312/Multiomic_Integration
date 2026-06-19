# REDIAL multi-omics integration

Two-tool primary stack plus two mediation arms for the REDIAL pediatric B-ALL
cohort.

## Variable roles (lock these in)

| Role | Variables |
|---|---|
| **Grouping** (MOFA structure) | `subtype` (overall and MRD-stratified runs) |
| **Exposures** (predictors of factors / mediation upstream) | `ancestry` (Graf), `adi_q` (quartile), germline risk-variant burden |
| **Outcomes** (functions of factors / mediation downstream) | `mrd_pos`, `relapse_category` (no / early / intermediate / late), `Surv(os_time, os_event)` |
| **Residualized noise** (regressed out pre-MOFA) | `age`, `sex`, `blast%` |

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
 group = subtype                 MRD-neg / MRD-pos runs            quasi-mediation:
 outcomes:                       group = subtype                   exposure → omics
   - OS Cox PH                   outcomes:                            cluster → relapse
   - relapse 4-level multinom      - OS Cox PH                     exposures: germline,
   - nested relapsed→timing        - relapse 4-level                  ADI≥3, AMI ancestry
   - MRD logit                     - MRD comparisons               
 exposures:                                                        
   - ancestry, ADI on factors                                      04_hima_survival.R
                                                                   feature-level mediation
                                                                   exposure → omics
                                                                      features → OS Cox
                                                                   FDR-controlled per modality
                                                                  
                                  99_diablo_optional.R
                                  supervised confirmatory
                                  Y = MOFA risk label or
                                      relapse_category
                                  (run only if MOFA is weak)
```

## Tools and why

| Script | Tool | Why |
|---|---|---|
| 01, 02 | **MOFA2** | Probabilistic factor model; handles block-missing modalities, mixed likelihoods, small-n strata; multi-group framework absorbs subtype heterogeneity so we avoid nested per-subtype refits. |
| 03 | **LUCIDus 3.x** | Quasi-mediation with latent omics clusters between exposure and outcome; bootstrap inference; g-computation for causal effects. |
| 04 | **HIMA** | Feature-level high-dimensional mediation with FDR; `hima_survival` for OS Cox outcome. Complements LUCIDus by naming *which* features mediate, where LUCIDus names *which clusters*. |
| 99 | **mixOmics DIABLO** | Optional supervised confirmatory — runs only if MOFA risk-class separation is weak. |

## Modalities

- methylation array (M-values, gaussian)
- bulk RNA-seq (vst-transformed counts, gaussian)
- tumor metabolome (log-scaled, gaussian)
- RNA-derived CNV (log2 ratios, gaussian)
- germline risk-variant carrier status (used as **exposure**, not as MOFA view)

## Run order

```r
source("data/synthetic_redial.R")   # one-time toy data generation
source("R/utils.R")                  # auto-sourced by all 0x_ scripts
source("01_mofa_overall.R")
source("02_mofa_mrd_stratified.R")
source("03_lucidus_mediation.R")
source("04_hima_survival.R")
source("05_plots.R")                 # forest plots, KM, factor scatters, heatmaps
source("06_annotation.R")            # CpG→gene, RNA GSEA, CNV→B-ALL drivers
# optional:
source("99_diablo_optional.R")
```

## Outputs

- `results/*.csv` — per-script tidy tables (Cox, multinomial, ordinal, logistic, association, HIMA mediators, loading summaries)
- `results/*.rds` — MOFA fits, LUCIDus fits, raw HIMA outputs (for re-plotting)
- `results/plots/*.pdf` — forest plots, KM curves (overall + by ancestry), factor scatters, ancestry × ADI risk heatmap, HIMA top-mediator bars
- `results/06_*` — annotated loadings (methylation CpG→gene, RNA GSEA Hallmark, CNV B-ALL driver overlay)

## Package install (one-time)

```r
install.packages(c("BiocManager", "survival", "glmnet", "ggplot2",
                   "dplyr", "tidyr", "nnet", "MASS",
                   "LUCIDus", "mixOmics", "HIMA",
                   "msigdbr"))
BiocManager::install(c("MOFA2", "fgsea", "minfi",
                       "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"))
```
