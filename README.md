# REDIAL multi-omics integration

Two-tool primary stack plus two mediation arms for the REDIAL pediatric B-ALL
cohort.

## Variable roles (lock these in)

| Role | Variables |
|---|---|
| **Grouping** (MOFA structure) | `subtype` (overall and MRD-stratified runs) |
| **Exposures (for LUCIDus / HIMA — NOT in MOFA)** | `ancestry` (Graf), `adi_q` (quartile), **germline**: multi-SNP matrix over the curated risk-variant panel (primary), PRS, or burden count |
| **Outcomes** (functions of factors / mediation downstream) | `mrd_pos`, `relapse_category` (no / early / intermediate / late), `Surv(os_time, os_event)` |
| **Residualized noise** (regressed out pre-MOFA) | `age`, `sex`, `blast%` |

### Germline encoding for LUCIDus G

LUCIDus' `G` matrix accepts continuous OR categorical OR mixed columns. The
"categorical" element in LUCIDus is the latent omics cluster `K`, not the
exposure — so SNPs can enter as a multi-column dosage/carrier matrix.
Three encodings are provided in `03_lucidus_mediation.R`:

- `mode = "snps"` — N × 10 multi-column risk-variant carrier matrix (**primary**)
- `mode = "prs"` — single weighted polygenic risk score column
- `mode = "burden"` — single count of risk-panel carrier variants

HIMA-survival expects a scalar `X` per run, so the germline arm in `04` uses
the PRS column.

### MOFA does not cluster

MOFA produces *factor scores* per sample. Sample clustering (high-risk vs
low-risk label) is a downstream **k-means on factor scores**, oriented by the
clinical outcome of the run:

- Script 01 (overall): clusters oriented by **OS**.
- Script 02 (MRD-negative): clusters oriented by **relapse y/n** — the story
  we care about is "which factor-loaded pathways drive relapse despite MRD
  negativity."
- Script 02 (MRD-positive): clusters oriented by **OS**.

### Modalities (kept as separate MOFA views)

- methylation array (gaussian M-values)
- bulk RNA-seq (gaussian, vst-transformed)
- tumor metabolome (gaussian, log-scaled)
- **CNV from RNA-seq** — kept as its own MOFA view, distinct from the raw RNA
  view, so the integration can ask whether copy-number-driven and
  expression-driven signal are jointly or independently associated with
  outcomes.

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
source("07_convergence.R")           # MOFA factors ↔ LUCIDus clusters + driver tagging
source("08_drug_repurposing.R")      # DGIdb + Open Targets + LINCS L1000 signature query
# optional:
source("99_diablo_optional.R")
```

## Pediatric B-ALL molecular subtypes (canonical, used in `subtype` column)

| Risk tier | Subtype | Token used in `redial$clinical$subtype` |
|---|---|---|
| Favorable | ETV6::RUNX1 (t(12;21)) | `ETV6_RUNX1` |
| Favorable | High hyperdiploid (>50 chr) | `HighHyperdiploid` |
| Favorable | DUX4-rearranged / ERG-deregulated | `DUX4r` |
| Intermediate | TCF3::PBX1 (t(1;19)) | `TCF3_PBX1` |
| Intermediate | MEF2D-rearranged | `MEF2Dr` |
| Intermediate | ZNF384-rearranged | `ZNF384r` |
| Intermediate | NUTM1-rearranged | `NUTM1r` |
| Intermediate | PAX5alt (P80R, other) | `PAX5alt` |
| Intermediate | IKZF1 N159Y | `IKZF1_N159Y` |
| Intermediate | ETV6::RUNX1-like | `ETV6_RUNX1_like` |
| Adverse | BCR::ABL1 (Ph+, t(9;22)) | `BCR_ABL1` |
| Adverse | Ph-like, CRLF2-r | `Ph_like_CRLF2` |
| Adverse | Ph-like, JAK2 / EPOR | `Ph_like_JAK2` |
| Adverse | Ph-like, ABL-class | `Ph_like_ABLclass` |
| Adverse | KMT2A-rearranged (formerly MLL) | `KMT2Ar` |
| Adverse | iAMP21 | `iAMP21` |
| Very high | Low hypodiploid / near-haploid | `LowHypodiploid` |
| Residual | B-other / NOS | `Bother` |

Use these tokens in `redial$clinical$subtype` for consistent grouping. Real
data: call subtypes from RNA-seq using ALLCatchR or ALLSorts, then merge.

## Driver discovery / sanity check

`07_convergence.R` answers three questions:

1. **Convergence**: do the MOFA factors driving OS (Cox PH) and relapse y/n
   (logistic) correlate with the LUCIDus latent omics clusters mediating
   exposure → relapse? Same biology, two angles → confidence boost.
2. **Sanity check**: in the high-risk-driving MOFA factors, which top-loading
   features are on the canonical pediatric B-ALL driver list (IKZF1, PAX5,
   CDKN2A/B, KMT2A, JAK2, CRLF2, EBF1, BCR/ABL1, …)? Recovery of known
   drivers in the toy / real data confirms the pipeline is working.
3. **Candidate novel drivers**: top-loading features NOT on the canonical
   list → flagged as `candidate_novel`. Outputs include per-(factor, modality)
   tables of known vs candidate, ranked by absolute loading.

Cross-tool corroboration: candidate features that ALSO show up as
high-coefficient LUCIDus Z-coefs in the relapse-associated cluster are the
strongest novel-driver candidates — `07_mofa_lucid_feature_overlap.csv`.

To take a candidate forward you'd want orthogonal evidence: variant calls
(GATK/Strelka2 on RNA-seq or WGS), fusion calls (Arriba, STAR-Fusion), CNV
breakpoints from paired DNA, and ideally functional follow-up.

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
                   "msigdbr", "httr", "jsonlite"))
BiocManager::install(c("MOFA2", "fgsea", "minfi",
                       "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                       "signatureSearch", "signatureSearchData"))
```
