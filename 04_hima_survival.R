## High-dimensional mediation for OS time-to-event, multi-omic.
##
## HIMA::hima_survival: exposure -> high-dim mediator features -> Cox(OS).
## FDR-controlled per-feature mediation effects, complementary to LUCIDus'
## latent-cluster quasi-mediation.
##
## Three exposures (each run independently):
##   (a) germline risk-variant burden        (continuous)
##   (b) ADI high-deprivation (q >= 3)       (binary)
##   (c) Amerindigenous ancestry indicator   (binary)
##
## Run separately per modality so the per-feature mediator list is
## interpretable within each omics layer. Complete-case within each modality;
## HIMA does not natively handle block-missing mediators.

source("R/utils.R")
if (!requireNamespace("HIMA", quietly = TRUE)) stop("install.packages('HIMA')")

redial <- load_redial()
cl <- redial$clinical

## HIMA-survival takes a SCALAR exposure X. For germline, use the PRS (single
## column) — for a multi-SNP entry, see LUCIDus (03), which is the natural
## home for jointly-modeled SNP exposures.
exposures <- list(
  germline_prs = redial$germline$prs,
  adi_high     = as.integer(cl$adi_q >= 3),
  ancestry_AMI = as.integer(cl$ancestry == "AMI")
)

modalities <- redial$omics  # methylation / rna / metabolome / cnv

cov_xm <- model.matrix(~ age + sex + subtype, data = cl)[, -1, drop = FALSE]
cov_my <- cov_xm  # same baseline adjustments on both arms

run_hima_one <- function(exposure_name, x_vec, mod_name, M) {
  keep <- complete.cases(M) & !is.na(x_vec)
  if (sum(keep) < 30) {
    message(sprintf("skip %s x %s: n=%d after complete-case",
                    exposure_name, mod_name, sum(keep)))
    return(NULL)
  }
  fit <- tryCatch(
    HIMA::hima_survival(
      X      = x_vec[keep],
      M      = M[keep, , drop = FALSE],
      OT     = cl$os_time[keep],
      status = cl$os_event[keep],
      COV.XM = cov_xm[keep, , drop = FALSE],
      COV.MY = cov_my[keep, , drop = FALSE],
      verbose = FALSE
    ),
    error = function(e) { message("HIMA failed: ", e$message); NULL }
  )
  if (is.null(fit) || nrow(fit) == 0) return(NULL)
  fit$exposure <- exposure_name
  fit$modality <- mod_name
  fit
}

all_results <- list()
for (en in names(exposures)) {
  for (mn in names(modalities)) {
    tag <- sprintf("%s__%s", en, mn)
    cat("HIMA-survival:", tag, "\n")
    out <- run_hima_one(en, exposures[[en]], mn, modalities[[mn]])
    if (!is.null(out)) all_results[[tag]] <- out
  }
}

if (length(all_results) > 0) {
  combined <- do.call(rbind, lapply(names(all_results), function(k) {
    df <- all_results[[k]]
    df$tag <- k
    df
  }))
  write.csv(combined,
            file.path(PATHS$results, "04_hima_survival_all.csv"),
            row.names = FALSE)
  cat(sprintf("HIMA-survival: %d significant mediators across %d exposure x modality pairs\n",
              nrow(combined), length(all_results)))
} else {
  cat("HIMA-survival: no significant mediators across any exposure x modality pair\n")
}

saveRDS(all_results, file.path(PATHS$results, "04_hima_survival_raw.rds"))
