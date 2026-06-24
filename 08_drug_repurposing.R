## Drug repurposing on top high-risk-driving MOFA factors.
##
## Two layers, both open-source:
##   Target-level    : DGIdb (REST) + Open Targets (GraphQL) per candidate gene
##                     from 07_drivers_in_high_risk_factors.csv. Merged into a
##                     per-target druggability + cancer-evidence table.
##   Signature-level : signatureSearch (Bioconductor) queries LINCS L1000 with
##                     the RNA factor-loading vector to rank perturbagens whose
##                     expression profile REVERSES the high-risk signature.
##
## Both layers degrade gracefully — missing packages or unreachable APIs print
## a warning and skip rather than erroring out.

source("R/utils.R")
suppressPackageStartupMessages({ library(dplyr) })

DGIDB_URL <- "https://dgidb.org/api/v2/interactions.json"
OT_URL    <- "https://api.platform.opentargets.org/api/v4/graphql"
EFO_BALL  <- "EFO_0000220"   # B-cell acute lymphoblastic leukemia

DRUGDIR <- file.path(PATHS$results, "drugs")
dir.create(DRUGDIR, showWarnings = FALSE, recursive = TRUE)

## ---- Pull candidate genes from 07 -----------------------------------------
drv_file <- file.path(PATHS$results, "07_drivers_in_high_risk_factors.csv")
if (!file.exists(drv_file)) stop("Run 07_convergence.R first.")
drv <- read.csv(drv_file)
gene_drv <- drv |> filter(modality %in% c("rna", "cnv"))
candidates <- unique(gene_drv$feature)
cat("Querying drug evidence for", length(candidates), "candidate genes\n")

## ---- (A) DGIdb -------------------------------------------------------------
query_dgidb <- function(genes) {
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    message("[DGIdb] install httr + jsonlite to enable"); return(NULL)
  }
  ## DGIdb v2 REST: comma-delimited genes, batch <= ~200 to be safe.
  out <- list()
  chunks <- split(genes, ceiling(seq_along(genes) / 100))
  for (chunk in chunks) {
    r <- tryCatch(httr::GET(DGIDB_URL, query = list(genes = paste(chunk, collapse = ","))),
                  error = function(e) { message("[DGIdb] ", e$message); NULL })
    if (is.null(r) || httr::status_code(r) != 200) next
    j <- jsonlite::fromJSON(httr::content(r, "text", encoding = "UTF-8"), flatten = TRUE)
    if (length(j$matchedTerms) == 0) next
    df <- j$matchedTerms
    if (!"interactions" %in% names(df)) next
    rows <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
      ix <- df$interactions[[i]]
      if (is.null(ix) || nrow(ix) == 0) return(NULL)
      data.frame(gene = df$searchTerm[i],
                 drug = ix$drugName,
                 interaction_types = sapply(ix$interactionTypes, paste, collapse = "|"),
                 score = ix$score)
    }))
    if (!is.null(rows)) out[[length(out) + 1]] <- rows
  }
  if (length(out)) do.call(rbind, out) else NULL
}

dgi <- query_dgidb(candidates)
if (!is.null(dgi)) {
  write.csv(dgi, file.path(DRUGDIR, "08_dgidb_interactions.csv"), row.names = FALSE)
  cat("DGIdb: ", nrow(dgi), "drug-gene interactions across",
      length(unique(dgi$gene)), "genes\n")
}

## ---- (B) Open Targets ------------------------------------------------------
query_opentargets <- function(genes, disease_efo = EFO_BALL) {
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    message("[OpenTargets] install httr + jsonlite to enable"); return(NULL)
  }
  ## Resolve gene symbols to Ensembl IDs first, then pull associatedDiseases
  ## scores and knownDrugs.
  ens_q <- '
    query ($q: String!) {
      search(queryString: $q, entityNames: ["target"]) {
        hits { id name entity }
      }
    }'
  drug_q <- '
    query ($ens: String!) {
      target(ensemblId: $ens) {
        approvedSymbol
        knownDrugs { count rows { drug { id name } phase mechanismOfAction } }
        associatedDiseases(efoIds: ["%s"]) { rows { disease { name id } score } }
      }
    }'
  rows <- list()
  for (g in genes) {
    e <- tryCatch(
      httr::POST(OT_URL,
                 body = list(query = ens_q, variables = list(q = g)),
                 encode = "json"),
      error = function(e) NULL)
    if (is.null(e) || httr::status_code(e) != 200) next
    j <- jsonlite::fromJSON(httr::content(e, "text", encoding = "UTF-8"), flatten = TRUE)
    hits <- j$data$search$hits
    if (is.null(hits) || nrow(hits) == 0) next
    ens <- hits$id[1]
    d <- tryCatch(
      httr::POST(OT_URL,
                 body = list(query = sprintf(drug_q, disease_efo),
                             variables = list(ens = ens)),
                 encode = "json"),
      error = function(e) NULL)
    if (is.null(d) || httr::status_code(d) != 200) next
    jd <- jsonlite::fromJSON(httr::content(d, "text", encoding = "UTF-8"), flatten = TRUE)
    tg <- jd$data$target
    if (is.null(tg)) next
    kd <- tg$knownDrugs$rows
    if (!is.null(kd) && nrow(kd) > 0) {
      rows[[length(rows) + 1]] <- data.frame(
        gene = g, ensembl = ens,
        drug = kd$drug.name, phase = kd$phase,
        mechanism = kd$mechanismOfAction,
        BALL_assoc_score = if (length(tg$associatedDiseases$rows))
          tg$associatedDiseases$rows$score[1] else NA
      )
    }
    Sys.sleep(0.1)  # courtesy rate limit
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}

ot <- query_opentargets(candidates, disease_efo = EFO_BALL)
if (!is.null(ot)) {
  write.csv(ot, file.path(DRUGDIR, "08_opentargets_knowndrugs.csv"), row.names = FALSE)
  cat("Open Targets:", nrow(ot), "known-drug rows across",
      length(unique(ot$gene)), "genes\n")
}

## ---- Per-target merged druggability table ---------------------------------
merge_targets <- function(dgi, ot) {
  parts <- list()
  if (!is.null(dgi)) {
    parts$dgi <- dgi |> group_by(gene) |>
      summarise(dgidb_n_drugs = dplyr::n_distinct(drug),
                dgidb_top_drugs = paste(head(unique(drug), 5), collapse = ", "),
                .groups = "drop")
  }
  if (!is.null(ot)) {
    parts$ot <- ot |> group_by(gene) |>
      summarise(ot_n_known_drugs = dplyr::n_distinct(drug),
                ot_max_phase = suppressWarnings(max(phase, na.rm = TRUE)),
                ot_top_drugs = paste(head(unique(drug), 5), collapse = ", "),
                ot_BALL_assoc = suppressWarnings(max(BALL_assoc_score, na.rm = TRUE)),
                .groups = "drop")
  }
  Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), parts)
}
merged <- merge_targets(dgi, ot)
if (!is.null(merged)) {
  merged <- merge(merged,
                  unique(gene_drv[, c("feature", "factor", "modality", "classification")]),
                  by.x = "gene", by.y = "feature", all.x = TRUE)
  write.csv(merged, file.path(DRUGDIR, "08_target_druggability.csv"), row.names = FALSE)
}

## ---- (C) Signature-level: signatureSearch on LINCS L1000 ------------------
## Requires the signatureSearchData reference (~large download, one-time):
##   BiocManager::install(c("signatureSearch", "signatureSearchData"))
sig_search <- function() {
  if (!requireNamespace("signatureSearch", quietly = TRUE)) {
    message("[signatureSearch] not installed — skip"); return(invisible(NULL))
  }
  mofa <- readRDS(file.path(PATHS$results, "mofa_overall.rds"))
  W <- MOFA2::get_weights(mofa, as.data.frame = FALSE)
  if (is.null(W$rna)) { message("[signatureSearch] no RNA view"); return(NULL) }

  cox  <- read.csv(file.path(PATHS$results, "01_cox_factors_OS.csv"))
  rel  <- read.csv(file.path(PATHS$results, "01_relapse_yn.csv"))
  ranked_factors <- unique(c(
    head(cox[order(cox$p), "factor"], 3),
    head(rel[order(rel$p), "factor"], 3)
  ))

  for (fac in ranked_factors) {
    v <- W$rna[, fac]
    v <- v[!is.na(v) & !is.na(names(v))]
    if (length(v) < 200) next
    up   <- names(sort(v, decreasing = TRUE))[1:150]
    down <- names(sort(v, decreasing = FALSE))[1:150]
    res <- tryCatch({
      qs <- signatureSearch::qSig(query = list(upset = up, downset = down),
                                  gess_method = "LINCS",
                                  refdb = "lincs")
      signatureSearch::gess_lincs(qs)
    }, error = function(e) { message("[signatureSearch ", fac, "] ", e$message); NULL })
    if (is.null(res)) next
    df <- as.data.frame(signatureSearch::result(res))
    df$factor <- fac
    write.csv(df,
              file.path(DRUGDIR, sprintf("08_lincs_query_%s.csv", fac)),
              row.names = FALSE)
  }
}
sig_search()

cat("Drug repurposing outputs written to", DRUGDIR, "\n")
