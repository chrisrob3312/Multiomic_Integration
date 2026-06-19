## Biology annotation of MOFA factor loadings.
##
## For each factor, extract the top-loading features per modality and:
##   - methylation : map CpG ID -> nearest gene (Illumina EPIC/450k annotation)
##   - rna         : run GSEA against MSigDB Hallmark / Reactome (fgsea + msigdbr)
##   - metabolome  : look up KEGG pathway membership (if metabolite IDs map)
##   - cnv         : map region/gene name -> cytoband and B-ALL driver list
##
## On toy data (feature IDs "f1, f2, ..."), each block falls through with a
## friendly skip message. On REAL REDIAL data, swap in your feature ID
## namespace (cg-IDs for methylation, HGNC/Ensembl for RNA, HMDB/KEGG for
## metabolome, HGNC for CNV-by-gene).

source("R/utils.R")

mofa <- readRDS(file.path(PATHS$results, "mofa_overall.rds"))
TOPN <- 50  # number of top loadings per modality per factor to annotate

## ---- Pull top loadings per modality x factor ------------------------------
get_loadings <- function(mofa) {
  W <- MOFA2::get_weights(mofa, as.data.frame = FALSE)
  do.call(rbind, lapply(names(W), function(mod) {
    mat <- W[[mod]]
    do.call(rbind, lapply(colnames(mat), function(fac) {
      v <- mat[, fac]
      ord <- order(-abs(v))[seq_len(min(TOPN, length(v)))]
      data.frame(modality = mod, factor = fac,
                 feature = rownames(mat)[ord],
                 loading = v[ord], abs_loading = abs(v[ord]))
    }))
  }))
}
loadings <- get_loadings(mofa)
write.csv(loadings, file.path(PATHS$results, "06_top_loadings.csv"), row.names = FALSE)

## ---- Methylation: CpG -> gene mapping -------------------------------------
annot_methylation <- function(df) {
  if (!any(grepl("^cg", df$feature))) {
    message("[methylation] feature IDs don't look like CpG IDs — skipping annotation")
    return(NULL)
  }
  pkg <- "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"
  if (!requireNamespace(pkg, quietly = TRUE)) {
    pkg <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("[methylation] install IlluminaHumanMethylationEPICanno (or 450k) to annotate")
      return(NULL)
    }
  }
  ann <- minfi::getAnnotation(get(pkg))
  keep <- intersect(df$feature, rownames(ann))
  if (!length(keep)) return(NULL)
  out <- merge(df, as.data.frame(ann[keep, c("chr", "pos", "UCSC_RefGene_Name",
                                              "UCSC_RefGene_Group", "Relation_to_Island")]),
               by.x = "feature", by.y = "row.names", all.x = TRUE)
  write.csv(out, file.path(PATHS$results, "06_methylation_annot.csv"), row.names = FALSE)
  out
}

## ---- RNA: GSEA on factor loadings -----------------------------------------
annot_rna_gsea <- function(mofa) {
  if (!requireNamespace("fgsea", quietly = TRUE) ||
      !requireNamespace("msigdbr", quietly = TRUE)) {
    message("[rna] install fgsea + msigdbr to run GSEA on RNA factor loadings")
    return(NULL)
  }
  W <- MOFA2::get_weights(mofa, as.data.frame = FALSE)
  if (is.null(W$rna)) return(NULL)
  Wrna <- W$rna
  if (!any(grepl("^[A-Z][A-Z0-9]+$", rownames(Wrna)))) {
    message("[rna] feature IDs don't look like HGNC symbols — GSEA may be empty")
  }
  hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  pathways <- split(hallmark$gene_symbol, hallmark$gs_name)
  res <- do.call(rbind, lapply(colnames(Wrna), function(fac) {
    ranks <- Wrna[, fac]; ranks <- ranks[!is.na(ranks)]
    if (length(ranks) < 50) return(NULL)
    g <- fgsea::fgsea(pathways, ranks, minSize = 10, maxSize = 500)
    g$factor <- fac
    as.data.frame(g[, c("factor", "pathway", "NES", "pval", "padj", "size")])
  }))
  if (!is.null(res))
    write.csv(res, file.path(PATHS$results, "06_rna_gsea_hallmark.csv"), row.names = FALSE)
  res
}

## ---- Metabolome: KEGG pathway lookup (skeleton) ---------------------------
annot_metabolome <- function(df) {
  if (!any(grepl("^HMDB|^C[0-9]{5}", df$feature))) {
    message("[metabolome] feature IDs don't look like HMDB/KEGG — skipping")
    return(NULL)
  }
  ## Real implementation: query KEGGREST or MetaboAnalystR.
  ## Toy fall-through: just persist top loadings tagged.
  out <- df
  write.csv(out, file.path(PATHS$results, "06_metabolome_top.csv"), row.names = FALSE)
  out
}

## ---- CNV: cytoband + B-ALL driver overlay ---------------------------------
BALL_DRIVERS <- c("IKZF1", "PAX5", "CDKN2A", "CDKN2B", "ETV6", "RUNX1",
                  "TP53", "TCF3", "PBX1", "MLL", "KMT2A", "BCR", "ABL1",
                  "JAK2", "CRLF2", "IL7R", "CREBBP", "EBF1", "RB1", "BTG1")
annot_cnv <- function(df) {
  df$is_BALL_driver <- df$feature %in% BALL_DRIVERS
  write.csv(df, file.path(PATHS$results, "06_cnv_top_with_drivers.csv"), row.names = FALSE)
  df
}

## ---- Run per modality -----------------------------------------------------
res <- list(
  methylation = annot_methylation(subset(loadings, modality == "methylation")),
  rna         = annot_rna_gsea(mofa),
  metabolome  = annot_metabolome(subset(loadings, modality == "metabolome")),
  cnv         = annot_cnv(subset(loadings, modality == "cnv"))
)

## ---- Cross-modality per-factor biology summary ----------------------------
summary_per_factor <- loadings |>
  dplyr::group_by(factor, modality) |>
  dplyr::summarise(
    n_top         = dplyr::n(),
    max_abs_load  = max(abs_loading, na.rm = TRUE),
    top_features  = paste(head(feature[order(-abs_loading)], 5), collapse = ", "),
    .groups = "drop"
  )
write.csv(summary_per_factor,
          file.path(PATHS$results, "06_factor_biology_summary.csv"),
          row.names = FALSE)

cat("Annotation written to", PATHS$results, "\n")
