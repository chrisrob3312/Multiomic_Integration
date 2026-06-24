## One-shot package install for the REDIAL multi-omics pipeline.
## Designed to run on an interactive HPC node inside screen/tmux:
##
##   screen -S redial-install
##   Rscript setup/install_packages.R 2>&1 | tee setup/install.log
##   # ctrl-a d to detach; screen -r redial-install to reattach
##
## Idempotent — already-installed packages are skipped. Failures are logged
## but do not stop subsequent installs. Total wall time on a fresh R is
## ~30-90 min depending on Bioconductor build cache; signatureSearchData adds
## a ~7-15 GB LINCS HDF5 download on top.
##
## Flags (set as env vars before running):
##   REDIAL_SKIP_LINCS=1     skip signatureSearchData (huge; do later if you
##                           don't need drug signature search yet)
##   REDIAL_USER_LIB=/path   override library install path (defaults to
##                           the user's personal library)
##   REDIAL_NCPUS=N          parallel CRAN install workers (default detect)

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = max(1800, getOption("timeout", 60)))   # for large downloads

ncpus <- as.integer(Sys.getenv("REDIAL_NCPUS",
                               max(1, parallel::detectCores() - 1)))
options(Ncpus = ncpus)

user_lib <- Sys.getenv("REDIAL_USER_LIB", "")
if (nzchar(user_lib)) {
  dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
  .libPaths(c(user_lib, .libPaths()))
}

skip_lincs <- nzchar(Sys.getenv("REDIAL_SKIP_LINCS"))

cat("R:        ", R.version.string, "\n")
cat("libPaths: ", paste(.libPaths(), collapse = "; "), "\n")
cat("Ncpus:    ", ncpus, "\n")
cat("Skip LINCS reference: ", skip_lincs, "\n\n")

CRAN_PKGS <- c(
  "BiocManager",
  "survival", "glmnet",
  "ggplot2", "dplyr", "tidyr",
  "nnet", "MASS",
  "httr", "jsonlite",
  "msigdbr",
  "LUCIDus", "mixOmics", "HIMA"
)

BIOC_PKGS <- c(
  "MOFA2",
  "fgsea",
  "minfi",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "signatureSearch"
)
LINCS_DATA_PKG <- "signatureSearchData"

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

safe_install <- function(pkg, installer) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[%s] OK   %s (already installed)\n", ts(), pkg))
    return(invisible(TRUE))
  }
  t0 <- Sys.time()
  cat(sprintf("[%s] INST %s ...\n", ts(), pkg))
  ok <- tryCatch({ installer(pkg); TRUE },
                 error = function(e) {
                   cat(sprintf("[%s] FAIL %s: %s\n", ts(), pkg, conditionMessage(e)))
                   FALSE
                 })
  if (ok) {
    cat(sprintf("[%s] DONE %s (%.1fs)\n",
                ts(), pkg, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  invisible(ok)
}

## ---- CRAN ------------------------------------------------------------------
cat("\n=== CRAN packages ===\n")
cran_install <- function(p) install.packages(p, dependencies = TRUE)
for (p in CRAN_PKGS) safe_install(p, cran_install)

## ---- BiocManager bootstrap & Bioconductor ---------------------------------
cat("\n=== Bioconductor packages ===\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  stop("BiocManager install failed — cannot continue with Bioconductor packages.")
}
bioc_install <- function(p) BiocManager::install(p, update = FALSE, ask = FALSE)
for (p in BIOC_PKGS) safe_install(p, bioc_install)

if (skip_lincs) {
  cat(sprintf("\n[%s] SKIP %s (REDIAL_SKIP_LINCS set)\n", ts(), LINCS_DATA_PKG))
} else {
  cat("\n=== signatureSearchData (LINCS reference — large, slow) ===\n")
  safe_install(LINCS_DATA_PKG, bioc_install)
}

## ---- MOFA2 Python backend warm-up (basilisk) -------------------------------
cat("\n=== MOFA2 python backend (basilisk) ===\n")
if (requireNamespace("MOFA2", quietly = TRUE)) {
  ok <- tryCatch({
    invisible(MOFA2::get_default_training_options(
      MOFA2::create_mofa(list(view1 = matrix(rnorm(20), 5, 4)))))
    cat(sprintf("[%s] DONE MOFA2 basilisk python env\n", ts()))
    TRUE
  }, error = function(e) {
    cat(sprintf("[%s] WARN MOFA2 python warm-up failed: %s\n", ts(), conditionMessage(e)))
    cat("     Will be retried on first MOFA fit.\n"); FALSE
  })
}

## ---- Summary ---------------------------------------------------------------
cat("\n=== Summary ===\n")
all_pkgs <- c(CRAN_PKGS, BIOC_PKGS, if (!skip_lincs) LINCS_DATA_PKG)
status <- sapply(all_pkgs, requireNamespace, quietly = TRUE)
cat(sprintf("%-50s %s\n", "package", "status"))
for (i in seq_along(all_pkgs))
  cat(sprintf("%-50s %s\n", all_pkgs[i], if (status[i]) "OK" else "MISSING"))

cat(sprintf("\n[%s] DONE. %d/%d packages installed.\n",
            ts(), sum(status), length(all_pkgs)))
if (any(!status)) {
  cat("Missing packages — rerun this script or install manually:\n")
  cat(paste("  -", all_pkgs[!status], collapse = "\n"), "\n")
  quit(save = "no", status = 1)
}
quit(save = "no", status = 0)
