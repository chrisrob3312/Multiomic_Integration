#!/usr/bin/env bash
## Wrapper for setup/install_packages.R — meant for an interactive HPC node
## inside screen/tmux. Logs to setup/install.log and exits non-zero on any
## package failure so you can grep the log on reattach.
##
##   screen -S redial-install
##   bash setup/install.sh
##   # ctrl-a d to detach; screen -r redial-install to reattach
##
## Optional env vars (see setup/install_packages.R for full list):
##   REDIAL_SKIP_LINCS=1            skip the ~7-15 GB LINCS reference download
##   REDIAL_USER_LIB=/path/to/lib   install into a user-writable library
##   REDIAL_NCPUS=8                 cap parallel CRAN install workers

set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="${SCRIPT_DIR}/install.log"

command -v Rscript >/dev/null || {
  echo "Rscript not on PATH — load your R module first (e.g., module load R/4.4)."
  exit 127
}

echo "Logging to $LOG"
echo "Started: $(date)" | tee -a "$LOG"
Rscript "${SCRIPT_DIR}/install_packages.R" 2>&1 | tee -a "$LOG"
rc=${PIPESTATUS[0]}
echo "Finished: $(date) (exit ${rc})" | tee -a "$LOG"
exit $rc
