#!/bin/bash
# ===========================================================================
# One-time setup of the `mr_env` conda environment for the MR pipeline.
# Run this ONCE on a login node (it needs internet to reach CRAN + GitHub):
#
#     bash setup_conda_env.sh            # creates env "mr_env"
#     bash setup_conda_env.sh my_env     # custom name
#
# After it finishes, the SLURM script just does `conda activate mr_env`.
# Works with conda or mamba (set CONDA=mamba to use mamba).
# ===========================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="${1:-mr_env}"
CONDA="${CONDA:-conda}"

echo "== Building conda env '${ENV_NAME}' from ${HERE}/environment.yml =="

# 1. Create (or update) the environment from the yml.
if "${CONDA}" env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "-- Env '${ENV_NAME}' exists; updating."
  "${CONDA}" env update -n "${ENV_NAME}" -f "${HERE}/environment.yml" --prune
else
  "${CONDA}" env create -n "${ENV_NAME}" -f "${HERE}/environment.yml"
fi

# 2. Activate it (works inside a non-interactive script).
source "$("${CONDA}" info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

# 3. Install the GitHub-only R packages into this env's library.
#    TwoSampleMR pulls in ieugwasr automatically, but we install it explicitly
#    so the MRCIEU version is used.
echo "== Installing GitHub R packages (TwoSampleMR, ieugwasr, genetics.binaRies, MR-PRESSO, MVMR) =="
Rscript -e '
options(repos = "https://cloud.r-project.org", Ncpus = max(1, parallel::detectCores()))
pkgs <- c("MRCIEU/ieugwasr", "MRCIEU/TwoSampleMR",
          "MRCIEU/genetics.binaRies", "rondolab/MR-PRESSO",
          "WSpiller/MVMR")   # MVMR: conditional F / Q / qhet for mr_mvmr.R
for (p in pkgs) remotes::install_github(p, upgrade = "never")
'

# 4. Verify everything the pipeline needs is importable.
echo "== Verifying =="
Rscript -e '
req <- c("optparse","dplyr","data.table","TwoSampleMR","ieugwasr","MRPRESSO","genetics.binaRies")
ok  <- vapply(req, requireNamespace, logical(1), quietly = TRUE)
for (p in req) cat(sprintf("  [%s] %-16s %s\n", if (ok[[p]]) "ok " else "MISS", p,
                           if (ok[[p]]) as.character(packageVersion(p)) else ""))
if (!all(ok)) { cat("\nSome packages failed to install.\n"); quit(status = 1) }
cat(sprintf("\nPLINK on PATH: %s\n", Sys.which("plink")))
cat("Environment ready. In your SLURM script use: conda activate '"${ENV_NAME}"'\n")
'
