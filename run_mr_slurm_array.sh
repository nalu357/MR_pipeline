#!/bin/bash
# Submit independent GWAS--GWAS MR commands concurrently:
#   sbatch --array=1-$(wc -l < mr_commands.txt)%12 run_mr_slurm_array.sh mr_commands.txt
#
# Each non-empty line of the manifest is one fully quoted `Rscript ...` command.
# All array tasks must share the same --cache_dir on cluster storage.  The %12
# cap is deliberate: tune it to the available memory / PLINK reference I/O.
# The manifest is trusted shell input owned by the analyst; do not use an
# untrusted manifest.
#SBATCH --job-name=mr_array
#SBATCH --output=mr_array_%A_%a.out
#SBATCH --error=mr_array_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --partition=icelake
#SBATCH --account=MRC-EPID-SL0-CPU

set -euo pipefail

MANIFEST=${1:?Usage: sbatch --array=1-N run_mr_slurm_array.sh /path/to/mr_commands.txt}
TASK=${SLURM_ARRAY_TASK_ID:?This script must be submitted as a Slurm array.}
COMMAND=$(sed -n "${TASK}p" "${MANIFEST}")
if [[ -z "${COMMAND}" || "${COMMAND}" =~ ^[[:space:]]*# ]]; then
  echo "No command on manifest line ${TASK}; exiting."
  exit 0
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mr_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

echo "[$(date -Is)] array task ${TASK}: ${COMMAND}"
bash -lc "${COMMAND}"
