#!/bin/bash
#SBATCH --job-name=mr_mvmr
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --partition=icelake        # <-- change to your cluster's partition
#SBATCH --account=MRC-EPID-SL0-CPU # <-- change to your project/account

# ===========================================================================
# Multivariable MR:
#   asthma + smoking     --> POI / ANM
#   asthma + birthweight --> POI / ANM
#
# Each `run` is ONE MVMR model (the two exposures modelled jointly) against the
# outcomes. Instruments are pooled across the two exposures and re-clumped to
# independence (r2<0.01, 1 Mb) before harmonisation. Also emits the univariate
# IVW for each exposure->outcome and a univariable-vs-multivariable forest.
# ===========================================================================
set -uo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mr_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

PIPE=/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline
SCRIPT="${PIPE}/pipeline_scripts/mr_mvmr.R"
LDREF=/rfs/project/rfs-mpB3sSsgAn4/Studies/1000GP_Phase3_genotyped/superpopulation/Phase3_eur
RESP_INFO="${PIPE}/manifests/read_me_respiratory_traits.xlsx"
REPR_INFO="${PIPE}/manifests/read_me_reproductive_traits.xlsx"
OUT="${PIPE}/output/mvmr/"; mkdir -p "${OUT}"

cd "${PIPE}"

run() {  # non-fatal: one model failing won't stop the others
  echo "=== [$(date +%T)] MVMR $* ==="
  Rscript "${SCRIPT}" "$@" \
    --out_info "${REPR_INFO}" --outcomes poi,anm \
    --ld_ref "${LDREF}" --out_prefix "${OUT}" \
    --clump_r2 0.01 --clump_kb 1000 \
    --egger || echo "WARNING: MVMR $* failed; continuing."
}

# asthma is in the respiratory manifest; smoking too; bweight is in the reproductive manifest.
run --exp_info "${RESP_INFO}"               --exposures asthma,smoking
run --exp_info "${RESP_INFO}","${REPR_INFO}" --exposures asthma,bweight

# add --exclude_mhc to either run for the MHC-excluded sensitivity
# add --dry_run to print the resolved model + planned runs and exit
echo "MVMR complete. Results in ${OUT}"
