#!/bin/bash
#SBATCH --job-name=mr_resp_menopause
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --partition=icelake        # <-- change to your cluster's partition
#SBATCH --account=MRC-EPID-SL0-CPU # <-- change to your project/account

# ===========================================================================
# Analysis 2:  respiratory conditions  -->  POI / ANM
#
# One direction only (respiratory traits are exposures). Every respiratory trait
# in the manifest is used as an exposure (IBF is auto-skipped: Downloaded=no);
# outcomes are restricted to poi and anm. Full 4-config sensitivity grid.
# Each respiratory trait uses its own ivs.file for the lenient config.
# ===========================================================================
set -uo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mr_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# ---- Paths (SET THESE) ------------------------------------------------------
PIPE=/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline
SCRIPT="${PIPE}/pipeline_scripts/mr_grid.R"
LDREF=/rfs/project/rfs-mpB3sSsgAn4/Studies/1000GP_Phase3_genotyped/superpopulation/Phase3_eur

# Info files (.xlsx or .csv). Copy the manifests somewhere the job can read them.
RESP_INFO="${PIPE}/manifests/read_me_respiratory_traits.xlsx"
REPR_INFO="${PIPE}/manifests/read_me_reproductive_traits.xlsx"

OUT="${PIPE}/output/respiratory_menopause/"; mkdir -p "${OUT}"

cd "${PIPE}"

Rscript "${SCRIPT}" \
  --exp_info "${RESP_INFO}" \
  --out_info "${REPR_INFO}" --outcomes poi,anm \
  --ld_ref "${LDREF}" \
  --out_prefix "${OUT}" \
  --sensitivity_grid \
  --proxies --proxy_r2 0.8 --proxy_kb 1000 \
  --presso_nbdist 5000
  # add --exposures asthma,COPD,... to restrict which respiratory traits run
  # add --dry_run to print the resolved traits + grid and exit
  # add --no_presso if MR-PRESSO runtime becomes a problem

echo "Analysis 2 complete. Results in ${OUT}"
