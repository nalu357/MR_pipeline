#!/bin/bash
#SBATCH --job-name=mr_grid
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
# Manifest-driven bidirectional MR grid.
#
# One Rscript call replaces the hand-written per-trait column/N/IVs arrays: the
# exposure and outcome traits (and all their columns, sample sizes, build and
# instrument lists) are read from "info files". See the README section
# "Trait info files (manifests)" for the expected columns.
#
#   --exp_info / --out_info : one or more info CSVs (comma-separated) defining
#                             the exposure set and the outcome set.
#   --exposures / --outcomes: optional comma lists to restrict to a subset of
#                             Short names (omit = all traits in the file).
#   --bidirectional         : also run outcome-as-exposure for every pair.
#   --sensitivity_grid      : run all 4 configs (lenient/strict x +/-MHC),
#                             writing '', 'rigid_', 'noMHC_', 'noMHC_rigid_'.
#   --data_dir              : base dir for relative file.name / ivs.file entries.
#
# Tip: add --dry_run to print the resolved traits + the exact grid that would
# run (and nothing else) before committing a long job.
# ===========================================================================
set -uo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mr_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# ---- Paths (SET THESE) ------------------------------------------------------
PIPE=/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline
SCRIPT="${PIPE}/pipeline_scripts/mr_grid.R"
LDREF=/rfs/project/rfs-mpB3sSsgAn4/Studies/1000GP_Phase3_genotyped/superpopulation/Phase3_eur
DATA=/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/autoimmune_traits/data   # base for file.name / ivs.file

# Info files. Here the same manifest supplies both sides; use two files if your
# exposure and outcome traits live in different tables.
EXP_INFO="${PIPE}/read_me_autoimmune_traits.csv"
OUT_INFO="${PIPE}/read_me_menopause_traits.csv"

OUT="${PIPE}/output/autoimmune/"; mkdir -p "${OUT}"

cd "${PIPE}"

Rscript "${SCRIPT}" \
  --exp_info "${EXP_INFO}" \
  --out_info "${OUT_INFO}" \
  --data_dir "${DATA}" \
  --bidirectional \
  --sensitivity_grid \
  --ld_ref "${LDREF}" \
  --out_prefix "${OUT}" \
  --proxies --proxy_r2 0.8 --proxy_kb 1000 \
  --presso_nbdist 5000
  # --exposures asthma,EO      # optional: restrict the exposure set
  # --outcomes  POI,ANM        # optional: restrict the outcome set
  # --dry_run                  # print the plan and exit

echo "Grid complete. Results (+ plots-ready *_full_mr_results.csv) in ${OUT}"
