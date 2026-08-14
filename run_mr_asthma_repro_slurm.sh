#!/bin/bash
#SBATCH --job-name=mr_asthma_repro
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
# Analysis 1:  asthma  -->  all reproductive traits (EXCEPT ANM / POI, already run)
#
# One direction only (asthma is the exposure). Full 4-config sensitivity grid
# (lenient ivs.file / strict re-clump x +/-MHC). Traits, columns, N, build and
# the asthma instrument list all come from the info files (manifests).
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

OUT="${PIPE}/output/asthma_repro/"; mkdir -p "${OUT}"

# All reproductive Short names EXCEPT anm and poi (bweight needs its se.col filled
# in the manifest, otherwise its runs are skipped with a warning).
REPRO_OUTCOMES=endo,aam,pcos,cycle,inf_all,inf_anat,inf_anov,ovar_canc,endo_canc,bweight,ufib,neb

cd "${PIPE}"

Rscript "${SCRIPT}" \
  --exp_info "${RESP_INFO}" --exposures asthma \
  --out_info "${REPR_INFO}" --outcomes "${REPRO_OUTCOMES}" \
  --ld_ref "${LDREF}" \
  --out_prefix "${OUT}" \
  --sensitivity_grid \
  --proxies --proxy_r2 0.8 --proxy_kb 1000 \
  --presso_nbdist 5000
  # add --dry_run to print the resolved traits + grid and exit
  # add --no_presso if MR-PRESSO runtime becomes a problem

echo "Analysis 1 complete. Results in ${OUT}"
