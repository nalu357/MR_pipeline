#!/bin/bash
#SBATCH --job-name=mr_pipeline
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1                  # single R process (NOT multi-node/MPI)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6          # data.table uses these threads
#SBATCH --mem=64G                  # PLINK proxy reserves ~30G (--proxy_mem) + R holds the sumstats
#SBATCH --partition=icelake        # <-- change to your cluster's partition name
#SBATCH --account=MRC-EPID-SL0-CPU # <-- change to your project/account code

set -euo pipefail

# ---- Environment ------------------------------------------------------------
# Activate the conda env built by setup_conda_env.sh (has R + all packages +
# MR-PRESSO + PLINK). If you use a different name, change mr_env below.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mr_env

# Let data.table use the allocated cores.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

PIPE=/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline
LDREF=/rfs/project/rfs-mpB3sSsgAn4/Studies/1000GP_Phase3_genotyped/superpopulation/Phase3_eur
OUT=${PIPE}/output/
ASTHMA=/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/autoimmune_traits/data/asthma_hg19_lifted_with_n.tsv.gz
ANM=/rfs/project/rfs-ic9ZbslBEP4/Datasets/Projects/ANM_2025/g2g/anm_ofh_ukb_23andMe_metal/anm_ofh_ukb_23andMe_metal_gwas_funcannoQC.txt
ANM_IVS=/rfs/project/rfs-ic9ZbslBEP4/Datasets/Projects/ANM_2025/g2g/anm_ofh_ukb_23andMe_metal/p01_signal_selection/anm_ofh_ukb_23andMe_metal_all_independent_signals_GWASinfo.txt

cd "${PIPE}"

################################################################################
# ----------------------- Less strict clump with MHC ---------------------------
################################################################################
# ANM --> asthma  (pre-defined ANM signals as instruments; full ANM GWAS for proxies)
Rscript "${PIPE}/pipeline_scripts/mr_pipeline.R" \
  --exp_gwas "${ANM}" \
  --exp_name ANM_ofh_meta \
  --exp_snp SNP --exp_chr CHR --exp_pos BP \
  --exp_ea ALLELE1 --exp_nea ALLELE0 --exp_beta BETA --exp_se SE \
  --exp_p P_BOLT_LMM --exp_eaf A1FREQ --exp_n_total 679746 \
  --exp_ivs "${ANM_IVS}" \
  --out_gwas "${ASTHMA}" \
  --out_name asthma \
  --out_snp rsid --out_chr '#CHR' --out_pos POS --out_n n --out_ncase N_case \
  --out_ea ALT --out_nea REF --out_beta inv_var_meta_beta --out_se inv_var_meta_sebeta \
  --out_p inv_var_meta_p --out_eaf all_meta_AF \
  --ld_ref "${LDREF}" \
  --out_prefix "${OUT}" \
  --proxies --proxy_r2 0.8 --proxy_kb 1000 \
  --presso_nbdist 5000

# Asthma --> ANM
Rscript "${PIPE}/pipeline_scripts/mr_pipeline.R" \
  --exp_gwas "${ASTHMA}" \
  --exp_name asthma \
  --exp_snp rsid --exp_chr '#CHR' --exp_pos POS \
  --exp_ea ALT --exp_nea REF --exp_beta inv_var_meta_beta --exp_se inv_var_meta_sebeta \
  --exp_p inv_var_meta_p --exp_eaf all_meta_AF \
  --exp_n n --exp_ncase N_case \
  --out_gwas "${ANM}" \
  --out_name ANM_ofh_meta \
  --out_snp SNP --out_chr CHR --out_pos BP \
  --out_ea ALLELE1 --out_nea ALLELE0 --out_beta BETA --out_se SE \
  --out_p P_BOLT_LMM_INF --out_eaf A1FREQ --out_n_total 679746 --out_type continuous \
  --ld_ref "${LDREF}" \
  --out_prefix "${OUT}" \
  --clump_r2 0.2 --clump_kb 10000 \
  --proxies --proxy_r2 0.8 --proxy_kb 1000 \
  --presso_nbdist 10000

echo "All MR runs complete."
