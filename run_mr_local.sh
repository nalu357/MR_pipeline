#!/bin/bash
# ===========================================================================
# Run the full 4-config MR grid WITHOUT SLURM (e.g. in an RStudio terminal).
#
#   bash run_mr_local.sh
#
# For a long run that survives the terminal closing, background it:
#   nohup bash run_mr_local.sh > run_mr_local.log 2>&1 &
#   tail -f run_mr_local.log
# ===========================================================================
set -uo pipefail   # no -e, so one failing config doesn't abort the whole grid

# ---- Environment ------------------------------------------------------------
# If your RStudio session's own R already has all the packages (TwoSampleMR,
# ieugwasr, MRPRESSO, R.utils, genetics.binaRies), you can COMMENT OUT the two
# conda lines below. Otherwise activate the conda env built by setup_conda_env.sh.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mr_env
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"   # data.table threads

# ---- Paths ------------------------------------------------------------------
PIPE=/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline
SCRIPT="${PIPE}/pipeline_scripts/mr_pipeline.R"
LDREF=/rfs/project/rfs-mpB3sSsgAn4/Studies/1000GP_Phase3_genotyped/superpopulation/Phase3_eur
OUT="${PIPE}/output/"
ASTHMA=/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/autoimmune_traits/data/asthma_hg19_lifted_with_n.tsv.gz
ANM=/rfs/project/rfs-ic9ZbslBEP4/Datasets/Projects/ANM_2025/g2g/anm_ofh_ukb_23andMe_metal/anm_ofh_ukb_23andMe_metal_gwas_funcannoQC.txt
ANM_IVS=/rfs/project/rfs-ic9ZbslBEP4/Datasets/Projects/ANM_2025/g2g/anm_ofh_ukb_23andMe_metal/p01_signal_selection/anm_ofh_ukb_23andMe_metal_all_independent_signals_GWASinfo.txt

# Common column mappings (arrays keep each run readable; '#CHR' passes safely)
ANM_EXP=(--exp_snp SNP --exp_chr CHR --exp_pos BP
         --exp_ea ALLELE1 --exp_nea ALLELE0 --exp_beta BETA --exp_se SE
         --exp_p P_BOLT_LMM --exp_eaf A1FREQ --exp_n_total 679746)
ANM_OUT=(--out_snp SNP --out_chr CHR --out_pos BP
         --out_ea ALLELE1 --out_nea ALLELE0 --out_beta BETA --out_se SE
         --out_p P_BOLT_LMM_INF --out_eaf A1FREQ --out_n_total 679746 --out_type continuous)
ASTHMA_EXP=(--exp_snp rsid --exp_chr '#CHR' --exp_pos POS
            --exp_ea ALT --exp_nea REF --exp_beta inv_var_meta_beta --exp_se inv_var_meta_sebeta
            --exp_p inv_var_meta_p --exp_eaf all_meta_AF --exp_n n --exp_ncase N_case)
ASTHMA_OUT=(--out_snp rsid --out_chr '#CHR' --out_pos POS --out_n n --out_ncase N_case
            --out_ea ALT --out_nea REF --out_beta inv_var_meta_beta --out_se inv_var_meta_sebeta
            --out_p inv_var_meta_p --out_eaf all_meta_AF)
PROXY=(--proxies --proxy_r2 0.8 --proxy_kb 1000)

cd "${PIPE}"

################################################################################
# ----------------------- Less strict clump, with MHC  (output/) ---------------
################################################################################
# ANM --> asthma  (pre-defined ANM signals as instruments; full ANM GWAS for proxies)
Rscript "${SCRIPT}" --exp_gwas "${ANM}" --exp_name ANM_ofh_meta "${ANM_EXP[@]}" \
  --exp_ivs "${ANM_IVS}" \
  --out_gwas "${ASTHMA}" --out_name asthma "${ASTHMA_OUT[@]}" \
  --ld_ref "${LDREF}" --out_prefix "${OUT}" "${PROXY[@]}" --presso_nbdist 5000 \
  || echo "WARNING: [lenient+MHC] ANM->asthma failed; continuing."

# Asthma --> ANM
Rscript "${SCRIPT}" --exp_gwas "${ASTHMA}" --exp_name asthma "${ASTHMA_EXP[@]}" \
  --out_gwas "${ANM}" --out_name ANM_ofh_meta "${ANM_OUT[@]}" \
  --ld_ref "${LDREF}" --out_prefix "${OUT}" --clump_r2 0.2 --clump_kb 10000 "${PROXY[@]}" --presso_nbdist 10000 \
  || echo "WARNING: [lenient+MHC] asthma->ANM failed; continuing."

################################################################################
# ----------------------- More strict clump, with MHC  (output/rigid_) ---------
################################################################################
# ANM --> asthma  (full ANM GWAS re-clumped at r2<0.001)
Rscript "${SCRIPT}" --exp_gwas "${ANM}" --exp_name ANM_ofh_meta "${ANM_EXP[@]}" \
  --out_gwas "${ASTHMA}" --out_name asthma "${ASTHMA_OUT[@]}" \
  --ld_ref "${LDREF}" --out_prefix "${OUT}rigid_" --clump_r2 0.001 --clump_kb 10000 "${PROXY[@]}" --presso_nbdist 5000 \
  || echo "WARNING: [strict+MHC] ANM->asthma failed; continuing."

# Asthma --> ANM
Rscript "${SCRIPT}" --exp_gwas "${ASTHMA}" --exp_name asthma "${ASTHMA_EXP[@]}" \
  --out_gwas "${ANM}" --out_name ANM_ofh_meta "${ANM_OUT[@]}" \
  --ld_ref "${LDREF}" --out_prefix "${OUT}rigid_" --clump_r2 0.001 --clump_kb 10000 "${PROXY[@]}" --presso_nbdist 10000 \
  || echo "WARNING: [strict+MHC] asthma->ANM failed; continuing."

################################################################################
# --------------------- Less strict clump, without MHC  (output/noMHC_) --------
################################################################################
# ANM --> asthma
Rscript "${SCRIPT}" --exp_gwas "${ANM}" --exp_name ANM_ofh_meta "${ANM_EXP[@]}" \
  --exp_ivs "${ANM_IVS}" \
  --out_gwas "${ASTHMA}" --out_name asthma "${ASTHMA_OUT[@]}" \
  --ld_ref "${LDREF}" --out_prefix "${OUT}noMHC_" "${PROXY[@]}" --presso_nbdist 5000 --analysis_exclude_mhc \
  || echo "WARNING: [lenient-noMHC] ANM->asthma failed; continuing."

# Asthma --> ANM
Rscript "${SCRIPT}" --exp_gwas "${ASTHMA}" --exp_name asthma "${ASTHMA_EXP[@]}" \
  --out_gwas "${ANM}" --out_name ANM_ofh_meta "${ANM_OUT[@]}" \
  --ld_ref "${LDREF}" --out_prefix "${OUT}noMHC_" --clump_r2 0.2 --clump_kb 10000 "${PROXY[@]}" --presso_nbdist 10000 --analysis_exclude_mhc \
  || echo "WARNING: [lenient-noMHC] asthma->ANM failed; continuing."

################################################################################
# --------------------- More strict clump, without MHC  (output/noMHC_rigid_) --
################################################################################
# ANM --> asthma
Rscript "${SCRIPT}" --exp_gwas "${ANM}" --exp_name ANM_ofh_meta "${ANM_EXP[@]}" \
  --out_gwas "${ASTHMA}" --out_name asthma "${ASTHMA_OUT[@]}" \
  --ld_ref "${LDREF}" --out_prefix "${OUT}noMHC_rigid_" --clump_r2 0.001 --clump_kb 10000 "${PROXY[@]}" --presso_nbdist 5000 --analysis_exclude_mhc \
  || echo "WARNING: [strict-noMHC] ANM->asthma failed; continuing."

# Asthma --> ANM
Rscript "${SCRIPT}" --exp_gwas "${ASTHMA}" --exp_name asthma "${ASTHMA_EXP[@]}" \
  --out_gwas "${ANM}" --out_name ANM_ofh_meta "${ANM_OUT[@]}" \
  --ld_ref "${LDREF}" --out_prefix "${OUT}noMHC_rigid_" --clump_r2 0.001 --clump_kb 10000 "${PROXY[@]}" --presso_nbdist 10000 --analysis_exclude_mhc \
  || echo "WARNING: [strict-noMHC] asthma->ANM failed; continuing."

echo "All MR runs complete."
