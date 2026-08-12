#!/bin/bash
#SBATCH --job-name=mr_eos
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G
#SBATCH --partition=icelake        # <-- change to your cluster's partition
#SBATCH --account=MRC-EPID-SL0-CPU # <-- change to your project/account

# Eosinophil count (EO) <-> asthma / ANM / POI, bidirectional, full 4-config grid
# (lenient/strict clumping x with/without MHC). Results go to output/eos/.
set -uo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mr_env
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# ---- Paths ------------------------------------------------------------------
PIPE=/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline
SCRIPT="${PIPE}/pipeline_scripts/mr_pipeline.R"
LDREF=/rfs/project/rfs-mpB3sSsgAn4/Studies/1000GP_Phase3_genotyped/superpopulation/Phase3_eur
OUT_EOS="${PIPE}/output/eos/"; mkdir -p "${OUT_EOS}"
NBDIST=5000                        # MR-PRESSO simulations (see runtime note at end)

# >>> SET THESE: full path to the eosinophil GWAS + its signal list <<<
# EOS = FULL eosinophil sumstats (needed so proxies can be found).
# EOS_IVS = the signal-ID list (eosinophils_count_signals.txt). The IDs in it must
#           match the --exp_snp column below AND be rsIDs (for LD clumping vs the
#           1000G reference and harmonisation with asthma/ANM/POI, which use rsIDs).
EOS=/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/SET_ME/eosinophils_count.tsv
EOS_IVS=/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/SET_ME/eosinophils_count_signals.txt

ASTHMA=/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/autoimmune_traits/data/asthma_hg19_lifted_with_n.tsv.gz
ANM=/rfs/project/rfs-ic9ZbslBEP4/Datasets/Projects/ANM_2025/g2g/anm_ofh_ukb_23andMe_metal/anm_ofh_ukb_23andMe_metal_gwas_funcannoQC.txt
ANM_IVS=/rfs/project/rfs-ic9ZbslBEP4/Datasets/Projects/ANM_2025/g2g/anm_ofh_ukb_23andMe_metal/p01_signal_selection/anm_ofh_ukb_23andMe_metal_all_independent_signals_GWASinfo.txt
POI=/rfs/project/rfs-ic9ZbslBEP4/Datasets/Projects/POI/g2g/poi_ofh_ukb_23andMe_metal/poi_ofh_ukb_23andMe_metal_gwas_funcannoQC.txt
POI_IVS=/rfs/project/rfs-ic9ZbslBEP4/Datasets/Projects/POI/g2g/poi_ofh_ukb_23andMe_metal/p01_signal_selection/poi_ofh_ukb_23andMe_metal_all_independent_signals_GWASinfo.txt

# ---- Column mappings (arrays keep '#CHR' safe from being read as a comment) --
# Eosinophils: continuous, b37, N=408112
EO_EXP=(--exp_snp variant_id --exp_chr chromosome --exp_pos base_pair_location
        --exp_ea effect_allele --exp_nea other_allele --exp_beta beta --exp_se standard_error
        --exp_p p_value --exp_eaf effect_allele_frequency --exp_n_total 408112)
EO_OUT=(--out_snp variant_id --out_chr chromosome --out_pos base_pair_location
        --out_ea effect_allele --out_nea other_allele --out_beta beta --out_se standard_error
        --out_p p_value --out_eaf effect_allele_frequency --out_n_total 408112 --out_type continuous)
# Asthma: binary
ASTHMA_EXP=(--exp_snp rsid --exp_chr '#CHR' --exp_pos POS --exp_ea ALT --exp_nea REF
            --exp_beta inv_var_meta_beta --exp_se inv_var_meta_sebeta --exp_p inv_var_meta_p
            --exp_eaf all_meta_AF --exp_n n --exp_ncase N_case)
ASTHMA_OUT=(--out_snp rsid --out_chr '#CHR' --out_pos POS --out_ea ALT --out_nea REF
            --out_beta inv_var_meta_beta --out_se inv_var_meta_sebeta --out_p inv_var_meta_p
            --out_eaf all_meta_AF --out_n n --out_ncase N_case)
# ANM: continuous. Exposure p = P_BOLT_LMM; outcome p = P_BOLT_LMM_INF (populated for ANM)
ANM_EXP=(--exp_snp SNP --exp_chr CHR --exp_pos BP --exp_ea ALLELE1 --exp_nea ALLELE0
         --exp_beta BETA --exp_se SE --exp_p P_BOLT_LMM --exp_eaf A1FREQ --exp_n_total 679746)
ANM_OUT=(--out_snp SNP --out_chr CHR --out_pos BP --out_ea ALLELE1 --out_nea ALLELE0
         --out_beta BETA --out_se SE --out_p P_BOLT_LMM_INF --out_eaf A1FREQ
         --out_n_total 679746 --out_type continuous)
# POI: binary. Both exposure and outcome p = P_BOLT_LMM (INF column is empty for POI)
POI_EXP=(--exp_snp SNP --exp_chr CHR --exp_pos BP --exp_ea ALLELE1 --exp_nea ALLELE0
         --exp_beta BETA --exp_se SE --exp_p P_BOLT_LMM --exp_eaf A1FREQ
         --exp_n_total 619464 --exp_ncase_total 22254)
POI_OUT=(--out_snp SNP --out_chr CHR --out_pos BP --out_ea ALLELE1 --out_nea ALLELE0
         --out_beta BETA --out_se SE --out_p P_BOLT_LMM --out_eaf A1FREQ
         --out_n_total 619464 --out_ncase_total 22254)
PROXY=(--proxies --proxy_r2 0.8 --proxy_kb 1000)

cd "${PIPE}"

run() {  # run <label> <Rscript args...>   (non-fatal: one failure won't stop the grid)
  local label="$1"; shift
  echo "=== [$(date +%T)] ${label} ==="
  Rscript "${SCRIPT}" "$@" || echo "WARNING: ${label} failed; continuing."
}

# grid_clump : exposure is clumped (lenient r2=0.2, strict r2=0.001). Used when the
#              exposure has no signal list (eosinophils, asthma).
# The 4 configs: prefix | selection args | mhc flag
grid_clump() {   # $1 exp_gwas $2 exp_name $3 EXP_ARR(name) $4 out_gwas $5 out_name $6 OUT_ARR(name) $7 tag
  local eg="$1" en="$2" ea="$3" og="$4" on="$5" oa="$6" tag="$7"
  eval "local EXPA=(\"\${$ea[@]}\")"; eval "local OUTA=(\"\${$oa[@]}\")"
  local specs=( "|--clump_r2 0.2 --clump_kb 10000|"
                "rigid_|--clump_r2 0.001 --clump_kb 10000|"
                "noMHC_|--clump_r2 0.2 --clump_kb 10000|--analysis_exclude_mhc"
                "noMHC_rigid_|--clump_r2 0.001 --clump_kb 10000|--analysis_exclude_mhc" )
  for s in "${specs[@]}"; do IFS='|' read -r pref sel mhc <<< "$s"
    run "${tag} [${pref:-lenient+MHC}]" \
      --exp_gwas "$eg" --exp_name "$en" "${EXPA[@]}" \
      --out_gwas "$og" --out_name "$on" "${OUTA[@]}" \
      --ld_ref "$LDREF" --out_prefix "${OUT_EOS}${pref}" $sel "${PROXY[@]}" --presso_nbdist "$NBDIST" $mhc
  done
}

# grid_ivs : exposure uses a pre-defined signal list for the lenient config, and
#            r2<0.001 re-clumping for the strict config (ANM, POI).
grid_ivs() {     # $1 exp_gwas $2 exp_name $3 EXP_ARR(name) $4 ivs_file $5 out_gwas $6 out_name $7 OUT_ARR(name) $8 tag
  local eg="$1" en="$2" ea="$3" ivs="$4" og="$5" on="$6" oa="$7" tag="$8"
  eval "local EXPA=(\"\${$ea[@]}\")"; eval "local OUTA=(\"\${$oa[@]}\")"
  local specs=( "|--exp_ivs ${ivs}|"
                "rigid_|--clump_r2 0.001 --clump_kb 10000|"
                "noMHC_|--exp_ivs ${ivs}|--analysis_exclude_mhc"
                "noMHC_rigid_|--clump_r2 0.001 --clump_kb 10000|--analysis_exclude_mhc" )
  for s in "${specs[@]}"; do IFS='|' read -r pref sel mhc <<< "$s"
    run "${tag} [${pref:-lenient+MHC}]" \
      --exp_gwas "$eg" --exp_name "$en" "${EXPA[@]}" \
      --out_gwas "$og" --out_name "$on" "${OUTA[@]}" \
      --ld_ref "$LDREF" --out_prefix "${OUT_EOS}${pref}" $sel "${PROXY[@]}" --presso_nbdist "$NBDIST" $mhc
  done
}

# ---- EO <-> asthma  (eos exposure uses its signal list for lenient) ---------
grid_ivs   "$EOS"    EO     EO_EXP     "$EOS_IVS" "$ASTHMA" asthma ASTHMA_OUT "eos->asthma"
grid_clump "$ASTHMA" asthma ASTHMA_EXP            "$EOS"    EO     EO_OUT     "asthma->eos"

# ---- EO <-> ANM -------------------------------------------------------------
grid_ivs   "$EOS" EO           EO_EXP  "$EOS_IVS" "$ANM" ANM_ofh_meta ANM_OUT "eos->ANM"
grid_ivs   "$ANM" ANM_ofh_meta ANM_EXP "$ANM_IVS" "$EOS" EO           EO_OUT  "ANM->eos"

# ---- EO <-> POI -------------------------------------------------------------
grid_ivs   "$EOS" EO           EO_EXP  "$EOS_IVS" "$POI" POI_ofh_meta POI_OUT "eos->POI"
grid_ivs   "$POI" POI_ofh_meta POI_EXP "$POI_IVS" "$EOS" EO           EO_OUT  "POI->eos"

echo "All EO MR runs complete. Results in ${OUT_EOS}"
