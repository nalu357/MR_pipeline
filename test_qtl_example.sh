Rscript pipeline_scripts/mr_pipeline.R \
 --qtl \
 --exp_gwas /lustre/groups/itg/teams/zeggini/projects/AD_T2D/data/AD_merged.txt \
 --out_gwas /lustre/groups/itg/teams/zeggini/projects/AD_T2D/data/T2D_european.txt\
 --ld_ref /lustre/groups/itg/shared/referenceData/1kG/EUR_1kg_v3/EUR \
 --out_prefix output/ \
 --plink_bin /lustre/groups/itg/teams/zeggini/projects/MR_pipeline/plink \
 \
 --exp_name "AD" \
 --exp_snp "rsID" \
 --exp_beta "beta" \
 --exp_se "se" \
 --exp_ea "EA" \
 --exp_nea "NEA" \
 --exp_p "pval" \
 --exp_n "N" \
 --exp_ncase "Ncases" \
 --exp_type "binary" \
 --exp_chr "CHR.b37" \
 --exp_pos "POS.b37" \
 \
 --out_name "T2D" \
 --out_snp "rsID" \
 --out_beta "beta" \
 --out_se "se" \
 --out_ea "EA" \
 --out_nea "NEA" \
 --out_p "pval" \
 --out_eaf "EAF" \
 --out_n "N" \
 --out_ncase "Ncases" \
 --out_type "binary" \
 --out_chr "CHR.b37" \
 --out_pos "POS.b37" \
 \
 --clump_p 5e-8 \
 --clump_r2 0.001 \
 --clump_kb 10000 \
 --f_stat 10 \
 --steiger \
 --presso
 
 
  # --exp_eaf "freq_effect" \