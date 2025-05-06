Rscript pipeline_scripts/mr_pipeline.R \
  --exposure_dir /lustre/groups/itg/teams/zeggini/projects/MR_pipeline/test_exposure \
  --outcome_dir /lustre/groups/itg/teams/zeggini/projects/MR_pipeline/test_outcome \
  --ld_ref /lustre/groups/itg/shared/referenceData/1kG/EUR_1kg_v3/EUR \
  --out_prefix test_output/ \
  --plink_bin /lustre/groups/itg/teams/zeggini/projects/MR_pipeline/plink \
  \
  --exp_snp "rsID" \
  --exp_beta "beta" \
  --exp_se "se" \
  --exp_ea "EA" \
  --exp_nea "NEA" \
  --exp_p "pval" \
  --exp_eaf "EAF" \
  --exp_n "N" \
  --exp_ncase "Ncases" \
  --exp_type "binary" \
  --exp_chr "CHR.b37" \
  --exp_pos "POS.b37" \
  \
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
  --presso FALSE