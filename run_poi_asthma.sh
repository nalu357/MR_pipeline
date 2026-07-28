cd /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline
#Rscript install_dependencies.R
  
# POI --> asthma  
Rscript /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline/pipeline_scripts/mr_pipeline.R \
  --exp_gwas /rfs/project/rfs-ic9ZbslBEP4/Datasets/Projects/POI/g2g/poi_finngen_metal_v2/p01_signal_selection/poi_finngen_metal_v2_all_independent_signals_GWASinfo.txt \
  --exp_name POI \
  --exp_snp SNP --exp_chr CHR --exp_pos BP \
  --exp_ea ALLELE1 --exp_nea ALLELE0 --exp_beta BETA --exp_se SE \
  --exp_p P_BOLT_LMM --exp_eaf A1FREQ --exp_n_total 739450 \
  --out_gwas /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/autoimmune_traits/data/asthma_hg19_lifted_with_n.tsv.gz \
  --out_name asthma \
  --out_snp rsid --out_chr '#CHR' --out_pos POS --out_n n \
  --out_ea ALT --out_nea REF --out_beta inv_var_meta_beta --out_se inv_var_meta_sebeta \
  --out_p inv_var_meta_p --out_eaf all_meta_AF \
  --ld_ref /rfs/project/rfs-mpB3sSsgAn4/Studies/1000GP_Phase3_genotyped/superpopulation/Phase3_eur \
  --out_prefix /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline/output/ \
  --clump_r2 0.001 --clump_kb 10000
  
# Asthma --> POI
Rscript /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline/pipeline_scripts/mr_pipeline.R \
  --exp_gwas /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/autoimmune_traits/data/asthma_hg19_lifted_with_n.tsv.gz \
  --exp_name asthma \
  --exp_snp rsid --exp_chr '#CHR' --exp_pos POS \
  --exp_ea ALT --exp_nea REF --exp_beta inv_var_meta_beta --exp_se inv_var_meta_sebeta \
  --exp_p inv_var_meta_p --exp_eaf all_meta_AF --exp_n n \
  --out_gwas /rfs/project/rfs-ic9ZbslBEP4/Datasets/Projects/POI/g2g/poi_finngen_metal_v2/poi_finngen_metal_v2_gwas_funcannoQC.txt \
  --out_name POI \
  --out_snp SNP --out_chr CHR --out_pos BP \
  --out_ea ALLELE1 --out_nea ALLELE0 --out_beta BETA --out_se SE \
  --out_p P_BOLT_LMM --out_eaf A1FREQ --out_n_total 739450 \
  --ld_ref /rfs/project/rfs-mpB3sSsgAn4/Studies/1000GP_Phase3_genotyped/superpopulation/Phase3_eur \
  --out_prefix /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline/output/ \
  --clump_r2 0.001 --clump_kb 10000
  

# POI --> asthma  
Rscript /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline/pipeline_scripts/mr_pipeline.R \
  --exp_gwas /rfs/project/rfs-ic9ZbslBEP4/Datasets/Projects/POI/g2g/poi_finngen_metal_v2/p01_signal_selection/poi_finngen_metal_v2_all_independent_signals_GWASinfo.txt \
  --exp_name POI \
  --exp_snp SNP --exp_chr CHR --exp_pos BP \
  --exp_ea ALLELE1 --exp_nea ALLELE0 --exp_beta BETA --exp_se SE \
  --exp_p P_BOLT_LMM --exp_eaf A1FREQ --exp_n_total 739450 \
  --out_gwas /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/autoimmune_traits/data/asthma_hg19_lifted_with_n.tsv.gz \
  --out_name asthma \
  --out_snp rsid --out_chr '#CHR' --out_pos POS --out_n n \
  --out_ea ALT --out_nea REF --out_beta inv_var_meta_beta --out_se inv_var_meta_sebeta \
  --out_p inv_var_meta_p --out_eaf all_meta_AF \
  --ld_ref /rfs/project/rfs-mpB3sSsgAn4/Studies/1000GP_Phase3_genotyped/superpopulation/Phase3_eur \
  --out_prefix /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline/output/noMHC_  \
  --clump_r2 0.001 --clump_kb 10000 --exclude_mhc
  
# Asthma --> POI
Rscript /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline/pipeline_scripts/mr_pipeline.R \
  --exp_gwas /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/autoimmune_traits/data/asthma_hg19_lifted_with_n.tsv.gz \
  --exp_name asthma \
  --exp_snp rsid --exp_chr '#CHR' --exp_pos POS \
  --exp_ea ALT --exp_nea REF --exp_beta inv_var_meta_beta --exp_se inv_var_meta_sebeta \
  --exp_p inv_var_meta_p --exp_eaf all_meta_AF --exp_n n \
  --out_gwas /rfs/project/rfs-ic9ZbslBEP4/Datasets/Projects/POI/g2g/poi_finngen_metal_v2/poi_finngen_metal_v2_gwas_funcannoQC.txt \
  --out_name POI \
  --out_snp SNP --out_chr CHR --out_pos BP \
  --out_ea ALLELE1 --out_nea ALLELE0 --out_beta BETA --out_se SE \
  --out_p P_BOLT_LMM --out_eaf A1FREQ --out_n_total 739450 \
  --ld_ref /rfs/project/rfs-mpB3sSsgAn4/Studies/1000GP_Phase3_genotyped/superpopulation/Phase3_eur \
  --out_prefix /rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline/output/noMHC_ \
  --clump_r2 0.001 --clump_kb 10000 --exclude_mhc