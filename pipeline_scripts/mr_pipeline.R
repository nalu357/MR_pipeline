#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# Portable start-up (no hardcoded machine paths).
# The script resolves its own location so it can be launched from anywhere
# (e.g. `Rscript /path/to/mr_pipeline.R ...` on any cluster) and still find
# its helper functions. An optional --lib_path lets you point at a custom R
# library without editing the script.
# ---------------------------------------------------------------------------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  getwd()
}
SCRIPT_DIR <- get_script_dir()

# Resolve the R library BEFORE loading packages (optparse runs later).
# Precedence: --lib_path (explicit) > repo-local Rpackages/ > R default.
# The repo-local library lets you `git clone`, install deps into Rpackages/,
# and run with zero configuration (see install_dependencies.R).
.raw_args <- commandArgs(trailingOnly = TRUE)
.libpath_idx <- which(.raw_args == "--lib_path")
.repo_lib_candidates <- c(file.path(dirname(SCRIPT_DIR), "Rpackages"),
                          file.path(SCRIPT_DIR, "Rpackages"))
.repo_lib <- .repo_lib_candidates[dir.exists(.repo_lib_candidates)][1]
if (length(.libpath_idx) == 1 && length(.raw_args) >= .libpath_idx + 1) {
  .user_lib <- .raw_args[.libpath_idx + 1]
  if (dir.exists(.user_lib)) {
    .libPaths(.user_lib)
    message(sprintf("Using R library path (--lib_path): %s", .user_lib))
  } else {
    warning(sprintf("--lib_path '%s' does not exist; ignoring.", .user_lib), call. = FALSE)
  }
} else if (!is.na(.repo_lib)) {
  .libPaths(.repo_lib)
  message(sprintf("Using repo-local R library: %s", .repo_lib))
}

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
})

option_list <- list(
  make_option("--exp_gwas", type="character", help="Path to Exposure GWAS summary statistics file (optional if --exposure_dir is set)."),
  make_option("--exposure_dir", type="character", default=NULL, help="Directory containing multiple exposure GWAS files (optional)."),
  make_option("--out_gwas", type="character", help="Path to Outcome GWAS summary statistics file (optional if --outcome_dir is set)."),
  make_option("--outcome_dir", type="character", default=NULL, help="Directory containing multiple outcome GWAS files (optional)."),
  make_option("--out_prefix", type="character", help="Prefix for output files (required)."),
  # Exposure column names
  make_option("--exp_name", type="character", default="Exposure", help="Name for the exposure trait"),
  make_option("--exp_snp", type="character", default="rs_id", help="Exposure: SNP ID column name"),
  make_option("--exp_beta", type="character", default="beta", help="Exposure: Beta column name"),
  make_option("--exp_se", type="character", default="standard_error", help="Exposure: Standard Error column name"),
  make_option("--exp_ea", type="character", default="effect_allele", help="Exposure: Effect Allele column name"),
  make_option("--exp_nea", type="character", default="other_allele", help="Exposure: Non-Effect Allele column name"),
  make_option("--exp_p", type="character", default="p_value", help="Exposure: P-value column name"),
  make_option("--exp_eaf", type="character", default="effect_allele_frequency", help="Exposure: Effect Allele Frequency column name"),
  make_option("--exp_n", type="character", default="n", help="Exposure: Sample Size column name"),
  make_option(c("--exp_ncase"), type="character", default="ncases", help="Exposure: Number of Cases column name (Required for binary outcome + Steiger) [default: %default]"),
  make_option(c("--exp_type"), type="character", default="binary", help="Exposure type ('binary' or 'continuous') [default: %default]"),
  make_option("--exp_chr", type="character", default="chr", help="Exposure: Chromosome column name"),
  make_option("--exp_pos", type="character", default="pos", help="Exposure: Position column name"),
  # Outcome column names
  make_option("--out_name", type="character", default="Outcome", help="Name for the outcome trait"),
  make_option("--out_snp", type="character", default="rs_id", help="Outcome: SNP ID column name"),
  make_option("--out_beta", type="character", default="beta", help="Outcome: Beta column name"),
  make_option("--out_se", type="character", default="standard_error", help="Outcome: Standard Error column name"),
  make_option("--out_ea", type="character", default="effect_allele", help="Outcome: Effect Allele column name"),
  make_option("--out_nea", type="character", default="other_allele", help="Outcome: Non-Effect Allele column name"),
  make_option("--out_p", type="character", default="p_value", help="Outcome: P-value column name"),
  make_option("--out_eaf", type="character", default="effect_allele_frequency", help="Outcome: Effect Allele Frequency column name"),
  make_option("--out_n", type="character", default="n", help="Outcome: Sample Size column name"),
  make_option("--out_ncase", type="character", default="ncases", help="Outcome: Number of Cases column name"),
  make_option("--out_type", type="character", default="binary", help="Outcome type ('binary' or 'continuous')"),
  make_option("--out_chr", type="character", default="chr", help="Outcome: Chromosome column name"),
  make_option("--out_pos", type="character", default="pos", help="Outcome: Position column name"),
  # Clumping parameters
  make_option("--clump_p", type="numeric", default=5e-8, help="P-value threshold for IV selection"),
  make_option("--clump_kb", type="numeric", default=10000, help="Kilobase radius for LD clumping"),
  make_option("--clump_r2", type="numeric", default=0.001, help="R-squared threshold for LD clumping"),
  make_option("--ld_ref", type="character", help="Path prefix to LD reference panel (PLINK bfile format)"),
  make_option("--plink_bin", type="character", default=NULL, help="Path to PLINK binary"),
  make_option("--skip_clump", action="store_true", default=FALSE, help="Skip LD clumping (ONLY for pre-clumped/QTL instruments already independent at r2<0.001). NOT for gene-mapping/COJO signal lists, which are only independent at looser r2 - see README."),
  make_option("--lib_path", type="character", default=NULL, help="Optional path to a custom R library (applied before packages load)."),
  # MHC handling (long-range LD + heavy pleiotropy). Default: keep but flag.
  make_option("--mhc_region", type="character", default="6:25000000-34000000",
              help="MHC region as CHR:START-END used to flag instruments [default: %default, GRCh37 extended MHC]. Set to match your GWAS build."),
  make_option("--exclude_mhc", action="store_true", default=FALSE,
              help="Drop instruments in --mhc_region instead of just flagging them [default: keep and flag]."),
  # Analysis parameters
  make_option("--f_stat", type="numeric", default=10, help="Minimum F-statistic for IVs"),
  # Output options
  make_option("--tmp_dir", type="character", default="./tmp_pipeline", help="Temporary directory for intermediate files"),
  # Sensitivity analysis options
  make_option("--steiger", type="logical", action="store_true", default=TRUE, help="Run Steiger filtering"),
  make_option("--no_steiger", action="store_false", dest="steiger"),
  make_option("--presso", type="logical", action="store_true", default=TRUE, help="Run MR-PRESSO outlier test"),
  make_option("--no_presso", action="store_false", dest="presso")
)
parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

# Source helper functions relative to this script's own location first, so the
# pipeline works regardless of the current working directory.
functions_candidates <- c(
  file.path(SCRIPT_DIR, "pipeline_functions.R"),
  "pipeline_functions.R",
  file.path("pipeline_scripts", "pipeline_functions.R")
)
functions_file <- functions_candidates[file.exists(functions_candidates)][1]
if (is.na(functions_file)) {
  stop(sprintf("Helper functions file 'pipeline_functions.R' not found in any of: %s",
               paste(functions_candidates, collapse=", ")), call.=FALSE)
}
source(functions_file)

if (!is.null(opt$presso) && opt$presso) {
  if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
    warning("MR-PRESSO analysis requested (--presso) but package 'MRPRESSO' is not installed. Skipping MR-PRESSO.", call. = FALSE)
    opt$presso <- FALSE
  } else {
    suppressPackageStartupMessages(library(MRPRESSO))
  }
}

if (!is.null(opt$exposure_dir)) {
  exposure_files <- list.files(
    opt$exposure_dir,
    pattern = "\\.(tsv|csv|txt|tsv\\.gz|csv\\.gz|txt\\.gz|zip|gz)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
} else if (!is.null(opt$exp_gwas)) {
  exposure_files <- opt$exp_gwas
} else {
  stop("You must provide either --exp_gwas or --exposure_dir.", call. = FALSE)
}

all_mr_results <- list()
for (exposure_file in exposure_files) {
  # Load and clump exposure ONCE per exposure
  exposure_name <- tools::file_path_sans_ext(basename(exposure_file))
  exp_col_args <- list(
    snp = opt$exp_snp, beta = opt$exp_beta, se = opt$exp_se,
    ea = opt$exp_ea, nea = opt$exp_nea, p = opt$exp_p,
    eaf = opt$exp_eaf, n = opt$exp_n, chr = opt$exp_chr, pos = opt$exp_pos
  )
  exposure_dat <- load_and_format_gwas(
    gwas_file = exposure_file,
    type = "exposure",
    col_args = exp_col_args,
    trait_name = exposure_name,
    tmp_dir = opt$tmp_dir
  )
  if (!inherits(exposure_dat, "data.frame") || nrow(exposure_dat) == 0) {
    stop("Exposure data loading and formatting failed or resulted in empty data frame. Exiting.", call. = FALSE)
  }
  message(sprintf("Successfully loaded and formatted exposure data for trait '%s'.", opt$exp_name))
  
  if (opt$skip_clump) {
    message("----- Skipping LD Clumping (--skip_clump specified) -----")
    warning("--skip_clump assumes the input instruments are ALREADY independent at MR standard (r2<0.001). ",
            "Gene-mapping / GCTA-COJO signal lists are typically only independent at r2<0.05 and are NOT ",
            "safe to use this way - correlated instruments understate IVW standard errors. If in doubt, ",
            "re-run WITHOUT --skip_clump to LD-clump at r2<0.001.", call. = FALSE)
    message(sprintf("Selecting IVs based on p < %g and F-statistic >= %f", opt$clump_p, opt$f_stat))
    exposure_ivs_dat <- exposure_dat %>%
      filter(pval.exposure < opt$clump_p)
    if (nrow(exposure_ivs_dat) == 0) {
      stop(sprintf("No SNPs found below the significance threshold p < %g.", opt$clump_p), call. = FALSE)
    }
    message(sprintf("Found %d SNPs below p-value threshold.", nrow(exposure_ivs_dat)))
    warn_non_rsid_instruments(exposure_ivs_dat)
    if (!"samplesize.exposure" %in% names(exposure_ivs_dat)) {
      message("Warning: Sample size column ('samplesize.exposure') not found. Cannot calculate F-statistic.")
      message("Warning: Proceeding without F-statistic filtering as sample size is missing.")
      exposure_ivs_dat$F_statistic <- NA
    } else {
      exposure_ivs_dat <- exposure_ivs_dat %>%
        mutate(F_statistic = (beta.exposure^2) / (se.exposure^2))
      rows_before_f <- nrow(exposure_ivs_dat)
      exposure_ivs_dat <- exposure_ivs_dat %>%
        filter(F_statistic >= opt$f_stat)
      rows_removed_f <- rows_before_f - nrow(exposure_ivs_dat)
      if (rows_removed_f > 0) {
        message(sprintf("Removed %d IVs with F-statistic < %f.", rows_removed_f, opt$f_stat))
      }
      if (nrow(exposure_ivs_dat) == 0) {
        stop(sprintf("No IVs remained after F-statistic filtering (F >= %f).", opt$f_stat), call.=FALSE)
      }
      message(sprintf("%d IVs remain after F-statistic filtering.", nrow(exposure_ivs_dat)))
    }
    exposure_ivs_dat <- flag_mhc_instruments(exposure_ivs_dat, opt$mhc_region, opt$exclude_mhc)
    message("----- Finished Selecting Instruments (Clumping Skipped) -----")
  } else {
    message("----- Selecting and Clumping Instruments (LD Clumping Enabled) -----")
    exposure_ivs_dat <- clump_and_filter_ivs(
      exposure_dat = exposure_dat,
      clump_p = opt$clump_p,
      clump_kb = opt$clump_kb,
      clump_r2 = opt$clump_r2,
      ld_ref = opt$ld_ref,
      plink_bin = opt$plink_bin,
      min_f_stat = opt$f_stat,
      mhc_region = opt$mhc_region,
      exclude_mhc = opt$exclude_mhc
    )
  }
  if (!inherits(exposure_ivs_dat, "data.frame") || nrow(exposure_ivs_dat) == 0) {
    stop("No exposure instruments (IVs) remained after clumping and filtering. Exiting.", call. = FALSE)
  }
  exposure_ivs <- exposure_ivs_dat$SNP
  message(sprintf("Final selected IVs count: %d", length(exposure_ivs)))
  # Write a fuller instrument table (not just SNP IDs) for transparency/reporting.
  iv_cols <- intersect(
    c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure",
      "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
      "F_statistic", "mhc"),
    names(exposure_ivs_dat))
  data.table::fwrite(
    data.table::as.data.table(exposure_ivs_dat)[, ..iv_cols],
    paste0(opt$out_prefix, exposure_name, "_exposure_ivs.tsv"),
    sep = "\t", na = "NA")
  
  if (!is.null(opt$outcome_dir)) {
    outcome_files <- list.files(
      opt$outcome_dir,
      pattern = "\\.(tsv|csv|txt|tsv\\.gz|csv\\.gz|txt\\.gz|zip|gz)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    if(length(outcome_files)==0) stop("The outcome directory is empty or contains an unknown data type.")
  } else if (!is.null(opt$out_gwas)) {
    outcome_files <- opt$out_gwas
  } else {
    stop("You must provide either --out_gwas or --outcome_dir.", call. = FALSE)
  }
  
  for (outcome_file in outcome_files) {
    outcome_name <- tools::file_path_sans_ext(basename(outcome_file))  # gsub(".*", "", outcome_file,)
    this_out_prefix <- paste0(opt$out_prefix, exposure_name, "_vs_", outcome_name)
    
    out_col_args <- list(
      snp = opt$out_snp, beta = opt$out_beta, se = opt$out_se,
      ea = opt$out_ea, nea = opt$out_nea, p = opt$out_p,
      eaf = opt$out_eaf, n = opt$out_n, ncase = opt$out_ncase,
      chr = opt$out_chr, pos = opt$out_pos
    )
    cur_mr_res <- run_mr_analysis(exposure_ivs_dat, outcome_file, outcome_name, out_col_args, this_out_prefix, opt)
    
    all_mr_results[[length(all_mr_results) + 1]] <- cur_mr_res
  }
}

all_mr_dt <- data.table::rbindlist(all_mr_results, fill = TRUE)
process_mr_results(all_mr_dt, opt)

message("All analyses complete.")

# all_mr_results <- list()
# for(i in c("T2D_merged_vs_PD_merged", "T2D_merged_vs_AD_merged", "LBD_merged_vs_PD_merged", "LBD_merged_vs_AD_merged")){
#   all_mr_results[[length(all_mr_results) + 1]] <- data.table::fread(paste0("test_output/", i, "_full_mr_results.csv"))
# }


