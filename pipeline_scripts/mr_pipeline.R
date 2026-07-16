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
  # Read exposure and select instruments once per exposure.
  # exposure_label = unique tag for output filenames; exposure_trait = the
  # label shown in the results (uses --exp_name when a single exposure is given).
  exposure_label <- tools::file_path_sans_ext(basename(exposure_file))
  exposure_trait <- if (length(exposure_files) == 1 && !is.null(opt$exp_name) && opt$exp_name != "Exposure") opt$exp_name else exposure_label
  exp_col_args <- list(
    snp = opt$exp_snp, beta = opt$exp_beta, se = opt$exp_se,
    ea = opt$exp_ea, nea = opt$exp_nea, p = opt$exp_p,
    eaf = opt$exp_eaf, n = opt$exp_n, chr = opt$exp_chr, pos = opt$exp_pos
  )
  # Read + clean the exposure (cheap; pre-filtered by p-value so format_data
  # later runs only on the selected instruments). Full summary stats are the
  # assumed input; --skip_clump handles pre-independent signal lists.
  exposure_raw <- read_gwas(
    gwas_file = exposure_file,
    type = "exposure",
    col_args = exp_col_args,
    trait_name = exposure_trait,
    tmp_dir = opt$tmp_dir,
    pval_thresh = opt$clump_p
  )
  message(sprintf("Successfully read exposure data for trait '%s'.", exposure_trait))

  # All instrument selection (p-filter -> clump/skip -> format -> F-stat -> MHC)
  # lives in select_instruments().
  exposure_ivs_dat <- select_instruments(
    exposure_raw = exposure_raw,
    trait_name = exposure_trait,
    clump_p = opt$clump_p,
    clump_kb = opt$clump_kb,
    clump_r2 = opt$clump_r2,
    ld_ref = opt$ld_ref,
    plink_bin = opt$plink_bin,
    min_f_stat = opt$f_stat,
    skip_clump = opt$skip_clump,
    mhc_region = opt$mhc_region,
    exclude_mhc = opt$exclude_mhc
  )
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
    paste0(opt$out_prefix, exposure_label, "_exposure_ivs.tsv"),
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
    outcome_label <- tools::file_path_sans_ext(basename(outcome_file))
    outcome_trait <- if (length(outcome_files) == 1 && !is.null(opt$out_name) && opt$out_name != "Outcome") opt$out_name else outcome_label
    this_out_prefix <- paste0(opt$out_prefix, exposure_label, "_vs_", outcome_label)

    out_col_args <- list(
      snp = opt$out_snp, beta = opt$out_beta, se = opt$out_se,
      ea = opt$out_ea, nea = opt$out_nea, p = opt$out_p,
      eaf = opt$out_eaf, n = opt$out_n, ncase = opt$out_ncase,
      chr = opt$out_chr, pos = opt$out_pos
    )
    cur_mr_res <- run_mr_analysis(exposure_ivs_dat, outcome_file, outcome_trait, out_col_args, this_out_prefix, opt)
    
    all_mr_results[[length(all_mr_results) + 1]] <- cur_mr_res
  }
}

all_mr_dt <- data.table::rbindlist(all_mr_results, fill = TRUE)
process_mr_results(all_mr_dt, opt)

message("All analyses complete.")

