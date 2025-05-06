#!/usr/bin/env Rscript

# Mendelian Randomization Pipeline Script
setwd("/lustre/groups/itg/teams/zeggini/projects/MR_pipeline")
.libPaths("Rpackages/")
#------------------- Libraries -----------------#
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr) # Needed for clumping function call logic below (plink binary finding)
})
# Conditionally load MRPRESSO if requested and available
# Note: opt object does not exist yet, this needs to be moved after parse_args
# We will handle library loading within the function call or after args are parsed.

#------------------- Command Line Arguments -----------------#

option_list <- list(
  # Input files
  make_option(c("--exp_gwas"), type="character", help="Path to Exposure GWAS summary statistics file (required). See README for format details."),
  make_option(c("--out_gwas"), type="character", help="Path to Outcome GWAS summary statistics file (required). See README for format details."),
  
  # Exposure column names (Map these to the actual headers in your exposure file)
  make_option(c("--exp_name"), type="character", default="Exposure", help="Name for the exposure trait"),
  make_option(c("--exp_snp"), type="character", default="rs_id", help="Exposure: SNP ID column name [default: %default]"),
  make_option(c("--exp_beta"), type="character", default="beta", help="Exposure: Beta column name [default: %default]"),
  make_option(c("--exp_se"), type="character", default="standard_error", help="Exposure: Standard Error column name [default: %default]"),
  make_option(c("--exp_ea"), type="character", default="effect_allele", help="Exposure: Effect Allele column name [default: %default]"),
  make_option(c("--exp_nea"), type="character", default="other_allele", help="Exposure: Non-Effect Allele column name [default: %default]"),
  make_option(c("--exp_p"), type="character", default="p_value", help="Exposure: P-value column name [default: %default]"),
  make_option(c("--exp_eaf"), type="character", default="effect_allele_frequency", help="Exposure: Effect Allele Frequency column name (Recommended) [default: %default]"),
  make_option(c("--exp_n"), type="character", default="n", help="Exposure: Sample Size column name (Recommended) [default: %default]"),
  make_option(c("--exp_ncase"), type="character", default="ncases", help="Exposure: Number of Cases column name (Required for binary outcome + Steiger) [default: %default]"),
  make_option(c("--exp_type"), type="character", default="binary", help="Exposure type ('binary' or 'continuous') [default: %default]"),
  make_option(c("--exp_chr"), type="character", default="chr", help="Exposure: Chromosome column name (Optional) [default: %default]"),
  make_option(c("--exp_pos"), type="character", default="pos", help="Exposure: Position column name (Optional) [default: %default]"),
  
  # Outcome column names (Map these to the actual headers in your outcome file)
  make_option(c("--out_name"), type="character", default="Outcome", help="Name for the outcome trait"),
  make_option(c("--out_snp"), type="character", default="rs_id", help="Outcome: SNP ID column name [default: %default]"),
  make_option(c("--out_beta"), type="character", default="beta", help="Outcome: Beta column name [default: %default]"),
  make_option(c("--out_se"), type="character", default="standard_error", help="Outcome: Standard Error column name [default: %default]"),
  make_option(c("--out_ea"), type="character", default="effect_allele", help="Outcome: Effect Allele column name [default: %default]"),
  make_option(c("--out_nea"), type="character", default="other_allele", help="Outcome: Non-Effect Allele column name [default: %default]"),
  make_option(c("--out_p"), type="character", default="p_value", help="Outcome: P-value column name [default: %default]"),
  make_option(c("--out_eaf"), type="character", default="effect_allele_frequency", help="Outcome: Effect Allele Frequency column name (Recommended) [default: %default]"),
  make_option(c("--out_n"), type="character", default="n", help="Outcome: Sample Size column name (Recommended) [default: %default]"),
  make_option(c("--out_ncase"), type="character", default="ncases", help="Outcome: Number of Cases column name (Required for binary outcome + Steiger) [default: %default]"),
  make_option(c("--out_type"), type="character", default="binary", help="Outcome type ('binary' or 'continuous') [default: %default]"),
  make_option(c("--out_chr"), type="character", default="chr", help="Outcome: Chromosome column name (Optional) [default: %default]"),
  make_option(c("--out_pos"), type="character", default="pos", help="Outcome: Position column name (Optional) [default: %default]"),
  
  # Clumping parameters
  make_option(c("--clump_p"), type="numeric", default=5e-8, help="P-value threshold for IV selection [default: %default]"),
  make_option(c("--clump_kb"), type="numeric", default=10000, help="Kilobase radius for LD clumping [default: %default]"),
  make_option(c("--clump_r2"), type="numeric", default=0.001, help="R-squared threshold for LD clumping [default: %default]"),
  make_option(c("--ld_ref"), type="character", help="Path prefix to LD reference panel (PLINK bfile format, required). Example: /path/to/1kg_eur"),
  make_option(c("--plink_bin"), type="character", default=NULL, help="Path to PLINK binary (if not in PATH or auto-detected)"),
  make_option(c("--skip_clump"), action="store_true", default=FALSE, help="Skip LD clumping step (use p-value & F-stat only, ignores --ld_ref, --clump_kb, --clump_r2)"), # New Option
  
  
  # Analysis parameters
  make_option(c("--f_stat"), type="numeric", default=10, help="Minimum F-statistic for IVs [default: %default]"),
  
  # Output options
  make_option(c("--out_prefix"), type="character", help="Prefix for output files (required). Example: results/exp_vs_out"),
  make_option(c("--tmp_dir"), type="character", default="./tmp_pipeline", help="Temporary directory for intermediate files [default: %default]"),
  
  # Sensitivity analysis options
  make_option(c("--steiger"), type="logical", action="store_true", default=TRUE, help="Run Steiger filtering [default: %default]"),
  make_option(c("--no_steiger"), action="store_false", dest="steiger"),
  make_option(c("--presso"), type="logical", action="store_true", default=TRUE, help="Run MR-PRESSO outlier test [default: %default]"),
  make_option(c("--no_presso"), action="store_false", dest="presso")
)

parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

#------------------- Load Libraries & Functions -----------------#
# Moved MRPRESSO loading here, after opt is defined
if (!is.null(opt$presso) && opt$presso) {
  if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
    warning("MR-PRESSO analysis requested (--presso) but package 'MRPRESSO' is not installed. Skipping MR-PRESSO.", call. = FALSE)
    opt$presso <- FALSE # Disable PRESSO for the rest of the script
  } else {
    suppressPackageStartupMessages(library(MRPRESSO))
    message("MRPRESSO package loaded.")
  }
}

# Source helper functions from external file
functions_file <- "pipeline_functions.R"
if (!file.exists(functions_file)) {
  # Try sourcing from pipeline_scripts/ if run from workspace root
  alt_path <- file.path("pipeline_scripts", functions_file)
  if(file.exists(alt_path)) {
    functions_file <- alt_path
  } else {
    stop(sprintf("Helper functions file '%s' or '%s' not found.", functions_file, alt_path), call.=FALSE)
  }
}
message(sprintf("Sourcing helper functions from: %s", functions_file))
source(functions_file)

#------------------- Validate Arguments -----------------#

# Check required arguments
required_missing <- c()
if (is.null(opt$exp_gwas)) required_missing <- c(required_missing, "--exp_gwas")
if (is.null(opt$out_gwas)) required_missing <- c(required_missing, "--out_gwas")
if (is.null(opt$out_prefix)) required_missing <- c(required_missing, "--out_prefix")
# ld_ref is only required if we are *not* skipping clumping
if (is.null(opt$ld_ref) && !opt$skip_clump) {
  required_missing <- c(required_missing, "--ld_ref (required unless --skip_clump)")
}

if (length(required_missing) > 0) {
  print_help(parser)
  stop(sprintf("Missing required arguments: %s. See README and --help.", paste(required_missing, collapse=", ")), call.=FALSE)
}

# Check file existence early
if (!file.exists(opt$exp_gwas)) stop(sprintf("Exposure GWAS file not found: %s", opt$exp_gwas), call. = FALSE)
if (!file.exists(opt$out_gwas)) stop(sprintf("Outcome GWAS file not found: %s", opt$out_gwas), call. = FALSE)
# Check LD ref files only if clumping is not skipped
if (!opt$skip_clump && !file.exists(paste0(opt$ld_ref, ".bed"))) {
  stop(sprintf("LD reference .bed file not found based on prefix: %s (provide --skip_clump to bypass)", opt$ld_ref), call. = FALSE)
}

if (is.null(opt$exp_gwas) || is.null(opt$out_gwas) || is.null(opt$ld_ref) || is.null(opt$out_prefix)) {
  print_help(parser)
  stop("Missing required arguments: --exp_gwas, --out_gwas, --ld_ref, --out_prefix. See README and --help.", call.=FALSE)
}

# Create output and temporary directories if they don't exist
out_dir <- dirname(opt$out_prefix)
if (!dir.exists(out_dir)) {
  message("Creating output directory: ", out_dir)
  dir.create(out_dir, recursive = TRUE)
}
if (!dir.exists(opt$tmp_dir)) {
  message("Creating temporary directory: ", opt$tmp_dir)
  dir.create(opt$tmp_dir, recursive = TRUE)
}

message("---------- Pipeline Arguments ----------")
message("Exposure GWAS: ", opt$exp_gwas)
message("Outcome GWAS: ", opt$out_gwas)
message("LD Reference: ", opt$ld_ref)
message("Output Prefix: ", opt$out_prefix)
message("Outcome Type: ", opt$out_type)
message("Run Steiger: ", opt$steiger)
message("Run MR-PRESSO: ", opt$presso)
message("Skip Clumping: ", opt$skip_clump) # Add message for the new option
message("----------------------------------------")

#------------------- Helper Functions -----------------#
# Definitions are now in pipeline_functions.R, sourced above


########################################################################################
#------------------------------------- Main Pipeline ----------------------------------#
########################################################################################

message("========== Starting MR Pipeline ==========")

# 1. Load and Format Exposure Data
message("[1/6] Loading and Formatting Exposure Data...")

# Group exposure column arguments for the helper function
exp_col_args <- list(
  snp = opt$exp_snp, beta = opt$exp_beta, se = opt$exp_se,
  ea = opt$exp_ea, nea = opt$exp_nea, p = opt$exp_p,
  eaf = opt$exp_eaf, n = opt$exp_n, chr = opt$exp_chr, pos = opt$exp_pos
)

# Call function defined in pipeline_functions.R
exposure_dat <- load_and_format_gwas(
  gwas_file = opt$exp_gwas,
  type = "exposure",
  col_args = exp_col_args,
  trait_name = opt$exp_name,
  tmp_dir = opt$tmp_dir
)

# Check if exposure data frame is valid after loading
if (!inherits(exposure_dat, "data.frame") || nrow(exposure_dat) == 0) {
  stop("Exposure data loading and formatting failed or resulted in empty data frame. Exiting.", call. = FALSE)
}
message(sprintf("Successfully loaded and formatted exposure data for trait '%s'.", opt$exp_name))


# 2. Select and Clump Instruments (IVs)
message("\n[2/6] Selecting and Clumping Instruments...")

if (opt$skip_clump) {
  message("----- Skipping LD Clumping (--skip_clump specified) -----")
  message(sprintf("Selecting IVs based on p < %g and F-statistic >= %f", opt$clump_p, opt$f_stat))
  
  # 1. Filter by p-value threshold
  exposure_ivs_dat <- exposure_dat %>%
    filter(pval.exposure < opt$clump_p)
  
  if (nrow(exposure_ivs_dat) == 0) {
    stop(sprintf("No SNPs found below the significance threshold p < %g.", opt$clump_p), call. = FALSE)
  }
  message(sprintf("Found %d SNPs below p-value threshold.", nrow(exposure_ivs_dat)))
  
  # Check if N is available
  if (!"samplesize.exposure" %in% names(exposure_ivs_dat)) {
    message("Warning: Sample size column ('samplesize.exposure') not found. Cannot calculate F-statistic.")
    # Decide if we proceed without F-stat filtering or stop? Let's proceed but warn heavily.
    message("Warning: Proceeding without F-statistic filtering as sample size is missing.")
    exposure_ivs_dat$F_statistic <- NA # Add column as NA
  } else {
    # 2. Calculate F-statistic (F = beta^2 / se^2)
    exposure_ivs_dat <- exposure_ivs_dat %>%
      mutate(F_statistic = (beta.exposure^2) / (se.exposure^2))
    
    # 3. Filter by F-statistic threshold
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
  message("----- Finished Selecting Instruments (Clumping Skipped) -----")
  
} else {
  # Call function defined in pipeline_functions.R (original logic)
  message("----- Selecting and Clumping Instruments (LD Clumping Enabled) -----")
  exposure_ivs_dat <- clump_and_filter_ivs(
    exposure_dat = exposure_dat,
    clump_p = opt$clump_p,
    clump_kb = opt$clump_kb,
    clump_r2 = opt$clump_r2,
    ld_ref = opt$ld_ref,
    plink_bin = opt$plink_bin,
    min_f_stat = opt$f_stat
  )
}

# Check if IVs remain
if (!inherits(exposure_ivs_dat, "data.frame") || nrow(exposure_ivs_dat) == 0) {
  stop("No exposure instruments (IVs) remained after clumping and filtering. Exiting.", call. = FALSE)
}

# Extract the list of IV SNPs for filtering the outcome data later
exposure_ivs <- exposure_ivs_dat$SNP
message(sprintf("Final selected IVs count: %d", length(exposure_ivs)))


# 3. Load and Format Outcome Data
message("\n[3/6] Loading and Formatting Outcome Data...")

# Group outcome column arguments
out_col_args <- list(
  snp = opt$out_snp, beta = opt$out_beta, se = opt$out_se,
  ea = opt$out_ea, nea = opt$out_nea, p = opt$out_p,
  eaf = opt$out_eaf, n = opt$out_n, chr = opt$out_chr, pos = opt$out_pos, ncase = opt$out_ncase
)

# Load the full outcome GWAS first using function from pipeline_functions.R
outcome_dat_full <- load_and_format_gwas(
  gwas_file = opt$out_gwas,
  type = "outcome",
  col_args = out_col_args,
  trait_name = opt$out_name,
  tmp_dir = opt$tmp_dir
)

# Check if outcome data frame is valid after loading
if (!inherits(outcome_dat_full, "data.frame") || nrow(outcome_dat_full) == 0) {
  stop("Outcome data loading and formatting failed or resulted in empty data frame. Exiting.", call. = FALSE)
}
message(sprintf("Successfully loaded and formatted full outcome data for trait '%s'.", opt$out_name))

# Filter outcome data to keep only the selected exposure IVs
message(sprintf("Filtering outcome data for %d selected exposure IVs...", length(exposure_ivs)))
outcome_dat <- outcome_dat_full %>%
  filter(SNP %in% exposure_ivs)

# Check if any IVs were found in the outcome data
n_outcome_ivs <- nrow(outcome_dat)
n_exposure_ivs <- length(exposure_ivs)

if (n_outcome_ivs == 0) {
  message(sprintf("WARNING: None of the %d selected exposure IVs were found in the formatted outcome data for '%s'.", n_exposure_ivs, opt$out_name))
  stop("Cannot proceed without overlapping IVs. Check input files and SNP identifiers.", call.=FALSE)
} else if (n_outcome_ivs < n_exposure_ivs) {
  missing_ivs <- setdiff(exposure_ivs, outcome_dat$SNP)
  message(sprintf("WARNING: Found %d out of %d exposure IVs in the outcome data.", n_outcome_ivs, n_exposure_ivs))
  message(sprintf("Missing %d IVs (examples): %s", length(missing_ivs), paste(head(missing_ivs), collapse=", ")))
  message("Proceeding with the %d overlapping IVs.", n_outcome_ivs)
} else {
  message(sprintf("Found all %d exposure IVs in the outcome data.", n_outcome_ivs))
}


# 4. Harmonize Data
message("\n[4/6] Harmonizing Exposure and Outcome Data...")

# Ensure both inputs are data frames
if (!is.data.frame(exposure_ivs_dat) || !is.data.frame(outcome_dat)) {
  stop("Input data for harmonization are not valid data frames.", call. = FALSE)
}

# Check inputs have overlapping SNPs
common_snps <- intersect(exposure_ivs_dat$SNP, outcome_dat$SNP)
if (length(common_snps) == 0) {
  stop("No common SNPs between exposure IV data and outcome data before harmonization. Check filtering steps.", call. = FALSE)
}
message(sprintf("Harmonizing data for %d common SNPs.", length(common_snps)))

harmonized_dat <- tryCatch({
  TwoSampleMR::harmonise_data(
    exposure_dat = exposure_ivs_dat,
    outcome_dat = outcome_dat,
    action = 2 # Default action
  )
}, error = function(e) {
  stop(sprintf("Error during data harmonization: %s", e$message), call. = FALSE)
})

# Check harmonization results
n_after_harmonize <- nrow(harmonized_dat)
n_before_harmonize <- length(common_snps)

if (n_after_harmonize == 0) {
  stop("Harmonization removed all SNPs. Check allele compatibility.", call. = FALSE)
}

message(sprintf("%d SNPs remain after harmonization.", n_after_harmonize))
n_removed_harmonize = sum(!harmonized_dat$mr_keep)
if (n_removed_harmonize > 0) {
  message(sprintf("Harmonization flagged %d SNPs for removal (e.g., incompatible alleles, palindromic with intermediate AF).", n_removed_harmonize))
  harmonized_dat <- harmonized_dat %>% filter(mr_keep == TRUE)
  if (nrow(harmonized_dat) == 0) {
    stop("Harmonization removed all SNPs after explicit filtering (mr_keep=FALSE).", call. = FALSE)
  }
  message(sprintf("Proceeding with %d SNPs where mr_keep = TRUE.", nrow(harmonized_dat)))
} else {
  message("No SNPs were removed during harmonization.")
}


# 5. Run MR Analysis
message("\n[5/6] Running Mendelian Randomization Analysis...")

# Define standard MR methods to run
mr_methods_list <- c("mr_ivw", "mr_weighted_median")
if (nrow(harmonized_dat) >= 3) {
  mr_methods_list <- c(mr_methods_list, "mr_egger_regression")
} else {
  message("Fewer than 3 SNPs, excluding MR Egger regression.")
}

# Run core MR methods
message("Running core MR methods: ", paste(mr_methods_list, collapse=", "))
mr_results <- tryCatch({
  TwoSampleMR::mr(harmonized_dat, method_list = mr_methods_list)
}, error = function(e) {
  stop(sprintf("Error running core TwoSampleMR methods: %s", e$message), call.=FALSE)
})

# Add formatted ORs and CIs
mr_results <- TwoSampleMR::generate_odds_ratios(mr_results)
mr_results_dt <- data.table::as.data.table(mr_results)

message("Running sensitivity analyses (Heterogeneity, Pleiotropy)...")
# Heterogeneity test (Cochran's Q)
het_results <- tryCatch({
  data.table::as.data.table(TwoSampleMR::mr_heterogeneity(harmonized_dat))
}, error = function(e) {
  warning(sprintf("Could not run heterogeneity tests: %s", e$message), call. = FALSE)
  NULL
})

# Pleiotropy test (MR-Egger Intercept)
plt_results <- NULL
if ("MR Egger" %in% mr_results_dt$method) {
  plt_results <- tryCatch({
    data.table::as.data.table(TwoSampleMR::mr_pleiotropy_test(harmonized_dat))
  }, error = function(e) {
    warning(sprintf("Could not run MR-Egger pleiotropy test: %s", e$message), call. = FALSE)
    NULL
  })
}

# Merge sensitivity results
if (!is.null(het_results)) {
  mr_results_dt <- merge(mr_results_dt, het_results[, .(id.exposure, id.outcome, method, Q, Q_df, Q_pval)],
                         by = c("id.exposure", "id.outcome", "method"), all.x = TRUE)
}
if (!is.null(plt_results)) {
  mr_results_dt[method == "MR Egger", `:=` (
    egger_intercept = plt_results$egger_intercept,
    egger_intercept_se = plt_results$se,
    egger_intercept_pval = plt_results$pval
  )]
}

# Steiger filtering (optional)
if (opt$steiger) {
  message("Running Steiger filtering...")
  required_steiger_cols <- c("samplesize.exposure", "samplesize.outcome")
  steiger_ok <- TRUE
  if (opt$out_type == "binary" && !"ncase.outcome" %in% names(harmonized_dat)) {
    warning("Steiger filtering requested for binary outcome, but 'ncase.outcome' column not found/provided. Skipping Steiger.", call. = FALSE)
    steiger_ok <- FALSE
  } else if (!all(required_steiger_cols %in% names(harmonized_dat))) {
    missing_cols <- setdiff(required_steiger_cols, names(harmonized_dat))
    warning(sprintf("Steiger filtering requested, but required columns missing: %s. Skipping Steiger.", paste(missing_cols, collapse=", ")), call. = FALSE)
    steiger_ok <- FALSE
  }
  
  if (steiger_ok) {
    steiger_results <- tryCatch({
      TwoSampleMR::steiger_filtering(harmonized_dat)
    }, error = function(e) {
      warning(sprintf("Steiger filtering failed: %s", e$message), call. = FALSE)
      NULL
    })
    
    if (!is.null(steiger_results)) {
      n_correct_direction <- sum(steiger_results$steiger_dir, na.rm = TRUE)
      prop_correct <- round(n_correct_direction / nrow(steiger_results) * 100, 1)
      message(sprintf("Steiger directionality: %d SNPs (%f%%) suggest correct causal direction (exposure -> outcome).", n_correct_direction, prop_correct))

      # steiger_results <- as.data.table(steiger_results)  
      # mr_results_dt <- merge(mr_results_dt, steiger_results[, .(SNP, steiger_pval)], by="SNP", all.x=TRUE)
      
      steiger_filtered_data <- steiger_results %>%
        filter(steiger_dir == TRUE & steiger_pval < 0.05)
      
      if (nrow(steiger_filtered_data) > 0 && nrow(steiger_filtered_data) < nrow(harmonized_dat)) {
        message(sprintf("Running IVW on %d SNPs passing Steiger filter.", nrow(steiger_filtered_data)))
        steiger_ivw <- tryCatch({
          TwoSampleMR::mr(steiger_filtered_data, method_list = c("mr_ivw"))
        }, error = function(e) {warning(sprintf("Failed IVW on Steiger-filtered SNPs: %s", e$message), call. = FALSE); NULL})
        if (!is.null(steiger_ivw)) {
          steiger_ivw <- TwoSampleMR::generate_odds_ratios(steiger_ivw)
          steiger_ivw_dt <- data.table::as.data.table(steiger_ivw); steiger_ivw_dt$method <- "mr_ivw_steiger"
          mr_results_dt <- rbind(mr_results_dt, steiger_ivw_dt, fill = TRUE)
          message("Added 'mr_ivw_steiger' results.")
        }
      } else if (nrow(steiger_filtered_data) == nrow(harmonized_dat)) { message("All SNPs passed Steiger filtering.")
      } else { message("No SNPs passed Steiger filtering criteria.") }
    }
  }
}

# MR-PRESSO (optional)
if (opt$presso) {
  message("Running MR-PRESSO outlier test...")
  if (nrow(harmonized_dat) < 4) {
    warning("Skipping MR-PRESSO: requires at least 4 SNPs.", call. = FALSE)
  } else if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
    warning("Skipping MR-PRESSO: Package not loaded.", call.=FALSE)
  } else {
    presso_results <- tryCatch({
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = as.data.frame(harmonized_dat),
        NbDistribution = 1000,
        SignifThreshold = 0.05
      )
    }, error = function(e) { warning(sprintf("MR-PRESSO failed: %s", e$message), call. = FALSE); NULL })
    
    if (!is.null(presso_results)) {
      message("MR-PRESSO finished.")
      global_pval <- presso_results$`MR-PRESSO results`$`Global Test`$Pvalue
      message(sprintf("MR-PRESSO Global P-value: %s", global_pval))
      
      presso_main <- data.table::as.data.table(presso_results$`Main MR results`)
      presso_main[, b := `Causal Estimate`]; presso_main[, se := Sd]; presso_main[, pval := `P-value`]
      presso_main[, method := ifelse(`MR Analysis` == "Raw", "mr_presso_raw", "mr_presso_corrected")]
      presso_main[, id.exposure := harmonized_dat$id.exposure[1]]; presso_main[, id.outcome := harmonized_dat$id.outcome[1]]
      presso_main[, exposure := harmonized_dat$exposure[1]]; presso_main[, outcome := harmonized_dat$outcome[1]]
      outlier_indices <- presso_results$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
      n_outliers <- if (is.null(outlier_indices) || any(is.na(outlier_indices))) 0 else length(outlier_indices)
      presso_main[method == "mr_presso_raw", nsnp := nrow(harmonized_dat)]
      presso_main[method == "mr_presso_corrected", nsnp := nrow(harmonized_dat) - n_outliers]
      presso_main[, or := exp(b)]; presso_main[, or_lci95 := exp(b - 1.96 * se)]; presso_main[, or_uci95 := exp(b + 1.96 * se)]
      cols_to_keep <- c("id.exposure", "id.outcome", "exposure", "outcome", "method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95")
      if (!is.null(presso_results$`MR-PRESSO results`$`Distortion Test`$Pvalue)) { presso_main[, distortion_pval := presso_results$`MR-PRESSO results`$`Distortion Test`$Pvalue]; cols_to_keep <- c(cols_to_keep, "distortion_pval") }
      if (!is.null(presso_results$`MR-PRESSO results`$`Global Test`$Pvalue)) { presso_main[, presso_global_pval := presso_results$`MR-PRESSO results`$`Global Test`$Pvalue]; cols_to_keep <- c(cols_to_keep, "presso_global_pval") }
      mr_results_dt <- rbind(mr_results_dt, presso_main[, ..cols_to_keep], fill = TRUE)
      message("Added MR-PRESSO results.")
    }
  }
}


# 6. Process and Save Results
message("\n[6/6] Processing and Saving Results...")

# Calculate I-squared
mr_results_dt[!is.na(Q) & !is.na(Q_df) & Q > Q_df, I_squared := (Q - Q_df) / Q]
mr_results_dt[!is.na(Q) & !is.na(Q_df) & Q <= Q_df, I_squared := 0]

# Define output file paths
harmonized_file <- paste0(opt$out_prefix, "_harmonized_data.rds")
ivs_file <- paste0(opt$out_prefix, "_exposure_ivs.csv")
results_file <- paste0(opt$out_prefix, "_full_mr_results.csv")

# Save harmonized data
message("Saving harmonized data to: ", harmonized_file)
saveRDS(harmonized_dat, file = harmonized_file)

# Save exposure IVs
message("Saving selected exposure IVs (with stats) to: ", ivs_file)
data.table::fwrite(as.data.table(exposure_ivs_dat), ivs_file, sep = ",", na = "NA")

# Save full results table
message("Saving full MR results to: ", results_file)
mr_results_dt <- mr_results_dt[order(match(method, c("mr_ivw", "mr_ivw_steiger", "mr_egger_regression", "mr_weighted_median", "mr_presso_raw", "mr_presso_corrected")))]
data.table::fwrite(mr_results_dt, results_file, sep = ",", na = "NA")

message("========================================")
message("Pipeline Finished Successfully!")
message("Output files generated with prefix: ", opt$out_prefix)
message(" -> Harmonized data (RDS): ", harmonized_file)
message(" -> Exposure IVs (CSV):    ", ivs_file)
message(" -> Full MR results (CSV): ", results_file)
message("========================================")