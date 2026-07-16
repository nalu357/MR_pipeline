# Mendelian Randomization Pipeline Helper Functions

# Note: The main script mr_pipeline.R should load necessary libraries.
# However, including them here can be useful if functions are sourced elsewhere.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
  library(ieugwasr)
})

#' Load and Format GWAS Summary Statistics
#'
#' Reads a GWAS summary statistics file, selects and renames columns
#' according to specified parameters, performs basic cleaning, and
#' prepares the data for TwoSampleMR analysis.
#'
#' @param gwas_file Path to the GWAS summary statistics file.
#' @param type Character string, either "exposure" or "outcome".
#' @param col_args List containing column name arguments from optparse (e.g., opt$exp_snp).
#'                 Expected names in list: snp, beta, se, ea, nea, p, eaf, n, chr, pos, ncase.
#' @param trait_name Name to assign to the trait.
#' @param tmp_dir Temporary directory for fread.
#'
#' @return A data frame formatted by TwoSampleMR::format_data.
#'
load_and_format_gwas <- function(gwas_file, type = "exposure", col_args,
                                 trait_name, tmp_dir = "./tmp_pipeline",
                                 keep_snps = NULL, pval_thresh = NULL) {
  
  message(sprintf("----- Loading and Formatting %s Data -----", toupper(type)))
  message(sprintf("Reading GWAS file: %s", gwas_file))
  
  # Define mapping from standard TwoSampleMR names to user-provided column names
  # Use provided names directly from col_args list
  col_map <- list(
    SNP = col_args$snp,
    beta = col_args$beta,
    se = col_args$se,
    effect_allele = col_args$ea,
    other_allele = col_args$nea,
    pval = col_args$p
  )
  # Add optional columns only if the argument was provided (i.e., not NULL)
  if (!is.null(col_args$eaf)) col_map$eaf <- col_args$eaf
  if (!is.null(col_args$n)) col_map$samplesize <- col_args$n
  if (!is.null(col_args$chr)) col_map$chr <- col_args$chr
  if (!is.null(col_args$pos)) col_map$pos <- col_args$pos
  if (type == "outcome" && !is.null(col_args$ncase)) {
    col_map$ncase <- col_args$ncase
  }
  
  # Check if any essential columns are missing in the map (should not happen with defaults)
  essential_std_names <- c("SNP", "beta", "se", "effect_allele", "other_allele", "pval")
  if (!all(essential_std_names %in% names(col_map))) {
    stop("Internal error: Not all essential columns mapped.", call. = FALSE)
  }
  
  # Columns to select from the file are the values in col_map
  select_cols_file <- unname(unlist(col_map))
  
  # Read using data.table, selecting only necessary columns
  dt <- tryCatch({
    # Using `select =` is efficient
    data.table::fread(gwas_file, select = select_cols_file, tmpdir = tmp_dir, data.table = TRUE)
  }, error = function(e) {
    stop(sprintf("Error reading GWAS file %s: %s. Check file format and specified column names.", gwas_file, e$message), call. = FALSE)
  })
  message(sprintf("Read %d rows, %d columns from %s", nrow(dt), ncol(dt), basename(gwas_file)))
  
  # Rename columns to standard TwoSampleMR names using the inverse map
  current_names <- names(dt)
  # Create inverse map: file_name -> standard_name
  rename_map_inv <- setNames(names(col_map), unlist(col_map))
  # Select only those columns that were actually read from file for renaming
  names_to_rename <- intersect(current_names, names(rename_map_inv))
  new_names <- rename_map_inv[names_to_rename]
  data.table::setnames(dt, old = names_to_rename, new = new_names)
  message("Renamed columns to standard names (SNP, beta, se, pval, etc.)")

  # --- Early SNP subset (outcomes) ---
  # We only ever need the instrument SNPs from an outcome GWAS. Filtering here,
  # before cleaning and TwoSampleMR::format_data, avoids formatting the entire
  # file (millions of rows) - which is slow and overflows R's protection stack
  # ("C stack usage is too close to the limit") on large GWAS.
  if (!is.null(keep_snps)) {
    rows_before <- nrow(dt)
    dt <- dt[SNP %in% keep_snps]
    message(sprintf("Subset to requested instruments: %d of %d rows retained (%d unique SNPs requested).",
                    nrow(dt), rows_before, length(unique(keep_snps))))
    if (nrow(dt) == 0) {
      stop(sprintf("None of the requested instrument SNPs were found in %s (check SNP-ID column / build).",
                   basename(gwas_file)), call. = FALSE)
    }
  }

  # For an exposure we only ever instrument on SNPs passing the selection
  # threshold, so pre-filter by p-value before format_data. This is the same
  # stack-overflow guard as keep_snps, for the exposure side. (Skipped when the
  # threshold is >=1, e.g. pre-clumped QTL inputs meant to be kept in full.)
  if (!is.null(pval_thresh) && pval_thresh < 1 && "pval" %in% names(dt)) {
    rows_before <- nrow(dt)
    pv <- suppressWarnings(as.numeric(dt$pval))
    dt <- dt[!is.na(pv) & pv < pval_thresh]
    message(sprintf("Pre-filtered to pval < %g: %d of %d rows retained.",
                    pval_thresh, nrow(dt), rows_before))
    if (nrow(dt) == 0) {
      stop(sprintf("No SNPs with pval < %g in %s.", pval_thresh, basename(gwas_file)), call. = FALSE)
    }
  }

  # --- Basic Cleaning ---
  initial_rows <- nrow(dt)
  message(sprintf("Initial rows: %d", initial_rows))
  
  # 1. Allele Cleaning: Ensure alleles are uppercase and only A/T/C/G
  if ("effect_allele" %in% names(dt)) dt[, effect_allele := toupper(effect_allele)]
  if ("other_allele" %in% names(dt)) dt[, other_allele := toupper(other_allele)]
  valid_alleles <- c("A", "T", "C", "G")
  rows_before <- nrow(dt)
  dt <- dt[effect_allele %in% valid_alleles & other_allele %in% valid_alleles]
  rows_removed <- rows_before - nrow(dt)
  if (rows_removed > 0) {
    message(sprintf("Removed %d rows with invalid alleles (non A/T/C/G).", rows_removed))
  }
  
  # 2. Missing Essential Data: Remove rows with NA in core numeric/ID cols
  essential_std_names_check <- intersect(essential_std_names, names(dt))
  rows_before <- nrow(dt)
  dt <- dt[complete.cases(dt[, ..essential_std_names_check])]
  rows_removed <- rows_before - nrow(dt)
  if (rows_removed > 0) {
    message(sprintf("Removed %d rows with missing essential data (%s).", rows_removed, paste(essential_std_names_check, collapse=", ")))
  }
  
  # 3. Convert relevant columns to numeric (handle potential non-numeric values)
  numeric_cols_std <- c("beta", "se", "pval", "eaf", "samplesize", "ncase", "pos")
  numeric_cols_present <- intersect(numeric_cols_std, names(dt))
  if (length(numeric_cols_present) > 0) {
    message("Converting numeric columns...")
    for(col in numeric_cols_present) {
      # Check if already numeric (fread might do this)
      if (!is.numeric(dt[[col]])) {
        # Suppress warnings during conversion, check for NAs introduced
        dt[, (col) := suppressWarnings(as.numeric(get(col)))]
      }
    }
    # Check for NAs introduced by conversion in essential columns
    rows_before <- nrow(dt)
    dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval)] # Re-check essential numerics
    rows_removed <- rows_before - nrow(dt)
    if (rows_removed > 0) {
      message(sprintf("Removed %d rows where essential numeric columns (beta, se, pval) became NA after conversion.", rows_removed))
    }
  }
  
  # Check for empty data table before format_data
  if(nrow(dt) == 0) {
    stop(sprintf("No valid SNPs remaining for %s after cleaning. Check input file and column specifications.", type), call.=FALSE)
  }
  message(sprintf("Rows after cleaning: %d", nrow(dt)))
  
  # --- Format using TwoSampleMR ---
  # Add trait name column for format_data
  dt[, trait := trait_name]
  
  # Dynamically create arguments for format_data, passing only columns that exist in dt
  format_args <- list(
    dat = as.data.frame(dt), # format_data often prefers data.frame
    type = type,
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval",
    phenotype_col = "trait"
  )
  optional_cols_map <- c(eaf = "eaf", samplesize = "samplesize", ncase = "ncase", chr = "chr", pos = "pos")
  for (std_name in names(optional_cols_map)) {
    tsmr_arg_name <- paste0(optional_cols_map[[std_name]], "_col")
    if (std_name %in% names(dt)) {
      format_args[[tsmr_arg_name]] <- std_name
    } else {
      # If the column wasn't found or read, don't pass the argument to format_data
      message(sprintf("Optional column '%s' (mapped from user arg '%s') not found in input file or removed during cleaning. Not passing to format_data.", std_name, col_map[[std_name]]))
    }
  }
  
  formatted_data <- tryCatch({
    do.call(TwoSampleMR::format_data, format_args)
  }, error = function(e) {
    stop(sprintf("Error during TwoSampleMR::format_data for %s: %s", type, e$message), call. = FALSE)
  })
  
  message(sprintf("TwoSampleMR::format_data completed. Final data has %d SNPs.", nrow(formatted_data)))
  if(nrow(formatted_data) == 0) {
    stop(sprintf("No SNPs remaining for %s after TwoSampleMR::format_data. Check for potential issues (e.g., inconsistent data types).", type), call.=FALSE)
  }
  message(sprintf("----- Finished Loading and Formatting %s Data -----", toupper(type)))
  
  return(formatted_data)
}


#' Parse an MHC region string "CHR:START-END" into its components.
#'
#' @param mhc_region Character, e.g. "6:25000000-34000000".
#' @return A list with numeric `chr`, `start`, `end`, or NULL if unparseable.
parse_mhc_region <- function(mhc_region) {
  if (is.null(mhc_region) || is.na(mhc_region) || !nzchar(mhc_region)) return(NULL)
  m <- regmatches(mhc_region, regexec("^\\s*(\\w+)\\s*:\\s*([0-9]+)\\s*-\\s*([0-9]+)\\s*$", mhc_region))[[1]]
  if (length(m) != 4) {
    warning(sprintf("Could not parse --mhc_region '%s' (expected CHR:START-END). MHC flagging skipped.", mhc_region), call. = FALSE)
    return(NULL)
  }
  list(chr = sub("^chr", "", m[2], ignore.case = TRUE), start = as.numeric(m[3]), end = as.numeric(m[4]))
}

#' Flag (and optionally drop) instruments falling in the MHC region.
#'
#' The MHC has long-range LD and pervasive pleiotropy, so it is standard MR
#' practice to at least flag it. Requires chr/pos columns in the exposure data.
#'
#' @param dat Formatted exposure data (with chr.exposure / pos.exposure).
#' @param mhc_region Character "CHR:START-END".
#' @param exclude_mhc Logical; if TRUE, drop MHC instruments instead of flagging.
#' @return `dat` with an added logical column `mhc` (and MHC rows removed if requested).
flag_mhc_instruments <- function(dat, mhc_region, exclude_mhc = FALSE) {
  dat$mhc <- FALSE
  region <- parse_mhc_region(mhc_region)
  if (is.null(region)) return(dat)
  if (!all(c("chr.exposure", "pos.exposure") %in% names(dat))) {
    warning("chr/pos not available in exposure data (map --exp_chr/--exp_pos); MHC flagging skipped.", call. = FALSE)
    return(dat)
  }
  chr_chr <- sub("^chr", "", as.character(dat$chr.exposure), ignore.case = TRUE)
  pos_num <- suppressWarnings(as.numeric(dat$pos.exposure))
  in_mhc <- !is.na(pos_num) & chr_chr == region$chr & pos_num >= region$start & pos_num <= region$end
  dat$mhc <- in_mhc
  n_mhc <- sum(in_mhc, na.rm = TRUE)
  if (n_mhc > 0) {
    if (exclude_mhc) {
      message(sprintf("MHC: dropping %d instrument(s) in %s (--exclude_mhc).", n_mhc, mhc_region))
      dat <- dat[!in_mhc, , drop = FALSE]
    } else {
      message(sprintf("MHC: flagged %d instrument(s) in %s (kept; column 'mhc'=TRUE). Consider a sensitivity analysis excluding them.", n_mhc, mhc_region))
    }
  }
  dat
}

#' Warn about instruments whose IDs are not rsIDs.
#'
#' LD clumping and harmonisation match on rsID against the reference panel .bim,
#' so non-rsID entries (e.g. indels named "13:60994514_GT_G") are typically
#' dropped silently. This surfaces that loss up front.
#'
#' @param dat Data frame with a SNP column.
#' @return Invisibly, the vector of non-rsID SNPs.
warn_non_rsid_instruments <- function(dat) {
  non_rsid <- dat$SNP[!grepl("^rs[0-9]+$", dat$SNP)]
  if (length(non_rsid) > 0) {
    message(sprintf("NOTE: %d of %d instrument IDs are not rsIDs (e.g. %s). These usually fail rsID-based LD matching/harmonisation and may be dropped.",
                    length(non_rsid), nrow(dat),
                    paste(utils::head(non_rsid, 3), collapse = ", ")))
  }
  invisible(non_rsid)
}

#' Select, Clump, and Filter Instruments (IVs)
#'
#' Filters exposure data by p-value, performs LD clumping, calculates F-statistic,
#' and filters IVs based on F-statistic threshold.
#'
#' @param exposure_dat Formatted exposure data frame from TwoSampleMR::format_data.
#' @param clump_p P-value threshold for selecting potential IVs.
#' @param clump_kb Clumping window size (kilobases).
#' @param clump_r2 Clumping r-squared threshold.
#' @param ld_ref Path prefix to LD reference panel (PLINK format).
#' @param plink_bin Path to PLINK binary (optional).
#' @param min_f_stat Minimum F-statistic threshold for IVs.
#'
#' @return A data frame containing the selected, clumped, and filtered IVs.
#'
clump_and_filter_ivs <- function(exposure_dat, clump_p, clump_kb, clump_r2, ld_ref, plink_bin = NULL, min_f_stat = 10,
                                  mhc_region = "6:25000000-34000000", exclude_mhc = FALSE) {

  message("----- Selecting and Clumping Instruments -----")
  n_input <- nrow(exposure_dat)
  
  # Ensure required columns exist
  required_cols <- c("SNP", "pval.exposure", "beta.exposure", "se.exposure", "exposure")
  if (!all(required_cols %in% names(exposure_dat))) {
    stop("Missing required columns in exposure data for clumping (SNP, pval.exposure, beta.exposure, se.exposure, exposure). Check formatting step.", call. = FALSE)
  }
  
  # Check if N is available for context, needed later for Steiger
  if (!"samplesize.exposure" %in% names(exposure_dat)) {
    message("Warning: Sample size column ('samplesize.exposure') not found in exposure data. F-statistic calculation will proceed, but interpretation may be less reliable, and Steiger filtering will fail later if requested.")
  }
  
  # 1. Filter by p-value threshold
  significant_snps <- exposure_dat %>%
    filter(pval.exposure < clump_p)
  
  if (nrow(significant_snps) == 0) {
    stop(sprintf("No SNPs found below the significance threshold p < %g.", clump_p), call. = FALSE)
  }
  message(sprintf("Found %d SNPs below p-value threshold %g.", nrow(significant_snps), clump_p))
  warn_non_rsid_instruments(significant_snps)

  # 2. Perform LD Clumping
  # ieugwasr::ld_clump requires a data frame with 'rsid' and 'pval' columns
  clump_input_df <- significant_snps %>%
    select(rsid = SNP, pval = pval.exposure) %>%
    # ld_clump can fail if pval is exactly 0, add a small amount
    mutate(pval = ifelse(pval == 0, .Machine$double.xmin, pval)) %>%
    as.data.frame() # ld_clump expects a data.frame
  
  # Determine PLINK binary path
  plink_path <- plink_bin
  if (is.null(plink_path)) {
    # Try to find PLINK using genetics.binaRies (helper for ieugwasr)
    # Suppress potential startup messages from genetics.binaRies
    suppressPackageStartupMessages({
      if (!requireNamespace("genetics.binaRies", quietly = TRUE)) {
        warning("Package 'genetics.binaRies' not found. Relying on PLINK being in PATH or provide --plink_bin.", call. = FALSE)
      } else {
        plink_path <- tryCatch(genetics.binaRies::get_plink_binary(), error = function(e) NULL)
        if (is.null(plink_path)) {
          warning("Could not automatically find PLINK binary using 'genetics.binaRies'. Relying on PLINK being in PATH or provide --plink_bin.", call. = FALSE)
        } else {
          message("Using PLINK binary found at: ", plink_path)
        }
      }
    })
  }
  
  message(sprintf("Performing LD clumping with kb=%d, r2=%f, p=%g using LD reference: %s",
                  clump_kb, clump_r2, clump_p, ld_ref))
  
  clumped_snps_df <- tryCatch({
    ieugwasr::ld_clump(
      dat = clump_input_df,
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      clump_p = 1, # Already filtered by p-value, keep all independent signals from the input set
      bfile = ld_ref,
      plink_bin = plink_path # Pass the determined path
    )
  }, error = function(e) {
    stop(sprintf("Error during LD clumping: %s. Check PLINK installation, LD reference path (%s), and inputs.", e$message, ld_ref), call. = FALSE)
  })
  
  if (nrow(clumped_snps_df) == 0) {
    stop("Clumping removed all SNPs. Check parameters and LD reference.", call. = FALSE)
  }
  message(sprintf("Identified %d independent instruments after clumping.", nrow(clumped_snps_df)))
  
  # Filter original exposure data to keep only clumped IVs
  exposure_ivs_dat <- exposure_dat %>%
    filter(SNP %in% clumped_snps_df$rsid)
  
  # 3. Calculate F-statistic (F = beta^2 / se^2)
  exposure_ivs_dat <- exposure_ivs_dat %>%
    mutate(F_statistic = (beta.exposure^2) / (se.exposure^2))
  
  # 4. Filter by F-statistic threshold
  rows_before_f <- nrow(exposure_ivs_dat)
  exposure_ivs_dat <- exposure_ivs_dat %>%
    filter(F_statistic >= min_f_stat)
  rows_removed_f <- rows_before_f - nrow(exposure_ivs_dat)
  
  if (rows_removed_f > 0) {
    message(sprintf("Removed %d IVs with F-statistic < %f.", rows_removed_f, min_f_stat))
  }
  
  if (nrow(exposure_ivs_dat) == 0) {
    stop(sprintf("No IVs remained after F-statistic filtering (F >= %f).", min_f_stat), call.=FALSE)
  }

  # 5. Flag / optionally drop MHC instruments (long-range LD + pleiotropy)
  exposure_ivs_dat <- flag_mhc_instruments(exposure_ivs_dat, mhc_region, exclude_mhc)

  message(sprintf(
    "Instrument attrition: %d input -> %d at p<%g -> %d after clumping (r2<%g, %dkb) -> %d after F>=%g%s.",
    n_input, nrow(significant_snps), clump_p, nrow(clumped_snps_df), clump_r2, clump_kb,
    nrow(exposure_ivs_dat), min_f_stat,
    if (!exclude_mhc && "mhc" %in% names(exposure_ivs_dat)) sprintf(" (incl. %d MHC flagged)", sum(exposure_ivs_dat$mhc)) else ""))
  message(sprintf("%d IVs remain after clumping and F-statistic filtering.", nrow(exposure_ivs_dat)))
  message("----- Finished Selecting and Clumping Instruments -----")

  return(exposure_ivs_dat)
}

#' Run Mendelian Randomization Analysis for an Exposure-Outcome Pair
#'
#' Loads, harmonizes, and analyzes a single exposure-outcome pair using the appropriate MR methods
#' (Wald ratio, IVW, weighted median, MR Egger, Steiger, MR-PRESSO) depending on the number of IVs.
#' Handles all harmonization, method selection, and sensitivity analyses, and saves harmonized data and full MR results to disk.
#'
#' @param exposure_ivs_dat data.table. The clumped and formatted exposure IVs.
#' @param outcome_file character. Path to the outcome GWAS file.
#' @param outcome_name character. Name of the outcome trait.
#' @param out_col_args list. Column mapping for the outcome GWAS.
#' @param out_prefix character. Prefix for output files.
#' @param opt list. List of pipeline options (parsed arguments).
#'
#' @return data.table of MR results for this exposure-outcome pair (invisible, but also written to disk).
run_mr_analysis <- function(exposure_ivs_dat, outcome_file, outcome_name, out_col_args, out_prefix, opt) {
  message(sprintf("Processing outcome: %s", outcome_file))
  outcome_dat <- load_and_format_gwas(
    gwas_file = outcome_file,
    type = "outcome",
    col_args = out_col_args,
    trait_name = outcome_name,
    tmp_dir = opt$tmp_dir,
    keep_snps = exposure_ivs_dat$SNP
  )
  if (!inherits(outcome_dat, "data.frame") || nrow(outcome_dat) == 0) {
    warning(sprintf("Outcome data loading failed for %s. Skipping.", outcome_file))
    return(NULL)
  }
  outcome_dat <- outcome_dat[outcome_dat$SNP %in% exposure_ivs_dat$SNP, ]
  if (nrow(outcome_dat) == 0) {
    warning(sprintf("No overlapping IVs found in outcome %s. Skipping.", outcome_file))
    return(NULL)
  }
  harmonized_dat <- tryCatch({
    TwoSampleMR::harmonise_data(
      exposure_dat = exposure_ivs_dat,
      outcome_dat = outcome_dat,
      action = 2
    )
  }, error = function(e) {
    warning(sprintf("Error during harmonization for %s: %s", outcome_file, e$message), call. = FALSE)
    NULL
  })
  if (is.null(harmonized_dat) || nrow(harmonized_dat) == 0) {
    warning(sprintf("Harmonization failed or removed all SNPs for %s. Skipping.", outcome_file))
    return(NULL)
  }
  
  n_snps <- nrow(harmonized_dat)
  if (n_snps == 1) {
    mr_methods_list <- c("mr_wald_ratio")
  } else if (n_snps == 2) {
    mr_methods_list <- c("mr_ivw", "mr_weighted_median")
  } else if (n_snps >= 3) {
    mr_methods_list <- c("mr_ivw", "mr_weighted_median", "mr_egger_regression")
  }
  
  mr_results <- tryCatch({
    TwoSampleMR::mr(harmonized_dat, method_list = mr_methods_list)
  }, error = function(e) {
    warning(sprintf("Error running core TwoSampleMR methods: %s", e$message), call.=FALSE)
    NULL
  })
  if (is.null(mr_results)) return(NULL)
  mr_results <- TwoSampleMR::generate_odds_ratios(mr_results)
  mr_results_dt <- data.table::as.data.table(mr_results)
  
  if (n_snps >= 3) {
    het_results <- tryCatch({
      data.table::as.data.table(TwoSampleMR::mr_heterogeneity(harmonized_dat))
    }, error = function(e) NULL)
    if (!is.null(het_results)) {
      mr_results_dt <- merge(mr_results_dt, het_results[, .(id.exposure, id.outcome, method, Q, Q_df, Q_pval)],
                             by = c("id.exposure", "id.outcome", "method"), all.x = TRUE)
    }
    
    plt_results <- NULL
    if ("MR Egger" %in% mr_results_dt$method) {
      plt_results <- tryCatch({
        data.table::as.data.table(TwoSampleMR::mr_pleiotropy_test(harmonized_dat))
      }, error = function(e) NULL)
    }
    if (!is.null(plt_results)) {
      mr_results_dt[method == "MR Egger", `:=` (
        egger_intercept = plt_results$egger_intercept,
        egger_intercept_se = plt_results$se,
        egger_intercept_pval = plt_results$pval
      )]
    }
  }
  
  if (opt$steiger) {
    steiger_results <- tryCatch({
      TwoSampleMR::steiger_filtering(harmonized_dat)
    }, error = function(e) NULL)
    if (!is.null(steiger_results)) {
      # steiger_results <- as.data.table(steiger_results)
      # mr_results_dt <- merge(mr_results_dt, steiger_results[, .(SNP, steiger_dir)], by = "SNP", all.x = TRUE)
      steiger_filtered_data <- steiger_results %>%
        filter(steiger_dir == TRUE & steiger_pval < 0.05)
      if (nrow(steiger_filtered_data) > 0 && nrow(steiger_filtered_data) < nrow(harmonized_dat)) {
        steiger_ivw <- tryCatch({
          TwoSampleMR::mr(steiger_filtered_data, method_list = c("mr_ivw"))
        }, error = function(e) NULL)
        if (!is.null(steiger_ivw)) {
          steiger_ivw <- TwoSampleMR::generate_odds_ratios(steiger_ivw)
          steiger_ivw_dt <- data.table::as.data.table(steiger_ivw); steiger_ivw_dt$method <- "mr_ivw_steiger"
          mr_results_dt <- rbind(mr_results_dt, steiger_ivw_dt, fill = TRUE)
        }
      }
    }
  }
  
  if (opt$presso && n_snps >= 4 && requireNamespace("MRPRESSO", quietly = TRUE)) {
    presso_results <- tryCatch({
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome", SdExposure = "se.exposure",
        OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
        data = as.data.frame(harmonized_dat),
        NbDistribution = 1000, SignifThreshold = 0.05
      )
    }, error = function(e) NULL)
    if (!is.null(presso_results)) {
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
    }
  }
  
  harmonized_file <- paste0(out_prefix, "_harmonized_data.rds")
  saveRDS(harmonized_dat, file = harmonized_file)
  
  results_file <- paste0(out_prefix, "_full_mr_results.csv")
  data.table::fwrite(mr_results_dt, results_file, sep = ",", na = "NA")
  message(sprintf("Saved results to: %s", results_file))
  return(mr_results_dt)
}

#' Process and Annotate Mendelian Randomization Results Across All Exposure-Outcome Pairs
#'
#' This function takes a combined data.table of MR results for all exposure-outcome pairs,
#' and performs post-processing to:
#' - Subset and combine IVW and Wald ratio results for each pair
#' - Flag pleiotropy (significant MR Egger intercept), heterogeneity (I2 > 0.5), and Steiger directionality discordance
#' - Compute direction concordance across sensitivity methods (MR Egger, Weighted Median, Steiger, MR-PRESSO)
#' - Add FDR-adjusted p-values
#' - Save a prettified, annotated summary table to disk
#'
#' @param all_mr_results data.table. Combined MR results for all exposure-outcome pairs.
#'
#' @return data.table. Prettified and annotated MR results, one row per exposure-outcome-method, and writes the summary to disk.
process_mr_results <- function(all_mr_results, opt) {
  all_mr_results <- as.data.table(all_mr_results)
  all_mr_results[, id.exposure:=NULL]
  all_mr_results[, id.outcome:=NULL]
  all_mr_results[, lo_ci:=NULL]
  all_mr_results[, up_ci:=NULL]	
  all_mr_results[, egger_intercept:=NULL]	
  all_mr_results[, egger_intercept_se:=NULL]	
  
  #for(i in nrow(all_mr_results)){
    cur_mr_result <- all_mr_results#[i]
    # Subset IVW and Wald ratio results
    cur_mr_result.IVW <- cur_mr_result[method == "Inverse variance weighted" & nsnp > 2]
    cur_mr_result.wald <- cur_mr_result[method == "Wald ratio"]
    
    # Prepare output
    all_mr_results.all <- data.table()
    
    # Loop over unique outcome/ancestry combinations
    for (out in unique(cur_mr_result$outcome)) {
     for (exp in unique(cur_mr_result[outcome==out, exposure])) {
        # MR Egger intercept FDR
        cur_mr_result.MREgger <- cur_mr_result[outcome == out & exposure == exp & method == "MR Egger"]
        # all_mr_results.MREgger[, egger_intercept_pval_fdr := p.adjust(egger_intercept_pval, method = "fdr")]
        MREgger.sig <- cur_mr_result.MREgger[egger_intercept_pval < 0.05]
        cur_mr_result.IVW$egger_intercept_pval <- cur_mr_result.MREgger$egger_intercept_pval
        
        # Heterogeneity flag
        #Compute the I2 statistics
        cur_mr_result.IVW$I2 <- (cur_mr_result.IVW$Q-cur_mr_result.IVW$Q_df)/cur_mr_result.IVW$Q
        Het.sig <- cur_mr_result.IVW[outcome == out & exposure == exp & I2 > 0.5]
        
        # Steiger directionality difference
        if("mr_ivw_steiger" %in% unique(cur_mr_result$method)){
          steiger_methods <- c("mr_ivw_steiger", "Inverse variance weighted")
          steiger_dt <- cur_mr_result[outcome == out & exposure == exp & method %in% steiger_methods, .(b)]
          Steiger.diff <- length(unique(sign(steiger_dt$b))) > 1
        }
        else Steiger.diff <- NA
        
        # Precompute sensitivity betas for all clusters
        get_beta <- function(dt, m) if (nrow(dt[method == m]) > 0) dt[method == m, b] else NA
        
        tmp.res <- cur_mr_result[outcome == out & exposure == exp]
        if (nrow(tmp.res) == 0) next
        
        # Add presso results
        if(any(grepl(cur_mr_result$method, pattern = "mr_presso"))){
          cur_mr_result.IVW$presso_global_pval <- cur_mr_result[method=="mr_presso_raw", presso_global_pval]
          cur_mr_result.IVW$presso_distortion_pval <- cur_mr_result[method=="mr_presso_raw", distortion_pval]
        }
        
        # Sensitivity betas
        if(any(grepl(cur_mr_result$method, pattern = "mr_presso")) & "mr_ivw_steiger" %in% unique(cur_mr_result$method)){
          beta.sensitivity <- list(
            MRPRESSO = if ("mr_presso_corrected" %in% tmp.res$method) {
              dval <- tmp.res[method == "mr_presso_corrected", distortion_pval]
              dval <- ifelse(dval == "<0.001", 1e-4, as.numeric(dval))
              if (!is.na(dval) && dval < 0.05) tmp.res[method == "mr_presso_corrected", b]
              else tmp.res[method == "mr_presso_raw", b]
            } else NA,
            MREgger = get_beta(tmp.res, "MR Egger"),
            WeightedMedian = get_beta(tmp.res, "Weighted median"),
            Steiger = get_beta(tmp.res, "mr_ivw_steiger")
          )
        }
        else if(!any(grepl(cur_mr_result$method, pattern = "mr_presso")) & "mr_ivw_steiger" %in% unique(cur_mr_result$method)){
          beta.sensitivity <- list(
            MREgger = get_beta(tmp.res, "MR Egger"),
            WeightedMedian = get_beta(tmp.res, "Weighted median"),
            Steiger = get_beta(tmp.res, "mr_ivw_steiger")
          )
        }
        else if(any(grepl(cur_mr_result$method, pattern = "mr_presso")) & !("mr_ivw_steiger" %in% unique(cur_mr_result$method))){
          beta.sensitivity <- list(
            MRPRESSO = if ("mr_presso_corrected" %in% tmp.res$method) {
              dval <- tmp.res[method == "mr_presso_corrected", distortion_pval]
              dval <- ifelse(dval == "<0.001", 1e-4, as.numeric(dval))
              if (!is.na(dval) && dval < 0.05) tmp.res[method == "mr_presso_corrected", b]
              else tmp.res[method == "mr_presso_raw", b]
            } else NA,
            MREgger = get_beta(tmp.res, "MR Egger"),
            WeightedMedian = get_beta(tmp.res, "Weighted median")
          )
        } else {
          beta.sensitivity <- list(
            MREgger = get_beta(tmp.res, "MR Egger"),
            WeightedMedian = get_beta(tmp.res, "Weighted median")
          )
        }
          
        # Direction concordance
        ivw_beta <- cur_mr_result.IVW[outcome == out & exposure == exp, b]
        Prop.SameDir <- sum(sign(unlist(beta.sensitivity)) != sign(ivw_beta), na.rm = TRUE)
        cur_mr_result.IVW[outcome == out & exposure == exp, DiffDirection := Prop.SameDir != 0]
        
        # Combine IVW and Wald for this outcome/pop
        cur_mr_result.out <- rbind(
          cur_mr_result.IVW[outcome == out & exposure == exp],
          cur_mr_result.wald[outcome == out & exposure == exp],
          fill = TRUE
        )
        
        # Add flags
        cur_mr_result.out[, FlagPleiotropy := nrow(MREgger.sig) > 0]
        cur_mr_result.out[, FlagHeterogeneity := nrow(Het.sig) > 0]
        cur_mr_result.out[, FlagSteiger := Steiger.diff]
        
        all_mr_results.all <- rbind(all_mr_results.all, cur_mr_result.out, fill = TRUE)
      }
    }
  #}
  
  # FDR-adjusted p-value
  all_mr_results.all[, distortion_pval:=NULL]
  all_mr_results.all[, p_value_fdr := p.adjust(pval, method = "fdr")]
  # all_mr_results.all[, p_value_plot := ifelse(p_value_fdr < 0.05 & DiffDirection == FALSE & FlagHeterogeneity == FALSE & FlagPleiotropy == FALSE, 0, 2)]
  results_file <- paste0(opt$out_prefix, "all_processed_mr_results.csv")
  data.table::fwrite(all_mr_results.all, results_file, sep = ",", na = "NA")
  
  return(all_mr_results.all)
}

