# Mendelian Randomization Pipeline Helper Functions

# Note: The main script mr_pipeline.R should load necessary libraries.
# However, including them here can be useful if functions are sourced elsewhere.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(TwoSampleMR)
  library(ieugwasr)
})

# Persistent cache helpers.  Keys include the source path, size and mtime plus
# every transformation input, so an edited source file or changed CLI mapping
# naturally misses the old cache.  This deliberately avoids hashing the entire
# (often multi-GB) source file just to decide whether to reuse a cache.
cache_token <- function(...) {
  x <- paste(vapply(list(...), function(z) paste(z, collapse = "\034"), character(1)), collapse = "\035")
  # A wide (64-bit) key: a collision would silently return the WRONG cached
  # object (e.g. another pair's harmonised data), so we avoid the previous
  # ~31-bit rolling hash. Prefer digest::xxhash64; fall back to a
  # dependency-free double rolling hash (two primes) when digest is absent.
  if (requireNamespace("digest", quietly = TRUE))
    return(digest::digest(x, algo = "xxhash64"))
  bytes <- utf8ToInt(enc2utf8(x)); h1 <- 0; h2 <- 0
  for (b in bytes) {
    h1 <- (h1 * 31  + b) %% 2147483647   # 2^31 - 1 (prime)
    h2 <- (h2 * 129 + b) %% 2147483629   # largest prime < 2^31
  }
  sprintf("%08x%08x", as.integer(h1), as.integer(h2))
}

source_fingerprint <- function(path) {
  info <- file.info(path)
  if (is.na(info$size)) stop(sprintf("Input file does not exist or cannot be read: %s", path), call. = FALSE)
  c(normalizePath(path, mustWork = FALSE), as.character(info$size),
    format(info$mtime, tz = "UTC", usetz = TRUE))
}

cache_path <- function(cache_dir, namespace, key) {
  if (is.null(cache_dir) || !nzchar(cache_dir)) return(NULL)
  dir <- file.path(cache_dir, namespace)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(key, ".rds"))
}

cache_read <- function(path, label) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  value <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.null(value)) {
    # Bump mtime so prune_cache() evicts by true recency-of-use (LRU), not just
    # recency-of-write: entries reused across a grid survive eviction.
    tryCatch(Sys.setFileTime(path, Sys.time()), error = function(e) NULL)
    message(sprintf("Cache hit (%s): %s", label, path))
  }
  value
}

cache_write <- function(value, path) {
  if (is.null(path)) return(invisible(FALSE))
  tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
  saveRDS(value, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    warning(sprintf("Could not write cache file: %s", path), call. = FALSE)
    return(invisible(FALSE))
  }
  invisible(TRUE)
}

#' Bound the on-disk cache with LRU eviction.
#'
#' The caches (gwas_subsets/, proxy_substitutions/, harmonized/) otherwise grow
#' without limit. When the total exceeds `max_gb`, the least-recently-USED files
#' (oldest mtime; cache_read() refreshes mtime on a hit) are deleted until the
#' cache is back under the cap. Call once at start-up, before the run reads or
#' writes anything, so a run never evicts the entries it just created.
#'
#' @param cache_dir cache directory (NULL/"" -> no-op).
#' @param max_gb size cap in GB; NULL/NA/<=0 means unlimited (no eviction).
prune_cache <- function(cache_dir, max_gb = 20) {
  if (is.null(cache_dir) || !nzchar(cache_dir) || !dir.exists(cache_dir)) return(invisible(FALSE))
  if (is.null(max_gb) || is.na(max_gb) || max_gb <= 0) return(invisible(FALSE))
  files <- list.files(cache_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(invisible(FALSE))
  info <- file.info(files)
  info <- info[!is.na(info$size), , drop = FALSE]
  total <- sum(info$size); cap <- max_gb * 1024^3
  if (total <= cap) return(invisible(FALSE))
  ord <- order(info$mtime)                         # oldest (least recently used) first
  freed <- 0; removed <- 0L; paths <- rownames(info)
  for (i in ord) {
    if (total - freed <= cap) break
    if (isTRUE(file.remove(paths[i]))) { freed <- freed + info$size[i]; removed <- removed + 1L }
  }
  message(sprintf("Cache prune: removed %d least-recently-used file(s) (%.2f GB) to keep '%s' under %g GB.",
                  removed, freed / 1024^3, cache_dir, max_gb))
  invisible(TRUE)
}

#' Read and Clean GWAS Summary Statistics
#'
#' Reads a GWAS file, selects/renames columns to standard names, optionally
#' subsets early (by SNP list and/or p-value threshold) to keep processing
#' cheap on very large files, and performs basic cleaning. It does NOT run
#' TwoSampleMR::format_data - that is done separately by format_gwas() on a
#' small selected set, which avoids formatting millions of rows (slow, and it
#' overflows R's protection stack).
#'
#' @param gwas_file Path to the GWAS summary statistics file.
#' @param type Character string, either "exposure" or "outcome".
#' @param col_args List containing column name arguments from optparse (e.g., opt$exp_snp).
#'                 Expected names in list: snp, beta, se, ea, nea, p, eaf, n, chr, pos, ncase.
#' @param trait_name Name to assign to the trait (added as a `trait` column).
#' @param tmp_dir Temporary directory for fread.
#' @param keep_snps Optional character vector; keep only these SNPs (used for outcomes).
#' @param pval_thresh Optional numeric; keep only rows with pval < this (early shrink for exposures).
#' @param n_total Optional numeric; constant total sample size assigned to every
#'                SNP when the file has no per-SNP N column (enables Steiger).
#'
#' @return A cleaned data.table with standard column names (and a `trait` column).
#'
read_gwas <- function(gwas_file, type = "exposure", col_args,
                      trait_name = NULL, tmp_dir = "./tmp_pipeline",
                      keep_snps = NULL, pval_thresh = NULL, n_total = NULL, ncase_total = NULL,
                      cache_dir = NULL, iv_keep = NULL) {

  gwas_cache_key <- cache_token(source_fingerprint(gwas_file), type, unlist(col_args, use.names = TRUE),
                                trait_name, sort(unique(keep_snps)), pval_thresh, n_total, ncase_total,
                                sort(unique(iv_keep)))
  gwas_cache_file <- cache_path(cache_dir, "gwas_subsets", gwas_cache_key)
  cached <- cache_read(gwas_cache_file, "cleaned GWAS subset")
  if (!is.null(cached)) return(data.table::as.data.table(cached))
  
  message(sprintf("----- Reading %s Data -----", toupper(type)))
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
  if (!is.null(col_args$ncase)) {   # read ncase for exposure or outcome (binary Steiger)
    col_map$ncase <- col_args$ncase
  }
  
  # Check if any essential columns are missing in the map (should not happen with defaults)
  essential_std_names <- c("SNP", "beta", "se", "effect_allele", "other_allele", "pval")
  if (!all(essential_std_names %in% names(col_map))) {
    stop("Internal error: Not all essential columns mapped.", call. = FALSE)
  }
  
  # Only select columns that actually exist in the file. Reading the header
  # first avoids noisy "Column name '...' not found" warnings for optional
  # columns the file doesn't have (e.g. a default n/ncase lookup on a file that
  # uses --*_n_total), and gives a clear error if an ESSENTIAL column is missing.
  header_cols <- tryCatch(names(data.table::fread(gwas_file, nrows = 0)),
                          error = function(e) stop(sprintf("Could not read header of %s: %s", basename(gwas_file), e$message), call. = FALSE))
  requested <- unname(unlist(col_map))
  essential_file_cols <- unname(unlist(col_map[essential_std_names]))
  ess_missing <- setdiff(essential_file_cols, header_cols)
  if (length(ess_missing) > 0) {
    stop(sprintf("Essential column(s) not found in %s: %s. Check the %s column-name mappings (snp/beta/se/ea/nea/p).",
                 basename(gwas_file), paste(ess_missing, collapse = ", "), type), call. = FALSE)
  }
  opt_missing <- setdiff(setdiff(requested, header_cols), essential_file_cols)
  if (length(opt_missing) > 0) {
    message(sprintf("Optional column(s) not in %s, skipped: %s.", basename(gwas_file), paste(opt_missing, collapse = ", ")))
  }
  select_cols_file <- intersect(requested, header_cols)

  # Read using data.table, selecting only the columns that exist
  dt <- tryCatch({
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
    keep <- !is.na(pv) & pv < pval_thresh
    # A pre-defined instrument list (--exp_ivs) is authoritative: never drop a
    # listed SNP just because it is sub-threshold in THIS GWAS. Keep it alongside
    # the genome-wide-significant pool (the latter still serves LD-proxy search).
    n_iv_sub <- 0L
    if (!is.null(iv_keep)) {
      in_list <- dt$SNP %in% iv_keep
      n_iv_sub <- sum(in_list & !keep)
      keep <- keep | in_list
    }
    dt <- dt[keep]
    message(sprintf("Pre-filtered to pval < %g: %d of %d rows retained%s.",
                    pval_thresh, nrow(dt), rows_before,
                    if (n_iv_sub > 0) sprintf(" (incl. %d sub-threshold --exp_ivs SNPs kept)", n_iv_sub) else ""))
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
  # Per-column NA counts, so we can name the culprit if rows get dropped (e.g. an
  # essential column mapped to a column that exists but is empty/NA in the file).
  na_by_col <- vapply(essential_std_names_check,
                      function(cc) sum(is.na(dt[[cc]]) | (is.character(dt[[cc]]) & dt[[cc]] == "")),
                      numeric(1))
  dt <- dt[complete.cases(dt[, ..essential_std_names_check])]
  rows_removed <- rows_before - nrow(dt)
  if (rows_removed > 0) {
    culprits <- na_by_col[na_by_col > 0]
    culprit_str <- if (length(culprits)) paste(sprintf("%s(mapped from '%s')=%d NA",
                       names(culprits), unlist(col_map[names(culprits)]), culprits), collapse = ", ") else "unknown"
    message(sprintf("Removed %d of %d rows with missing essential data. NA counts by column: %s.",
                    rows_removed, rows_before, culprit_str))
    if (nrow(dt) == 0) {
      fully_na <- names(na_by_col[na_by_col == rows_before])
      if (length(fully_na)) {
        stop(sprintf("Essential column(s) %s are entirely empty/NA in %s (mapped from %s). Check the %s p/beta/se column mappings - e.g. an INF vs non-INF p-value column.",
                     paste(fully_na, collapse = ", "), basename(gwas_file),
                     paste(sprintf("'%s'", unlist(col_map[fully_na])), collapse = ", "), type), call. = FALSE)
      }
    }
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
  
  # Constant total-N fallback: if the file had no per-SNP sample-size column,
  # assign the supplied total N to every SNP so F-statistics and Steiger have a
  # sample size to work with.
  if (!is.null(n_total)) {
    if (!"samplesize" %in% names(dt)) {
      dt[, samplesize := as.numeric(n_total)]
      message(sprintf("Assigned constant sample size N = %s to all SNPs (no per-SNP N column found).",
                      format(n_total, scientific = FALSE)))
    } else {
      message("Per-SNP sample-size column present; ignoring supplied total N.")
    }
  }

  # Constant total case-count fallback (binary traits): if there is no per-SNP
  # ncase column, assign the supplied total number of cases to every SNP so
  # Steiger's log-odds r2 (--steiger_logodds) has case counts to work with.
  if (!is.null(ncase_total)) {
    if (!"ncase" %in% names(dt)) {
      dt[, ncase := as.numeric(ncase_total)]
      message(sprintf("Assigned constant case count Ncase = %s to all SNPs (no per-SNP ncase column found).",
                      format(ncase_total, scientific = FALSE)))
    } else {
      message("Per-SNP ncase column present; ignoring supplied total case count.")
    }
  }

  if (!is.null(trait_name)) dt[, trait := trait_name]
  message(sprintf("----- Finished Reading %s Data -----", toupper(type)))
  result <- dt[]
  cache_write(result, gwas_cache_file)
  return(result)
}

#' Format a cleaned GWAS table with TwoSampleMR::format_data
#'
#' Runs format_data on an ALREADY-SMALL, cleaned table (selected instruments,
#' or an outcome subset to the instruments). Optional columns that are present
#' are passed through automatically.
#'
#' @param dt Cleaned data.table from read_gwas().
#' @param type "exposure" or "outcome".
#' @param trait_name Optional; sets/overwrites the `trait` column.
#' @return A data.frame from TwoSampleMR::format_data.
format_gwas <- function(dt, type = "exposure", trait_name = NULL) {
  dt <- data.table::as.data.table(dt)
  if (!is.null(trait_name)) dt[, trait := trait_name]
  if (!"trait" %in% names(dt)) dt[, trait := type]
  
  format_args <- list(
    dat = as.data.frame(dt), type = type,
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "effect_allele", other_allele_col = "other_allele",
    pval_col = "pval", phenotype_col = "trait"
  )
  optional_cols_map <- c(eaf = "eaf", samplesize = "samplesize", ncase = "ncase", chr = "chr", pos = "pos")
  for (std_name in names(optional_cols_map)) {
    if (std_name %in% names(dt)) {
      format_args[[paste0(optional_cols_map[[std_name]], "_col")]] <- std_name
    }
  }
  
  formatted <- tryCatch(
    do.call(TwoSampleMR::format_data, format_args),
    error = function(e) stop(sprintf("Error during TwoSampleMR::format_data for %s: %s", type, e$message), call. = FALSE)
  )
  if (nrow(formatted) == 0) {
    stop(sprintf("No SNPs remaining for %s after TwoSampleMR::format_data.", type), call. = FALSE)
  }
  message(sprintf("format_data completed for %s: %d SNPs.", type, nrow(formatted)))
  formatted
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
    n_blank <- sum(!nzchar(non_rsid))
    examples <- utils::head(unique(non_rsid[nzchar(non_rsid)]), 3)
    example_str <- if (length(examples) > 0) sprintf("e.g. %s", paste(examples, collapse = ", ")) else "all blank IDs"
    blank_str <- if (n_blank > 0) sprintf(" (%d are blank)", n_blank) else ""
    message(sprintf("NOTE: %d of %d instrument IDs are not rsIDs (%s)%s. These usually fail rsID-based LD matching/harmonisation and may be dropped.",
                    length(non_rsid), nrow(dat), example_str, blank_str))
  }
  invisible(non_rsid)
}

#' Read a list of pre-defined instrument SNP IDs from a file.
#'
#' Accepts a one-per-line list (no header) or a table with a SNP-like column
#' (SNP/rsid/rs_id/variant/variant_id/markername, case-insensitive; else the
#' first column).
#'
#' @param file Path to the instrument-list file.
#' @param col Optional explicit column name holding the SNP IDs; if NULL, a
#'            SNP-like column is auto-detected (else the first column is used).
#' @return Character vector of unique, non-empty SNP IDs.
read_iv_list <- function(file, col = NULL) {
  if (grepl("\\.xlsx?$", file, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop(sprintf("Instrument list '%s' is an Excel file but package 'readxl' is not installed.", basename(file)), call. = FALSE)
    dt <- data.table::as.data.table(readxl::read_excel(file, col_types = "text"))
  } else {
    dt <- data.table::fread(file, header = "auto")
  }
  if (!is.null(col)) {
    if (!col %in% names(dt)) {
      stop(sprintf("--exp_ivs_col '%s' not found in %s (columns: %s).",
                   col, basename(file), paste(names(dt), collapse = ", ")), call. = FALSE)
    }
    snp_col <- col
  } else {
    # Pick the column that actually holds SNP IDs. Excel exports often have many
    # blank/duplicated headers (readxl renames them ...2, ...3, ...), so choosing
    # the FIRST column is unreliable. Score every column by how many values look
    # like rsIDs (then chr:pos); fall back to a SNP-named column, then column 1.
    rs_score <- vapply(dt, function(v) sum(grepl("^rs[0-9]+$", trimws(as.character(v)))), integer(1))
    cp_score <- vapply(dt, function(v) sum(grepl("^([0-9]{1,2}|[XYxy])[:_][0-9]+", trimws(as.character(v)))), integer(1))
    if (max(rs_score) > 0) {
      snp_col <- names(dt)[which.max(rs_score)]
    } else if (max(cp_score) > 0) {
      snp_col <- names(dt)[which.max(cp_score)]
    } else {
      snp_col <- names(dt)[tolower(names(dt)) %in%
                             c("snp", "rsid", "rs_id", "variant", "variant_id", "markername")][1]
      if (is.na(snp_col)) snp_col <- names(dt)[1]
    }
    message(sprintf("read_iv_list: '%s' -> column '%s' (%d rsID-like values of %d rows).",
                    basename(file), snp_col, max(rs_score), nrow(dt)))
  }
  ids <- as.character(dt[[snp_col]])
  # Headerless single-column list: the "column name" is really the first ID.
  if (ncol(dt) == 1 && grepl("^(rs[0-9]+|[0-9]+:[0-9]+)", snp_col)) ids <- c(snp_col, ids)
  ids <- unique(ids)
  ids[!is.na(ids) & nzchar(ids)]
}

#' Read one or more trait "info files" (manifests) into a list of trait specs.
#'
#' A manifest is a CSV/TSV or Excel (.xlsx) table with one row per trait
#' describing where its GWAS lives and which columns to use, so the grid driver
#' (mr_grid.R) does not need per-trait command-line column flags. The expected
#' headers (case-insensitive, dots/spaces/underscores interchangeable) are:
#'   Short name, file.name, N, Ncases, Population prevalence, Genome build,
#'   ivs.file, Downloaded,
#'   chr.col, pos.col, snp.col, ea.col, nea.col, eaf.col, beta.col, se.col, pval.col
#'
#' Behaviour:
#'   * Numbers tolerate thousands separators and stray whitespace ("588,452 ").
#'   * `N` / `Ncases` may be either a NUMBER (applied as a constant sample size /
#'     case count) OR the NAME of a per-SNP column in the GWAS - if the value is
#'     non-numeric text it is taken as a column name.
#'   * A trait is BINARY when Ncases is populated, else CONTINUOUS.
#'   * Empty eaf.col / chr.col / pos.col are fine (dropped).
#'   * file.name / ivs.file: absolute paths used as-is, else resolved vs data_dir.
#'   * Rows with Downloaded == "no" or an empty file.name are skipped (warning).
#'   * A row missing a required per-row value (snp/ea/nea/beta/se/pval col) is
#'     kept but WARNED about - its MR runs will fail rather than silently mislead.
#'   * Short names are sanitised for safe filenames (e.g. "FEV1/FVC" -> "FEV1_FVC").
#'   * Genome build maps to an MHC region (b37 default, b38 lifted coordinates).
#'
#' @param paths character vector of manifest paths (.csv/.tsv/.txt/.xlsx, >= 1).
#' @param data_dir optional base directory to resolve relative file.name/ivs.file.
#' @param select optional character vector of Short names to keep (NULL = all).
#'   Requesting a name absent from the manifests is an error.
#' @param role "exposure"/"outcome", used only for clearer messages.
#' @return named list (keyed by sanitised Short name) of trait specs.
read_trait_manifest <- function(paths, data_dir = NULL, select = NULL, role = "trait") {
  norm <- function(x) gsub("[ ._]+", "", tolower(x))
  # Strip commas AND all whitespace (incl. non-breaking space U+00A0) before parsing.
  as_num <- function(x) {
    x <- gsub("[[:space:] ,]+", "", trimws(as.character(x)))
    suppressWarnings(ifelse(nzchar(x), as.numeric(x), NA_real_))
  }
  blank <- function(x) is.null(x) || is.na(x) || !nzchar(trimws(as.character(x)))
  # Filenames must not contain path separators or other awkward characters.
  safe_name <- function(x) gsub("[^A-Za-z0-9._+-]+", "_", trimws(x))
  resolve <- function(f) {
    if (blank(f)) return(NA_character_)
    f <- trimws(as.character(f))
    if (grepl("^(/|~|[A-Za-z]:)", f) || is.null(data_dir)) f else file.path(data_dir, f)
  }
  mhc_for_build <- function(build) {
    if (grepl("38", norm(build))) "6:25726063-33400644" else "6:25000000-34000000"
  }
  # A manifest N/Ncases cell is either a constant (number) or a column name (text).
  # Returns list(total = <numeric or NULL>, col = <colname or NULL>).
  n_field <- function(v) {
    if (blank(v)) return(list(total = NULL, col = NULL))
    num <- as_num(v)
    if (is.na(num)) list(total = NULL, col = trimws(as.character(v)))
    else list(total = num, col = NULL)
  }
  read_table <- function(p) {
    if (grepl("\\.xlsx?$", p, ignore.case = TRUE)) {
      if (!requireNamespace("readxl", quietly = TRUE))
        stop(sprintf("Manifest '%s' is an Excel file but package 'readxl' is not installed.", basename(p)), call. = FALSE)
      data.table::as.data.table(readxl::read_excel(p, col_types = "text"))
    } else {
      data.table::fread(p, header = TRUE, colClasses = "character", na.strings = c("", "NA"))
    }
  }

  specs <- list()
  for (p in paths) {
    if (!file.exists(p)) stop(sprintf("%s info file not found: %s", role, p), call. = FALSE)
    dt <- read_table(p)
    keys <- norm(names(dt))
    get <- function(row, header) {
      idx <- which(keys == norm(header))
      if (length(idx) == 0) return(NA_character_) else as.character(row[[idx[1]]])
    }
    need <- c("Short name", "file.name", "snp.col", "ea.col", "nea.col", "beta.col", "se.col", "pval.col")
    missing_cols <- need[!(norm(need) %in% keys)]
    if (length(missing_cols))
      stop(sprintf("%s info file '%s' is missing required column(s): %s",
                   role, basename(p), paste(missing_cols, collapse = ", ")), call. = FALSE)

    for (i in seq_len(nrow(dt))) {
      row <- as.list(dt[i])
      raw_name <- trimws(get(row, "Short name"))
      if (blank(raw_name)) next
      name <- safe_name(raw_name)
      if (!identical(name, raw_name))
        warning(sprintf("Manifest '%s': Short name '%s' sanitised to '%s' for filenames.",
                        basename(p), raw_name, name), call. = FALSE)
      downloaded <- norm(get(row, "Downloaded"))
      gwas <- resolve(get(row, "file.name"))
      if (identical(downloaded, "no") || is.na(gwas)) {
        warning(sprintf("Manifest '%s': skipping trait '%s' (Downloaded=%s, file.name=%s).",
                        basename(p), raw_name, get(row, "Downloaded"), get(row, "file.name")), call. = FALSE)
        next
      }
      # Required per-row values: warn (don't skip) so selection stays predictable
      # and the grid's per-pair tryCatch turns any failure into a clear, non-fatal one.
      req_missing <- c(snp = "snp.col", ea = "ea.col", nea = "nea.col",
                       beta = "beta.col", se = "se.col", p = "pval.col")
      req_missing <- names(req_missing)[vapply(unname(req_missing), function(h) blank(get(row, h)), logical(1))]
      if (length(req_missing))
        warning(sprintf("Manifest '%s': trait '%s' has empty %s - its MR runs will fail until filled.",
                        basename(p), raw_name, paste(req_missing, collapse = "/")), call. = FALSE)

      nf <- n_field(get(row, "N")); cf <- n_field(get(row, "Ncases"))
      col_args <- list(
        snp  = get(row, "snp.col"),  beta = get(row, "beta.col"), se = get(row, "se.col"),
        ea   = get(row, "ea.col"),   nea  = get(row, "nea.col"),  p  = get(row, "pval.col"),
        chr  = if (!blank(get(row, "chr.col"))) get(row, "chr.col") else NULL,
        pos  = if (!blank(get(row, "pos.col"))) get(row, "pos.col") else NULL,
        eaf  = if (!blank(get(row, "eaf.col"))) get(row, "eaf.col") else NULL,
        n    = nf$col, ncase = cf$col
      )
      spec <- list(
        name         = name,
        gwas         = gwas,
        ivs          = resolve(get(row, "ivs.file")),
        col_args     = col_args,
        n_total      = nf$total,
        ncase_total  = cf$total,
        prevalence   = as_num(get(row, "Population prevalence")),
        type         = if (is.null(cf$total) && is.null(cf$col)) "continuous" else "binary",
        build        = get(row, "Genome build"),
        mhc_region   = mhc_for_build(get(row, "Genome build"))
      )
      if (is.na(spec$prevalence)) spec$prevalence <- NULL   # drop -> stays unset downstream
      specs[[name]] <- spec   # later manifest wins on duplicate Short name
    }
  }
  if (length(specs) == 0) stop(sprintf("No usable %s traits found in: %s", role, paste(paths, collapse = ", ")), call. = FALSE)

  if (!is.null(select)) {
    select <- unique(gsub("[^A-Za-z0-9._+-]+", "_", trimws(select)))  # match sanitised keys
    select <- select[nzchar(select)]
    absent <- setdiff(select, names(specs))
    if (length(absent))
      stop(sprintf("Requested %s trait(s) not found in the info file(s): %s\nAvailable: %s",
                   role, paste(absent, collapse = ", "), paste(names(specs), collapse = ", ")), call. = FALSE)
    specs <- specs[select]
  }
  specs
}

#' Read + select instruments for one exposure (clump/IV-list once, reusable
#' across many outcomes). Encapsulates the per-exposure block shared by
#' mr_pipeline.R and mr_grid.R, and writes the <prefix><name>_exposure_ivs.tsv
#' transparency table.
#'
#' @param exposure_file path to the exposure GWAS.
#' @param exposure_name trait label (used in filenames + results).
#' @param exp_col_args list(snp,beta,se,ea,nea,p,eaf,n,ncase,chr,pos).
#' @param opt run-level options (clump_*, ld_ref, plink_bin, f_stat, skip_clump,
#'   mhc_region, exclude_mhc, tmp_dir, cache_dir, out_prefix, exp_n_total,
#'   exp_ncase_total).
#' @param iv_list optional pre-defined instrument IDs (restricts, no clumping).
#' @return list(ivs_dat = formatted exposure instrument data.frame,
#'   raw = cleaned full exposure), or NULL if no instruments survived.
prepare_exposure_instruments <- function(exposure_file, exposure_name, exp_col_args, opt,
                                         iv_list = NULL) {
  exposure_raw <- read_gwas(
    gwas_file = exposure_file, type = "exposure", col_args = exp_col_args,
    trait_name = exposure_name, tmp_dir = opt$tmp_dir, pval_thresh = opt$clump_p,
    n_total = opt$exp_n_total, ncase_total = opt$exp_ncase_total, cache_dir = opt$cache_dir,
    iv_keep = iv_list)   # keep listed instruments even if sub-threshold in this GWAS
  message(sprintf("Successfully read exposure data for trait '%s'.", exposure_name))

  exposure_ivs_dat <- select_instruments(
    exposure_raw = exposure_raw, trait_name = exposure_name,
    clump_p = opt$clump_p, clump_kb = opt$clump_kb, clump_r2 = opt$clump_r2,
    ld_ref = opt$ld_ref, plink_bin = opt$plink_bin, min_f_stat = opt$f_stat,
    skip_clump = opt$skip_clump, iv_list = iv_list,
    mhc_region = opt$mhc_region, exclude_mhc = opt$exclude_mhc)
  if (!inherits(exposure_ivs_dat, "data.frame") || nrow(exposure_ivs_dat) == 0) return(NULL)
  message(sprintf("Final selected IVs count: %d", nrow(exposure_ivs_dat)))

  iv_cols <- intersect(
    c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure",
      "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
      "F_statistic", "mhc"),
    names(exposure_ivs_dat))
  data.table::fwrite(
    data.table::as.data.table(exposure_ivs_dat)[, ..iv_cols],
    paste0(opt$out_prefix, exposure_name, "_exposure_ivs.tsv"), sep = "\t", na = "NA")

  list(ivs_dat = exposure_ivs_dat, raw = exposure_raw)
}

#' Resolve the PLINK binary path.
#'
#' Uses an explicit path if given, else tries genetics.binaRies, else relies on
#' PLINK being on PATH (returns NULL so ieugwasr falls back to PATH).
#'
#' @param plink_bin Optional explicit path to the PLINK executable.
#' @return A path string, or NULL.
resolve_plink <- function(plink_bin = NULL) {
  if (!is.null(plink_bin)) return(plink_bin)
  suppressPackageStartupMessages({
    if (requireNamespace("genetics.binaRies", quietly = TRUE)) {
      p <- tryCatch(genetics.binaRies::get_plink_binary(), error = function(e) NULL)
      if (!is.null(p)) { message("Using PLINK binary found at: ", p); return(p) }
    }
  })
  warning("Could not find PLINK via genetics.binaRies. Relying on PLINK being on PATH or provide --plink_bin.", call. = FALSE)
  NULL
}

#' Select exposure instruments (p-value filter -> LD clumping -> format -> F-stat -> MHC).
#'
#' Single, coherent home for ALL instrument selection. Works on the RAW cleaned
#' exposure (from read_gwas), so LD clumping - which only needs rsID + p-value -
#' is done BEFORE TwoSampleMR::format_data. format_data therefore only ever runs
#' on the handful of selected instruments.
#'
#' Two modes:
#'   * default: LD-clump candidate SNPs to r2<clump_r2 (MR-standard independence).
#'   * skip_clump=TRUE: keep all candidates as-is (ONLY for inputs already
#'     independent at r2<0.001, e.g. pre-clumped cis-QTL instruments).
#'
#' @param exposure_raw Cleaned exposure data.table from read_gwas().
#' @param trait_name Exposure trait name.
#' @param clump_p P-value threshold for candidate instruments.
#' @param clump_kb,clump_r2 LD clumping window / r-squared threshold.
#' @param ld_ref PLINK bfile prefix for the LD reference panel.
#' @param plink_bin Optional PLINK binary path.
#' @param min_f_stat Minimum per-SNP F-statistic (beta^2/se^2).
#' @param skip_clump If TRUE, skip LD clumping (inputs assumed already independent).
#' @param mhc_region,exclude_mhc MHC handling (see flag_mhc_instruments()).
#'
#' @return A TwoSampleMR-formatted exposure data.frame of instruments.
select_instruments <- function(exposure_raw, trait_name,
                               clump_p = 5e-8, clump_kb = 10000, clump_r2 = 0.001,
                               ld_ref = NULL, plink_bin = NULL, min_f_stat = 10,
                               skip_clump = FALSE, iv_list = NULL,
                               mhc_region = "6:25000000-34000000", exclude_mhc = FALSE) {
  
  message("----- Selecting Instruments -----")
  exposure_raw <- data.table::as.data.table(exposure_raw)
  n_input <- nrow(exposure_raw)
  
  # 1. Candidate instruments: p < clump_p (drives clumping / skip-clump; a
  #    pre-defined --exp_ivs list does NOT require the p threshold, see below).
  candidates <- exposure_raw[!is.na(pval) & pval < clump_p]
  if (is.null(iv_list)) {
    if (nrow(candidates) == 0) {
      stop(sprintf("No SNPs found below the significance threshold p < %g.", clump_p), call. = FALSE)
    }
    message(sprintf("Found %d SNPs below p-value threshold %g.", nrow(candidates), clump_p))
    warn_non_rsid_instruments(candidates)
  }

  # 2. Independence: pre-defined list, LD clump (default), or skip.
  if (!is.null(iv_list)) {
    # A pre-defined instrument list is AUTHORITATIVE: use every listed SNP that is
    # PRESENT in the exposure sumstats (with beta/se), regardless of whether it
    # reaches clump_p in this particular GWAS (a published instrument is often
    # sub-threshold in a smaller or lifted-over file). No clumping (the list is
    # assumed independent at r2<0.001); the F-statistic filter below still guards
    # against weak instruments.
    ivs_raw <- exposure_raw[SNP %in% iv_list]
    n_listed <- length(unique(iv_list)); n_present <- length(unique(ivs_raw$SNP))
    if (n_present == 0) {
      stop("None of the pre-defined instruments (--exp_ivs) are present in the exposure sumstats (check the SNP-ID column and genome build).", call. = FALSE)
    }
    warn_non_rsid_instruments(ivs_raw)
    pv <- suppressWarnings(as.numeric(ivs_raw$pval)); n_sub <- sum(is.na(pv) | pv >= clump_p)
    message(sprintf("Using %d of %d pre-defined instruments present in the exposure (%d not found in the sumstats).",
                    n_present, n_listed, n_listed - n_present))
    if (n_sub > 0) {
      message(sprintf("  Note: %d of the %d present have p >= %g here; kept because --exp_ivs is authoritative (F>=%g still enforced).",
                      n_sub, n_present, clump_p, min_f_stat))
    }
    selected_snps <- unique(ivs_raw$SNP)
    warning("--exp_ivs: instruments taken as given (no clumping). Ensure the list is independent at MR standard (r2<0.001).", call. = FALSE)
  } else if (skip_clump) {
    warning("--skip_clump: treating inputs as already independent at MR standard (r2<0.001). ",
            "Gene-mapping / GCTA-COJO signal lists are typically only independent at r2<0.05 and are ",
            "NOT safe this way (correlated instruments understate IVW SEs). If unsure, re-run WITHOUT ",
            "--skip_clump to LD-clump at r2<0.001.", call. = FALSE)
    selected_snps <- candidates$SNP
    message(sprintf("Skipping LD clumping; keeping all %d candidate SNPs.", length(selected_snps)))
  } else {
    if (is.null(ld_ref)) stop("LD clumping requested but --ld_ref was not provided.", call. = FALSE)
    clump_input_df <- data.frame(
      rsid = candidates$SNP,
      pval = ifelse(candidates$pval == 0, .Machine$double.xmin, candidates$pval),
      stringsAsFactors = FALSE
    )
    plink_path <- resolve_plink(plink_bin)
    message(sprintf("Performing LD clumping with kb=%d, r2=%g using LD reference: %s",
                    clump_kb, clump_r2, ld_ref))
    clumped_snps_df <- tryCatch(
      ieugwasr::ld_clump(dat = clump_input_df, clump_kb = clump_kb, clump_r2 = clump_r2,
                         clump_p = 1, bfile = ld_ref, plink_bin = plink_path),
      error = function(e) stop(sprintf("Error during LD clumping: %s. Check PLINK installation, LD reference path (%s), and inputs.", e$message, ld_ref), call. = FALSE)
    )
    if (nrow(clumped_snps_df) == 0) {
      stop("Clumping removed all SNPs. Check parameters and LD reference.", call. = FALSE)
    }
    selected_snps <- clumped_snps_df$rsid
    message(sprintf("Identified %d independent instruments after clumping.", length(selected_snps)))
  }
  
  # 3. Format ONLY the selected instruments (small set -> no stack blow-up).
  #    For --exp_ivs, ivs_raw was already taken from the full present list above;
  #    for clump / skip-clump it comes from the p<clump_p candidates.
  if (is.null(iv_list)) ivs_raw <- candidates[SNP %in% selected_snps]
  exposure_ivs_dat <- format_gwas(ivs_raw, type = "exposure", trait_name = trait_name)
  
  # 4. F-statistic (F = beta^2 / se^2; needs no sample size) and filter
  exposure_ivs_dat <- exposure_ivs_dat %>%
    mutate(F_statistic = (beta.exposure^2) / (se.exposure^2))
  rows_before_f <- nrow(exposure_ivs_dat)
  exposure_ivs_dat <- exposure_ivs_dat %>% filter(F_statistic >= min_f_stat)
  if (nrow(exposure_ivs_dat) < rows_before_f) {
    message(sprintf("Removed %d IVs with F-statistic < %g.", rows_before_f - nrow(exposure_ivs_dat), min_f_stat))
  }
  if (nrow(exposure_ivs_dat) == 0) {
    stop(sprintf("No IVs remained after F-statistic filtering (F >= %g).", min_f_stat), call. = FALSE)
  }
  
  # 5. Flag / optionally drop MHC instruments (long-range LD + pleiotropy)
  exposure_ivs_dat <- flag_mhc_instruments(exposure_ivs_dat, mhc_region, exclude_mhc)

  if (!is.null(iv_list)) {
    message(sprintf(
      "Instrument attrition: %d listed -> %d present in sumstats -> %d after F>=%g%s.",
      length(unique(iv_list)), length(selected_snps), nrow(exposure_ivs_dat), min_f_stat,
      if (!exclude_mhc && "mhc" %in% names(exposure_ivs_dat)) sprintf(" (incl. %d MHC flagged)", sum(exposure_ivs_dat$mhc)) else ""))
  } else {
    sel_desc <- if (skip_clump) "kept (no clump)" else sprintf("clumped (r2<%g, %dkb)", clump_r2, clump_kb)
    message(sprintf(
      "Instrument attrition: %d input -> %d candidates (p<%g) -> %d %s -> %d after F>=%g%s.",
      n_input, nrow(candidates), clump_p, length(selected_snps), sel_desc,
      nrow(exposure_ivs_dat), min_f_stat,
      if (!exclude_mhc && "mhc" %in% names(exposure_ivs_dat)) sprintf(" (incl. %d MHC flagged)", sum(exposure_ivs_dat$mhc)) else ""))
  }
  message("----- Finished Selecting Instruments -----")
  
  return(exposure_ivs_dat)
}

#' Find LD proxies for a set of query SNPs using a local PLINK reference.
#'
#' Runs PLINK `--r2` with `--ld-snp-list` against the LD reference to list, for
#' each query SNP, all reference SNPs with r2 >= proxy_r2 within proxy_kb.
#'
#' @param query_snps Character vector of SNPs to find proxies for.
#' @param ld_ref PLINK bfile prefix.
#' @param plink_bin Optional PLINK binary path (else resolved automatically).
#' @param proxy_r2 Minimum r2.
#' @param proxy_kb Search window in kb.
#' @param tmp_dir Directory for temporary files.
#' @return data.table(query, proxy, r2), excluding self-pairs; empty if none.
find_proxies <- function(query_snps, ld_ref, plink_bin = NULL,
                         proxy_r2 = 0.8, proxy_kb = 1000, tmp_dir = tempdir(),
                         plink_mem = 30000) {
  query_snps <- unique(query_snps[!is.na(query_snps) & nzchar(query_snps)])
  if (length(query_snps) == 0) return(data.table::data.table())
  if (is.null(ld_ref)) { warning("Proxy search needs --ld_ref; skipping.", call. = FALSE); return(data.table::data.table()) }
  plink <- resolve_plink(plink_bin)
  if (is.null(plink)) { warning("No PLINK binary for proxy search; skipping proxies.", call. = FALSE); return(data.table::data.table()) }

  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  snp_file <- tempfile(tmpdir = tmp_dir, fileext = ".snps")
  out_pref <- tempfile(tmpdir = tmp_dir)
  writeLines(query_snps, snp_file)
  ok <- tryCatch(
    system2(plink, args = c("--bfile", ld_ref, "--r2",
                            "--ld-snp-list", snp_file,
                            "--ld-window-kb", format(proxy_kb, scientific = FALSE),
                            "--ld-window", "99999",
                            "--ld-window-r2", format(proxy_r2, scientific = FALSE),
                            "--memory", format(plink_mem, scientific = FALSE),
                            "--threads", "1",
                            "--out", out_pref),
            stdout = FALSE, stderr = FALSE),
    error = function(e) 1L)
  ld_file <- paste0(out_pref, ".ld")
  if (!file.exists(ld_file)) { warning("Proxy search produced no LD output (check --ld_ref).", call. = FALSE); return(data.table::data.table()) }
  ld <- tryCatch(data.table::fread(ld_file), error = function(e) NULL)
  if (is.null(ld) || nrow(ld) == 0) return(data.table::data.table())
  res <- ld[SNP_A != SNP_B, .(query = SNP_A, proxy = SNP_B, r2 = R2)]
  res[r2 >= proxy_r2][order(query, -r2)]
}

#' Substitute LD proxies for instruments missing from the outcome.
#'
#' For each exposure instrument absent from the outcome, picks the highest-r2
#' proxy that is present in BOTH the (full) exposure and the outcome, and swaps
#' the instrument for that proxy (using the proxy's own effects in both data
#' sets, so standard harmonisation applies).
#'
#' @param exposure_ivs_dat Formatted exposure instruments.
#' @param outcome_raw Cleaned outcome data.table (already subset to the IVs).
#' @param exposure_full Cleaned raw exposure (from read_gwas) to source proxy exposure effects.
#' @param outcome_file,out_col_args Outcome file + column map (to read proxy outcome records).
#' @param opt Pipeline options.
#' @return list(exposure, outcome_raw, map) with proxies substituted (map may be NULL).
apply_proxies <- function(exposure_ivs_dat, outcome_raw, exposure_full, outcome_file, out_col_args, opt) {
  unchanged <- list(exposure = exposure_ivs_dat, outcome_raw = outcome_raw, map = NULL)
  exposure_ivs_dat <- data.table::as.data.table(exposure_ivs_dat)
  missing <- setdiff(exposure_ivs_dat$SNP, outcome_raw$SNP)
  if (length(missing) == 0) { message("Proxies: all instruments present in outcome; none needed."); return(unchanged) }
  proxy_cache_key <- cache_token(source_fingerprint(outcome_file), sort(missing),
                                 exposure_ivs_dat$SNP, exposure_ivs_dat$beta.exposure, exposure_ivs_dat$se.exposure, opt$ld_ref,
                                 opt$proxy_r2, opt$proxy_kb, opt$f_stat, opt$clump_p)
  proxy_cache_file <- cache_path(opt$cache_dir, "proxy_substitutions", proxy_cache_key)
  cached_proxy <- cache_read(proxy_cache_file, "proxy substitution")
  if (!is.null(cached_proxy)) return(cached_proxy)
  message(sprintf("Proxies: %d instrument(s) missing from outcome; searching LD reference (r2>=%g, %gkb).",
                  length(missing), opt$proxy_r2, opt$proxy_kb))
  if (is.null(exposure_full)) { warning("Proxy search needs the full exposure; skipping.", call. = FALSE); return(unchanged) }
  exposure_full <- data.table::as.data.table(exposure_full)

  prox <- find_proxies(missing, opt$ld_ref, opt$plink_bin, opt$proxy_r2, opt$proxy_kb, opt$tmp_dir,
                       plink_mem = if (!is.null(opt$proxy_mem)) opt$proxy_mem else 30000)
  if (nrow(prox) == 0) { message("Proxies: no LD partner (r2>=", opt$proxy_r2, ") found for any missing instrument."); cache_write(unchanged, proxy_cache_file); return(unchanged) }
  message(sprintf("Proxies:   LD reference returned %d partner(s), covering %d/%d missing instruments.",
                  nrow(prox), length(unique(prox$query)), length(missing)))
  # Proxy must have an exposure effect and not already be an instrument.
  prox <- prox[proxy %in% exposure_full$SNP & !proxy %in% exposure_ivs_dat$SNP]
  if (nrow(prox) == 0) { message("Proxies: candidates found but none usable from the exposure."); cache_write(unchanged, proxy_cache_file); return(unchanged) }
  message(sprintf("Proxies:   %d partner(s) present in the exposure (covering %d/%d missing).",
                  nrow(prox), length(unique(prox$query)), length(missing)))
  # Proxy must be present in the outcome (read outcome records for candidates).
  # read_gwas stop()s if none of keep_snps are found, which is a valid outcome
  # here (no proxy in the outcome), so catch it and fall back to unchanged.
  proxy_out <- tryCatch(
    read_gwas(outcome_file, type = "outcome", col_args = out_col_args,
              trait_name = if ("trait" %in% names(outcome_raw)) outcome_raw$trait[1] else NULL,
              tmp_dir = opt$tmp_dir, keep_snps = unique(prox$proxy),
              n_total = opt$out_n_total, ncase_total = opt$out_ncase_total,
              cache_dir = opt$cache_dir),
    error = function(e) NULL)
  if (is.null(proxy_out) || nrow(proxy_out) == 0) {
    message("Proxies: candidate proxies were not found in the outcome; no substitution.")
    cache_write(unchanged, proxy_cache_file)
    return(unchanged)
  }
  prox <- prox[proxy %in% proxy_out$SNP]
  if (nrow(prox) == 0) { message("Proxies: candidates found but none present in the outcome."); cache_write(unchanged, proxy_cache_file); return(unchanged) }
  message(sprintf("Proxies:   %d partner(s) present in BOTH exposure and outcome (covering %d/%d missing).",
                  nrow(prox), length(unique(prox$query)), length(missing)))
  # A replacement SNP is an instrument in its own right: it must pass the
  # same minimum F-statistic as the original IVs.  `exposure_full` is already
  # p-value filtered at --clump_p, so this retains the required strong
  # exposure association (p < 5e-8 by default) as well as checking strength
  # explicitly.  Filter before choosing the highest-r2 partner so a weak top
  # partner does not prevent use of a lower-r2, valid proxy.
  proxy_exp_raw <- exposure_full[SNP %in% prox$proxy]
  proxy_exp_raw[, proxy_f_stat := (beta ^ 2) / (se ^ 2)]
  valid_proxy_ids <- proxy_exp_raw[is.finite(proxy_f_stat) & proxy_f_stat >= opt$f_stat,
                                   unique(SNP)]
  n_before_f <- nrow(prox)
  prox <- prox[proxy %in% valid_proxy_ids]
  if (nrow(prox) == 0) {
    message(sprintf("Proxies: none of the candidates passed the proxy F-statistic threshold (F >= %g); no substitution.",
                    opt$f_stat))
    cache_write(unchanged, proxy_cache_file)
    return(unchanged)
  }
  if (nrow(prox) < n_before_f) {
    message(sprintf("Proxies: removed %d candidate partner(s) with proxy F-statistic < %g.",
                    n_before_f - nrow(prox), opt$f_stat))
  }

  # One eligible proxy per missing instrument (highest r2).
  prox <- prox[order(query, -r2)][, .SD[1], by = query]

  # Format proxy exposure records and force the exposure id to match the set.
  proxy_exp_fmt <- data.table::as.data.table(
    format_gwas(proxy_exp_raw[SNP %in% prox$proxy], type = "exposure",
                trait_name = exposure_ivs_dat$exposure[1]))
  proxy_exp_fmt[, `:=`(exposure = exposure_ivs_dat$exposure[1],
                       id.exposure = exposure_ivs_dat$id.exposure[1])]

  exposure_new <- rbind(exposure_ivs_dat[!SNP %in% prox$query], proxy_exp_fmt, fill = TRUE)
  outcome_raw_new <- rbind(outcome_raw, proxy_out[SNP %in% prox$proxy], fill = TRUE)
  n_sub <- nrow(prox); n_missing <- length(missing)
  message(sprintf("Proxies: SUBSTITUTED %d of %d missing instrument(s); %d could not be rescued (no proxy in both datasets). Mapping -> <prefix>_proxies.tsv",
                  n_sub, n_missing, n_missing - n_sub))
  result <- list(exposure = exposure_new, outcome_raw = outcome_raw_new,
                 map = prox[, .(instrument = query, proxy, r2)])
  cache_write(result, proxy_cache_file)
  result
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
run_mr_analysis <- function(exposure_ivs_dat, outcome_file, outcome_name, out_col_args, out_prefix, opt,
                            exposure_full = NULL) {
  message(sprintf("Processing outcome: %s", outcome_file))
  runtime_start <- Sys.time()
  runtime <- data.table::data.table(stage = character(), seconds = numeric())
  mark_runtime <- function(stage) {
    runtime <<- rbind(runtime, data.table::data.table(stage = stage,
      seconds = as.numeric(difftime(Sys.time(), runtime_start, units = "secs"))))
  }
  # exp(beta) is an odds ratio only for a binary (log-odds) outcome. For a
  # continuous outcome we report the beta and its CI (lo_ci/up_ci) and omit the
  # meaningless OR columns. Controlled by --out_type (default "binary").
  outcome_is_binary <- !identical(tolower(trimws(as.character(opt$out_type))), "continuous")
  message(sprintf("Outcome '%s' treated as %s: %s.", outcome_name,
                  if (outcome_is_binary) "binary" else "continuous",
                  if (outcome_is_binary) "reporting odds ratios (exp(beta))" else "reporting beta, no OR"))
  outcome_raw <- read_gwas(
    gwas_file = outcome_file,
    type = "outcome",
    col_args = out_col_args,
    trait_name = outcome_name,
    tmp_dir = opt$tmp_dir,
    keep_snps = exposure_ivs_dat$SNP,
    n_total = opt$out_n_total,
    ncase_total = opt$out_ncase_total,
    cache_dir = opt$cache_dir
  )
  mark_runtime("outcome_read")
  # Optional: substitute LD proxies for instruments absent from the outcome.
  # Never let a proxy problem abort the MR: on any error, proceed without proxies.
  if (isTRUE(opt$proxies)) {
    pr <- tryCatch(
      apply_proxies(exposure_ivs_dat, outcome_raw, exposure_full, outcome_file, out_col_args, opt),
      error = function(e) { warning(sprintf("Proxy substitution failed (%s); proceeding without proxies.", e$message), call. = FALSE); NULL })
    if (!is.null(pr)) {
      exposure_ivs_dat <- pr$exposure
      outcome_raw <- pr$outcome_raw
      if (!is.null(pr$map) && nrow(pr$map) > 0) {
        data.table::fwrite(pr$map, paste0(out_prefix, "_proxies.tsv"), sep = "\t", na = "NA")
      }
    }
  }
  mark_runtime("proxy_substitution")
  outcome_dat <- format_gwas(outcome_raw, type = "outcome", trait_name = outcome_name)
  if (!inherits(outcome_dat, "data.frame") || nrow(outcome_dat) == 0) {
    warning(sprintf("No overlapping IVs found in outcome %s. Skipping.", outcome_file))
    return(NULL)
  }
  harmonized_cache_key <- cache_token(source_fingerprint(outcome_file), outcome_name,
                                      exposure_ivs_dat$SNP, exposure_ivs_dat$beta.exposure,
                                      exposure_ivs_dat$se.exposure, exposure_ivs_dat$effect_allele.exposure,
                                      exposure_ivs_dat$other_allele.exposure, opt$out_snp, opt$out_beta,
                                      opt$out_se, opt$out_ea, opt$out_nea)
  harmonized_cache_file <- cache_path(opt$cache_dir, "harmonized", harmonized_cache_key)
  harmonized_dat <- cache_read(harmonized_cache_file, "harmonized data")
  if (is.null(harmonized_dat)) {
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
    if (!is.null(harmonized_dat)) cache_write(harmonized_dat, harmonized_cache_file)
  }
  # TwoSampleMR versions differ in whether non-standard exposure columns
  # survive harmonisation.  Restore the MHC flag by SNP when necessary so the
  # downstream no-MHC branch remains deterministic (including on old caches).
  if (!is.null(harmonized_dat) && "mhc" %in% names(exposure_ivs_dat) &&
      !any(c("mhc", "mhc.exposure") %in% names(harmonized_dat)) && "SNP" %in% names(harmonized_dat)) {
    harmonized_dat$mhc.exposure <- exposure_ivs_dat$mhc[
      match(harmonized_dat$SNP, exposure_ivs_dat$SNP)]
  }
  mark_runtime("harmonisation")
  if (is.null(harmonized_dat) || nrow(harmonized_dat) == 0) {
    warning(sprintf("Harmonization failed or removed all SNPs for %s. Skipping.", outcome_file))
    return(NULL)
  }
  
  # Use the SAME instrument set for every method. harmonise_data(action=2)
  # keeps palindromic/ambiguous SNPs in the frame but marks mr_keep=FALSE; the
  # core TwoSampleMR methods honour that, but mr_presso / steiger_filtering do
  # not. Restrict to mr_keep==TRUE here so PRESSO and Steiger analyse exactly
  # the same SNPs as IVW/Egger/median.
  analysis_dat <- if ("mr_keep" %in% names(harmonized_dat)) {
    harmonized_dat[harmonized_dat$mr_keep == TRUE, , drop = FALSE]
  } else {
    harmonized_dat
  }
  if (isTRUE(opt$analysis_exclude_mhc)) {
    mhc_col <- intersect(c("mhc", "mhc.exposure"), names(analysis_dat))[1]
    if (is.na(mhc_col)) {
      warning("--analysis_exclude_mhc requested but no MHC flag was retained after harmonisation; analysing all SNPs.", call. = FALSE)
    } else {
      before_mhc <- nrow(analysis_dat)
      is_mhc <- as.logical(analysis_dat[[mhc_col]])
      analysis_dat <- analysis_dat[is.na(is_mhc) | !is_mhc, , drop = FALSE]
      message(sprintf("MHC sensitivity: removed %d harmonised instrument(s) using '%s'.",
                      before_mhc - nrow(analysis_dat), mhc_col))
    }
  }
  if (nrow(analysis_dat) == 0) {
    warning(sprintf("No usable SNPs after harmonisation (all mr_keep==FALSE) for %s. Skipping.", outcome_file))
    return(NULL)
  }
  
  n_snps <- nrow(analysis_dat)
  if (n_snps == 1) {
    mr_methods_list <- c("mr_wald_ratio")
  } else if (n_snps == 2) {
    mr_methods_list <- c("mr_ivw", "mr_weighted_median")
  } else if (n_snps >= 3) {
    mr_methods_list <- c("mr_ivw", "mr_weighted_median", "mr_egger_regression")
  }
  
  mr_results <- tryCatch({
    TwoSampleMR::mr(analysis_dat, method_list = mr_methods_list)
  }, error = function(e) {
    warning(sprintf("Error running core TwoSampleMR methods: %s", e$message), call.=FALSE)
    NULL
  })
  if (is.null(mr_results)) return(NULL)
  mark_runtime("core_mr")
  mr_results <- TwoSampleMR::generate_odds_ratios(mr_results)
  mr_results_dt <- data.table::as.data.table(mr_results)
  if (!outcome_is_binary) mr_results_dt[, c("or", "or_lci95", "or_uci95") := NULL]
  if (n_snps >= 3) {
    het_results <- tryCatch({
      data.table::as.data.table(TwoSampleMR::mr_heterogeneity(analysis_dat))
    }, error = function(e) NULL)
    if (!is.null(het_results)) {
      mr_results_dt <- merge(mr_results_dt, het_results[, .(id.exposure, id.outcome, method, Q, Q_df, Q_pval)],
                             by = c("id.exposure", "id.outcome", "method"), all.x = TRUE)
    }
    
    plt_results <- NULL
    if ("MR Egger" %in% mr_results_dt$method) {
      plt_results <- tryCatch({
        data.table::as.data.table(TwoSampleMR::mr_pleiotropy_test(analysis_dat))
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
  
  # Optional: for a binary trait (case counts present), tell Steiger to use the
  # log-odds r2 (get_r_from_lor). This is OPT-IN (--steiger_logodds) because it
  # is only self-consistent when BOTH traits are binary: mixing a log-odds
  # (liability-scale) r2 for one trait with a continuous r2 for the other is
  # not calibrated and can swing Steiger to remove ~all or ~none of the SNPs.
  # Affects only Steiger, not the MR estimates.
  if (isTRUE(opt$steiger_logodds)) {
    if ("ncase.exposure" %in% names(analysis_dat) &&
        any(is.finite(suppressWarnings(as.numeric(analysis_dat$ncase.exposure))))) {
      analysis_dat$units.exposure <- "log odds"
      if (!is.null(opt$exp_prevalence)) analysis_dat$prevalence.exposure <- opt$exp_prevalence
      message("--steiger_logodds: binary exposure -> Steiger uses the log-odds r2 formula.")
    }
    if ("ncase.outcome" %in% names(analysis_dat) &&
        any(is.finite(suppressWarnings(as.numeric(analysis_dat$ncase.outcome))))) {
      analysis_dat$units.outcome <- "log odds"
      if (!is.null(opt$out_prevalence)) analysis_dat$prevalence.outcome <- opt$out_prevalence
      message("--steiger_logodds: binary outcome -> Steiger uses the log-odds r2 formula.")
    }
  }

  # Steiger needs a valid (non-NA) sample size for BOTH traits; without it the
  # r2 values are undefined (NaN) and the test is meaningless, so skip loudly.
  steiger_has_n <- all(c("samplesize.exposure", "samplesize.outcome") %in% names(analysis_dat)) &&
    any(is.finite(suppressWarnings(as.numeric(analysis_dat$samplesize.exposure)))) &&
    any(is.finite(suppressWarnings(as.numeric(analysis_dat$samplesize.outcome))))
  if (opt$steiger && !steiger_has_n) {
    message("Steiger filtering skipped: needs a sample size for BOTH traits (map --exp_n/--out_n or set --exp_n_total/--out_n_total).")
  } else if (opt$steiger) {
    steiger_results <- tryCatch({
      TwoSampleMR::steiger_filtering(analysis_dat)
    }, error = function(e) { message(sprintf("Steiger filtering failed (needs sample sizes for both traits): %s", e$message)); NULL })
    if (is.null(steiger_results)) {
      message("Steiger filtering produced no result (skipped). Check that both traits have a sample size (--*_n / --*_n_total).")
    } else {
      steiger_filtered_data <- steiger_results %>%
        filter(steiger_dir == TRUE & steiger_pval < 0.05)
      n_kept <- nrow(steiger_filtered_data); n_removed <- nrow(analysis_dat) - n_kept
      message(sprintf("Steiger filtering: %d of %d SNPs pass the exposure->outcome direction test (removed %d).",
                      n_kept, nrow(analysis_dat), n_removed))
      if (n_removed > 0.5 * nrow(analysis_dat)) {
        message(sprintf("WARNING: Steiger removed a large fraction (%d/%d = %.0f%%). This often reflects mis-specified trait variance (e.g. a BINARY trait analysed without case counts + prevalence, so its r^2 is under-estimated) rather than true reverse causation. Interpret 'mr_ivw_steiger' with caution.",
                        n_removed, nrow(analysis_dat), 100 * n_removed / nrow(analysis_dat)))
      }
      if (n_kept > 0 && n_removed > 0) {
        steiger_ivw <- tryCatch({
          TwoSampleMR::mr(steiger_filtered_data, method_list = c("mr_ivw"))
        }, error = function(e) NULL)
        if (!is.null(steiger_ivw)) {
          steiger_ivw <- TwoSampleMR::generate_odds_ratios(steiger_ivw)
          steiger_ivw_dt <- data.table::as.data.table(steiger_ivw); steiger_ivw_dt$method <- "mr_ivw_steiger"
          if (!outcome_is_binary) steiger_ivw_dt[, c("or", "or_lci95", "or_uci95") := NULL]
          mr_results_dt <- rbind(mr_results_dt, steiger_ivw_dt, fill = TRUE)
          message("Added 'mr_ivw_steiger' (IVW on Steiger-passing SNPs).")
        }
      } else if (n_removed == 0) {
        message("Steiger removed no SNPs; the Steiger-filtered IVW equals the main IVW, so no separate 'mr_ivw_steiger' row is added.")
      } else {
        message("Steiger removed ALL SNPs; no Steiger IVW computed.")
      }
    }
  }
  
  # MR-PRESSO is costly, so for broad GWAS-by-GWAS screening only run it as a
  # follow-up to a nominally significant IVW result.
  ivw_pval <- mr_results_dt[method == "Inverse variance weighted", pval]
  ivw_significant <- length(ivw_pval) > 0 && any(is.finite(ivw_pval) & ivw_pval < 0.05)
  if (opt$presso && n_snps >= 4 && requireNamespace("MRPRESSO", quietly = TRUE) && ivw_significant) {
    nbdist <- if (!is.null(opt$presso_nbdist)) opt$presso_nbdist else 1000
    message(sprintf("Running MR-PRESSO on %d SNPs (NbDistribution=%d)... cost scales ~n_snps^2 x NbDistribution, so this can be slow for large/high settings (use --clump_r2 0.001 for fewer, independent instruments, or --no_presso to skip).",
                    n_snps, nbdist))
    presso_t0 <- Sys.time()
    presso_results <- tryCatch({
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome", SdExposure = "se.exposure",
        OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
        data = as.data.frame(analysis_dat),
        NbDistribution = nbdist,
        SignifThreshold = 0.05
      )
    }, error = function(e) NULL)
    message(sprintf("MR-PRESSO finished in %.1f min.", as.numeric(difftime(Sys.time(), presso_t0, units = "mins"))))
    if (!is.null(presso_results)) {
      global_pval <- presso_results$`MR-PRESSO results`$`Global Test`$Pvalue
      message(sprintf("MR-PRESSO Global P-value: %s", global_pval))
      presso_main <- data.table::as.data.table(presso_results$`Main MR results`)
      presso_main[, b := `Causal Estimate`]; presso_main[, se := Sd]; presso_main[, pval := `P-value`]
      presso_main[, method := ifelse(`MR Analysis` == "Raw", "mr_presso_raw", "mr_presso_corrected")]
      presso_main[, id.exposure := analysis_dat$id.exposure[1]]; presso_main[, id.outcome := analysis_dat$id.outcome[1]]
      presso_main[, exposure := analysis_dat$exposure[1]]; presso_main[, outcome := analysis_dat$outcome[1]]
      outlier_indices <- presso_results$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
      n_outliers <- if (is.null(outlier_indices) || any(is.na(outlier_indices))) 0 else length(outlier_indices)
      presso_main[method == "mr_presso_raw", nsnp := nrow(analysis_dat)]
      presso_main[method == "mr_presso_corrected", nsnp := nrow(analysis_dat) - n_outliers]
      presso_main[, lo_ci := b - 1.96 * se]; presso_main[, up_ci := b + 1.96 * se]
      cols_to_keep <- c("id.exposure", "id.outcome", "exposure", "outcome", "method", "nsnp", "b", "se", "pval", "lo_ci", "up_ci")
      if (outcome_is_binary) {
        presso_main[, or := exp(b)]; presso_main[, or_lci95 := exp(b - 1.96 * se)]; presso_main[, or_uci95 := exp(b + 1.96 * se)]
        cols_to_keep <- c(cols_to_keep, "or", "or_lci95", "or_uci95")
      }
      if (!is.null(presso_results$`MR-PRESSO results`$`Distortion Test`$Pvalue)) { presso_main[, distortion_pval := presso_results$`MR-PRESSO results`$`Distortion Test`$Pvalue]; cols_to_keep <- c(cols_to_keep, "distortion_pval") }
      if (!is.null(presso_results$`MR-PRESSO results`$`Global Test`$Pvalue)) { presso_main[, presso_global_pval := presso_results$`MR-PRESSO results`$`Global Test`$Pvalue]; cols_to_keep <- c(cols_to_keep, "presso_global_pval") }
      mr_results_dt <- rbind(mr_results_dt, presso_main[, ..cols_to_keep], fill = TRUE)
    }
  } else if (opt$presso && n_snps >= 4 && requireNamespace("MRPRESSO", quietly = TRUE)) {
    ivw_label <- if (length(ivw_pval) == 0) "not available" else format(ivw_pval[1], digits = 4)
    message(sprintf("MR-PRESSO skipped: IVW p-value is %s (requires p < 0.05).", ivw_label))
  }
  mark_runtime("sensitivity_and_presso")
  
  harmonized_file <- paste0(out_prefix, "_harmonized_data.rds")
  saveRDS(harmonized_dat, file = harmonized_file)
  
  results_file <- paste0(out_prefix, "_full_mr_results.csv")
  data.table::fwrite(mr_results_dt, results_file, sep = ",", na = "NA")
  runtime[, total_seconds := seconds]
  data.table::fwrite(runtime, paste0(out_prefix, "_runtime.tsv"), sep = "\t")
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
  # Drop columns we don't carry forward (only if present). lo_ci/up_ci (the beta
  # CI) are kept, so continuous-outcome results retain a CI even without OR cols.
  drop_cols <- intersect(c("id.exposure", "id.outcome",
                           "egger_intercept", "egger_intercept_se"), names(all_mr_results))
  if (length(drop_cols)) all_mr_results[, (drop_cols) := NULL]
  
  # Several columns are only present depending on which methods/sensitivity
  # tests produced output for a given pair (e.g. MR-PRESSO with no outliers
  # omits distortion_pval; Wald-only pairs have no Egger/heterogeneity columns).
  # Add any that are missing as NA so the references below are always valid.
  ensure_cols <- c("Q", "Q_df", "Q_pval", "egger_intercept_pval",
                   "distortion_pval", "presso_global_pval")
  for (col in ensure_cols) if (!col %in% names(all_mr_results)) all_mr_results[, (col) := NA_real_]
  
  # Process each exposure-outcome pair INDEPENDENTLY. (The previous version
  # mutated the whole IVW table with a single pair's values inside the loop,
  # which recycled wrong flags across pairs and crashed when the number of
  # MR-PRESSO rows - PRESSO only runs when IVW p<0.05 - differed from the number
  # of IVW rows.)
  get_beta <- function(dt, m) { v <- dt[method == m, b]; if (length(v)) v[1] else NA_real_ }
  first_or_na <- function(v) if (length(v)) v[1] else NA_real_

  pairs <- unique(all_mr_results[, .(exposure, outcome)])
  out_rows <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    exp <- pairs$exposure[i]; out <- pairs$outcome[i]
    tmp.res <- all_mr_results[exposure == exp & outcome == out]
    if (nrow(tmp.res) == 0) next
    has_presso  <- any(grepl("mr_presso", tmp.res$method))
    has_steiger <- "mr_ivw_steiger" %in% tmp.res$method

    ivw  <- tmp.res[method == "Inverse variance weighted" & nsnp > 2]
    wald <- tmp.res[method == "Wald ratio"]

    # Pleiotropy: MR-Egger intercept p < 0.05
    egger_p <- first_or_na(tmp.res[method == "MR Egger", egger_intercept_pval])
    flag_pleiotropy <- isTRUE(is.finite(egger_p) && egger_p < 0.05)

    # Heterogeneity: I2 = (Q - Q_df)/Q > 0.5 on the IVW row
    if (nrow(ivw) > 0) {
      ivw[, I2 := (Q - Q_df) / Q]
      ivw[, egger_intercept_pval := egger_p]
      if (has_presso) {
        ivw[, presso_global_pval := first_or_na(tmp.res[method == "mr_presso_raw", presso_global_pval])]
        ivw[, presso_distortion_pval := first_or_na(tmp.res[method == "mr_presso_raw", distortion_pval])]
      }
    }
    flag_het <- nrow(ivw) > 0 && isTRUE(any(is.finite(ivw$I2) & ivw$I2 > 0.5))

    # Steiger directionality discordance (IVW vs Steiger-filtered IVW disagree in sign)
    steiger_diff <- NA
    if (has_steiger) {
      sgn <- sign(tmp.res[method %in% c("mr_ivw_steiger", "Inverse variance weighted"), b])
      sgn <- sgn[is.finite(sgn)]
      steiger_diff <- length(unique(sgn)) > 1
    }

    # Sensitivity betas -> direction concordance vs IVW
    presso_beta <- NA_real_
    if (has_presso) {
      if ("mr_presso_corrected" %in% tmp.res$method) {
        dval <- tmp.res[method == "mr_presso_corrected", distortion_pval]
        dval <- suppressWarnings(first_or_na(ifelse(dval == "<0.001", 1e-4, as.numeric(dval))))
        presso_beta <- if (!is.na(dval) && dval < 0.05) get_beta(tmp.res, "mr_presso_corrected")
                       else get_beta(tmp.res, "mr_presso_raw")
      }
    }
    beta.sensitivity <- c(
      if (has_presso) c(MRPRESSO = presso_beta) else NULL,
      MREgger = get_beta(tmp.res, "MR Egger"),
      WeightedMedian = get_beta(tmp.res, "Weighted median"),
      if (has_steiger) c(Steiger = get_beta(tmp.res, "mr_ivw_steiger")) else NULL
    )
    ivw_beta <- if (nrow(ivw) > 0) first_or_na(ivw$b) else NA_real_
    diff_dir <- sum(sign(unlist(beta.sensitivity)) != sign(ivw_beta), na.rm = TRUE) != 0

    pair.out <- rbind(ivw, wald, fill = TRUE)
    if (nrow(pair.out) == 0) next
    pair.out[, DiffDirection := diff_dir]
    pair.out[, FlagPleiotropy := flag_pleiotropy]
    pair.out[, FlagHeterogeneity := flag_het]
    pair.out[, FlagSteiger := steiger_diff]
    out_rows[[i]] <- pair.out
  }

  all_mr_results.all <- data.table::rbindlist(Filter(Negate(is.null), out_rows), fill = TRUE)
  if (nrow(all_mr_results.all) > 0) {
    if ("distortion_pval" %in% names(all_mr_results.all)) all_mr_results.all[, distortion_pval := NULL]
    all_mr_results.all[, p_value_fdr := p.adjust(pval, method = "fdr")]
  }
  results_file <- paste0(opt$out_prefix, "all_processed_mr_results.csv")
  data.table::fwrite(all_mr_results.all, results_file, sep = ",", na = "NA")

  return(all_mr_results.all)
}
