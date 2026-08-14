#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# mr_grid.R - manifest-driven bidirectional MR grid.
#
# Instead of hand-coding each trait's column mappings / sample sizes / build on
# the command line (as run_mr_slurm_eos.sh does with EO_EXP, ASTHMA_EXP, ...),
# you point this driver at one or more "info files" (manifests) describing the
# exposure traits and the outcome traits. It then runs the full
# exposure x outcome grid - both directions (--bidirectional) and, if asked, the
# 4-config sensitivity grid (--sensitivity_grid: lenient/strict x +/-MHC) - in a
# single command, reading every column/N/type/build detail from the manifest.
#
# See read_trait_manifest() in pipeline_functions.R for the manifest schema.
#
# Example:
#   Rscript mr_grid.R \
#     --exp_info autoimmune.csv --out_info menopause.csv \
#     --data_dir /path/to/gwas --ld_ref /path/to/ld/Phase3_eur \
#     --bidirectional --sensitivity_grid \
#     --out_prefix output/autoimmune/
# ---------------------------------------------------------------------------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  getwd()
}
SCRIPT_DIR <- get_script_dir()

# Resolve the R library BEFORE loading packages (same policy as mr_pipeline.R).
.raw_args <- commandArgs(trailingOnly = TRUE)
.libpath_idx <- which(.raw_args == "--lib_path")
.repo_lib_candidates <- c(file.path(dirname(SCRIPT_DIR), "Rpackages"),
                          file.path(SCRIPT_DIR, "Rpackages"))
.repo_lib <- .repo_lib_candidates[dir.exists(.repo_lib_candidates)][1]
if (length(.libpath_idx) == 1 && length(.raw_args) >= .libpath_idx + 1) {
  .user_lib <- .raw_args[.libpath_idx + 1]
  if (dir.exists(.user_lib)) { .libPaths(.user_lib); message(sprintf("Using R library path (--lib_path): %s", .user_lib)) }
  else warning(sprintf("--lib_path '%s' does not exist; ignoring.", .user_lib), call. = FALSE)
} else if (!is.na(.repo_lib)) {
  .libPaths(.repo_lib); message(sprintf("Using repo-local R library: %s", .repo_lib))
}

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
})

option_list <- list(
  # ---- Manifests + trait selection ----
  make_option("--exp_info", type = "character",
              help = "Comma-separated info file(s) describing EXPOSURE traits (required). See README 'Trait info files'."),
  make_option("--out_info", type = "character",
              help = "Comma-separated info file(s) describing OUTCOME traits (required)."),
  make_option("--data_dir", type = "character", default = NULL,
              help = "Base directory used to resolve relative file.name / ivs.file entries in the manifests."),
  make_option("--exposures", type = "character", default = NULL,
              help = "Comma-separated Short names to restrict the exposure set (default: all traits in --exp_info)."),
  make_option("--outcomes", type = "character", default = NULL,
              help = "Comma-separated Short names to restrict the outcome set (default: all traits in --out_info)."),
  make_option("--bidirectional", action = "store_true", default = FALSE,
              help = "Also run outcome-traits-as-exposures (reverse direction) for every pair."),
  make_option("--sensitivity_grid", action = "store_true", default = FALSE,
              help = "Run all 4 configs per pair: lenient (ivs.file) / strict (r2=0.001) x with / without MHC. Output prefixes: '', 'rigid_', 'noMHC_', 'noMHC_rigid_'."),
  make_option("--dry_run", action = "store_true", default = FALSE,
              help = "Print the resolved traits and the grid that WOULD run, then exit without running MR."),
  # ---- Output ----
  make_option("--out_prefix", type = "character",
              help = "Output path prefix (e.g. 'output/autoimmune/'); the config prefix and '<exp>_vs_<out>' are appended (required)."),
  # ---- Run-level knobs (shared with mr_pipeline.R) ----
  make_option("--ld_ref", type = "character", help = "Path prefix to LD reference panel (PLINK bfile) for clumping/proxies."),
  make_option("--plink_bin", type = "character", default = NULL, help = "Path to PLINK binary."),
  make_option("--clump_p", type = "numeric", default = 5e-8, help = "P-value threshold for IV selection [default: %default]."),
  make_option("--clump_kb", type = "numeric", default = 10000, help = "Kilobase radius for LD clumping [default: %default]."),
  make_option("--clump_r2", type = "numeric", default = 0.001,
              help = "R2 for LD clumping in a SINGLE (non-grid) run [default: %default]. In --sensitivity_grid the lenient/strict r2 are fixed at 0.2/0.001."),
  make_option("--exclude_mhc", action = "store_true", default = FALSE,
              help = "For a SINGLE (non-grid) run, drop MHC instruments. --sensitivity_grid controls MHC per config."),
  make_option("--mhc_region", type = "character", default = NULL,
              help = "Override the MHC region CHR:START-END. Default: chosen per trait from its Genome build (b37/b38)."),
  make_option("--proxies", action = "store_true", default = FALSE,
              help = "Substitute LD proxies for instruments absent from the outcome (needs --ld_ref)."),
  make_option("--proxy_r2", type = "numeric", default = 0.8, help = "Minimum r2 for an LD proxy [default: %default]."),
  make_option("--proxy_kb", type = "numeric", default = 1000, help = "Window (kb) to search for LD proxies [default: %default]."),
  make_option("--proxy_mem", type = "numeric", default = 30000, help = "Memory (MB) cap for the PLINK proxy search [default: %default]."),
  make_option("--f_stat", type = "numeric", default = 10, help = "Minimum F-statistic for IVs [default: %default]."),
  make_option("--presso_nbdist", type = "numeric", default = 1000, help = "MR-PRESSO simulations (NbDistribution) [default: %default]."),
  make_option("--tmp_dir", type = "character", default = "./tmp_pipeline", help = "Temporary directory [default: %default]."),
  make_option("--cache_dir", type = "character", default = "./mr_cache", help = "Reusable cache directory [default: %default]."),
  make_option("--no_cache", action = "store_true", default = FALSE, help = "Disable reusable on-disk caches."),
  make_option("--cache_max_gb", type = "numeric", default = 20, help = "Cap the on-disk cache at this many GB (least-recently-used files evicted at start-up); 0 = unlimited [default: %default]."),
  make_option("--lib_path", type = "character", default = NULL, help = "Optional custom R library (applied before packages load)."),
  make_option("--steiger", type = "logical", action = "store_true", default = TRUE, help = "Run Steiger filtering."),
  make_option("--no_steiger", action = "store_false", dest = "steiger"),
  make_option("--steiger_logodds", action = "store_true", default = FALSE,
              help = "Use log-odds r2 for binary traits in Steiger (only self-consistent when BOTH traits are binary)."),
  make_option("--presso", type = "logical", action = "store_true", default = TRUE, help = "Run MR-PRESSO when IVW p<0.05."),
  make_option("--no_presso", action = "store_false", dest = "presso")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (isTRUE(opt$no_cache)) opt$cache_dir <- NULL

# --skip_clump is never appropriate here (the grid decides selection per config).
opt$skip_clump <- FALSE
# --exp_ivs on the shared helpers comes from the manifest, not a global flag.
opt$exp_ivs <- NULL; opt$exp_ivs_col <- NULL

# Source helper functions relative to this script (works from any CWD).
functions_candidates <- c(
  file.path(SCRIPT_DIR, "pipeline_functions.R"),
  "pipeline_functions.R",
  file.path("pipeline_scripts", "pipeline_functions.R"))
functions_file <- functions_candidates[file.exists(functions_candidates)][1]
if (is.na(functions_file))
  stop(sprintf("Helper functions file 'pipeline_functions.R' not found in any of: %s",
               paste(functions_candidates, collapse = ", ")), call. = FALSE)
source(functions_file)

# Bound the reusable cache before doing any work (never evicts this run's data).
prune_cache(opt$cache_dir, opt$cache_max_gb)

if (isTRUE(opt$presso)) {
  if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
    warning("MR-PRESSO requested (--presso) but package 'MRPRESSO' is not installed. Skipping MR-PRESSO.", call. = FALSE)
    opt$presso <- FALSE
  } else suppressPackageStartupMessages(library(MRPRESSO))
}

# ---- Validate required args -------------------------------------------------
for (req in c("exp_info", "out_info", "out_prefix")) {
  if (is.null(opt[[req]])) stop(sprintf("--%s is required.", req), call. = FALSE)
}
split_csv <- function(x) if (is.null(x)) NULL else trimws(strsplit(x, ",")[[1]])

# ---- Read manifests ---------------------------------------------------------
exp_specs <- read_trait_manifest(split_csv(opt$exp_info), data_dir = opt$data_dir,
                                 select = split_csv(opt$exposures), role = "exposure")
out_specs <- read_trait_manifest(split_csv(opt$out_info), data_dir = opt$data_dir,
                                 select = split_csv(opt$outcomes), role = "outcome")

# ---- Config specs -----------------------------------------------------------
# use_ivs=TRUE  -> prefer the exposure's ivs.file (fallback: clump at clump_r2).
# use_ivs=FALSE -> force re-clump at clump_r2 (the "strict"/robustness configs).
if (isTRUE(opt$sensitivity_grid)) {
  configs <- list(
    list(prefix = "",             use_ivs = TRUE,  clump_r2 = 0.2,   exclude_mhc = FALSE),
    list(prefix = "rigid_",       use_ivs = FALSE, clump_r2 = 0.001, exclude_mhc = FALSE),
    list(prefix = "noMHC_",       use_ivs = TRUE,  clump_r2 = 0.2,   exclude_mhc = TRUE),
    list(prefix = "noMHC_rigid_", use_ivs = FALSE, clump_r2 = 0.001, exclude_mhc = TRUE))
} else {
  configs <- list(
    list(prefix = "", use_ivs = TRUE, clump_r2 = opt$clump_r2, exclude_mhc = isTRUE(opt$exclude_mhc)))
}

# Directions: forward always; reverse only with --bidirectional.
directions <- list(list(tag = "forward", exp = exp_specs, out = out_specs))
if (isTRUE(opt$bidirectional))
  directions[[length(directions) + 1]] <- list(tag = "reverse", exp = out_specs, out = exp_specs)

# ---- Show / dry-run the resolved plan --------------------------------------
spec_line <- function(s) sprintf("  %-16s type=%-10s build=%-4s ivs=%s\n    file=%s",
  s$name, s$type, ifelse(is.null(s$build) || is.na(s$build), "?", s$build),
  ifelse(is.na(s$ivs), "(none)", basename(s$ivs)), s$gwas)
cat("\n=== Exposure traits ===\n"); for (s in exp_specs) cat(spec_line(s), "\n")
cat("\n=== Outcome traits ===\n"); for (s in out_specs) cat(spec_line(s), "\n")
cat(sprintf("\n=== Grid ===\nconfigs: %s\ndirections: %s\n",
            paste(vapply(configs, function(c) ifelse(nzchar(c$prefix), c$prefix, "lenient"), ""), collapse = ", "),
            paste(vapply(directions, `[[`, "", "tag"), collapse = ", ")))

planned <- 0L
for (cfg in configs) for (d in directions) for (e in d$exp) for (o in d$out) if (e$name != o$name) planned <- planned + 1L
cat(sprintf("planned MR runs (skipping self-pairs): %d\n\n", planned))

if (isTRUE(opt$dry_run)) { message("--dry_run set; not running MR."); quit(status = 0) }

# Make sure the output directory exists.
out_par <- dirname(opt$out_prefix)
if (nzchar(out_par) && !dir.exists(out_par)) dir.create(out_par, recursive = TRUE, showWarnings = FALSE)

# ---- Run the grid -----------------------------------------------------------
for (cfg in configs) {
  cfg_prefix <- paste0(opt$out_prefix, cfg$prefix)
  results_cfg <- list()

  for (d in directions) {
    for (e in d$exp) {
      # Per-config, per-exposure option clone (instruments are exposure-side).
      opt_e <- opt
      opt_e$out_prefix          <- cfg_prefix
      opt_e$exclude_mhc         <- FALSE                 # flagging done via mhc_region; drop happens at analysis
      opt_e$analysis_exclude_mhc <- cfg$exclude_mhc
      opt_e$clump_r2            <- cfg$clump_r2
      opt_e$mhc_region          <- if (!is.null(opt$mhc_region)) opt$mhc_region else e$mhc_region
      opt_e$exp_n_total         <- e$n_total
      opt_e$exp_ncase_total     <- e$ncase_total
      opt_e$exp_prevalence      <- e$prevalence

      # Lenient configs use the trait's pre-defined instrument list; strict re-clumps.
      iv_list <- NULL
      if (isTRUE(cfg$use_ivs) && !is.na(e$ivs)) {
        if (file.exists(e$ivs)) {
          iv_list <- read_iv_list(e$ivs)
          message(sprintf("[%s%s] using %d pre-defined instruments from %s.",
                          cfg$prefix, e$name, length(iv_list), basename(e$ivs)))
        } else {
          warning(sprintf("[%s%s] ivs.file not found (%s); clumping at r2=%s instead.",
                          cfg$prefix, e$name, e$ivs, cfg$clump_r2), call. = FALSE)
        }
      }
      if (isTRUE(cfg$use_ivs) && is.null(iv_list) && is.na(e$ivs))
        message(sprintf("[%s%s] no ivs.file; lenient clump at r2=%s.", cfg$prefix, e$name, cfg$clump_r2))

      prep <- tryCatch(
        prepare_exposure_instruments(e$gwas, e$name, e$col_args, opt_e, iv_list = iv_list),
        error = function(err) { warning(sprintf("[%s%s] exposure prep failed: %s", cfg$prefix, e$name, err$message), call. = FALSE); NULL })
      if (is.null(prep)) { warning(sprintf("[%s%s] no instruments; skipping this exposure.", cfg$prefix, e$name), call. = FALSE); next }

      for (o in d$out) {
        if (e$name == o$name) next                       # no self-MR
        this_prefix <- paste0(cfg_prefix, e$name, "_vs_", o$name)

        opt_eo <- opt_e
        opt_eo$out_type        <- o$type
        opt_eo$out_n_total     <- o$n_total
        opt_eo$out_ncase_total <- o$ncase_total
        opt_eo$out_prevalence  <- o$prevalence
        # Column names used inside run_mr_analysis for the harmonise cache key.
        opt_eo$out_snp  <- o$col_args$snp;  opt_eo$out_beta <- o$col_args$beta
        opt_eo$out_se   <- o$col_args$se;   opt_eo$out_ea   <- o$col_args$ea
        opt_eo$out_nea  <- o$col_args$nea
        # Per-SNP n/ncase come from the manifest when N/Ncases named a column
        # (e.g. aam/neb); otherwise they are NULL and the *_total constants apply.
        out_col_args <- list(
          snp = o$col_args$snp, beta = o$col_args$beta, se = o$col_args$se,
          ea = o$col_args$ea, nea = o$col_args$nea, p = o$col_args$p,
          eaf = o$col_args$eaf, n = o$col_args$n, ncase = o$col_args$ncase,
          chr = o$col_args$chr, pos = o$col_args$pos)

        message(sprintf("\n=== [%s] %s%s -> %s ===", format(Sys.time(), "%T"), cfg$prefix, e$name, o$name))
        res <- tryCatch(
          run_mr_analysis(prep$ivs_dat, o$gwas, o$name, out_col_args, this_prefix, opt_eo,
                          exposure_full = prep$raw),
          error = function(err) { warning(sprintf("[%s%s -> %s] failed: %s", cfg$prefix, e$name, o$name, err$message), call. = FALSE); NULL })
        if (!is.null(res)) results_cfg[[length(results_cfg) + 1]] <- res
      }
    }
  }

  # One processed-summary per config (mirrors a single mr_pipeline.R invocation).
  if (length(results_cfg) > 0) {
    opt_cfg <- opt; opt_cfg$out_prefix <- cfg_prefix
    process_mr_results(data.table::rbindlist(results_cfg, fill = TRUE), opt_cfg)
  } else {
    warning(sprintf("Config '%s' produced no results.", ifelse(nzchar(cfg$prefix), cfg$prefix, "lenient")), call. = FALSE)
  }
}

message("All grid analyses complete.")
