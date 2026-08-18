#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# mr_mvmr.R - manifest-driven multivariable MR (MVMR).
#
# One invocation runs ONE MVMR model: all traits named in --exposures go in
# together, estimated jointly against each trait in --outcomes. Columns, N,
# type, build and instrument lists come from the same info files (manifests)
# used by mr_grid.R. Also emits the univariate IVW for each exposure->outcome
# and a univariable-vs-multivariable forest.
#
# Example:
#   Rscript mr_mvmr.R --exp_info respiratory.xlsx --out_info reproductive.xlsx \
#     --exposures asthma,smoking --outcomes poi,anm \
#     --ld_ref /path/Phase3_eur --out_prefix output/mvmr/
# ---------------------------------------------------------------------------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  getwd()
}
SCRIPT_DIR <- get_script_dir()

.raw_args <- commandArgs(trailingOnly = TRUE)
.libpath_idx <- which(.raw_args == "--lib_path")
.repo_lib_candidates <- c(file.path(dirname(SCRIPT_DIR), "Rpackages"), file.path(SCRIPT_DIR, "Rpackages"))
.repo_lib <- .repo_lib_candidates[dir.exists(.repo_lib_candidates)][1]
if (length(.libpath_idx) == 1 && length(.raw_args) >= .libpath_idx + 1) {
  .user_lib <- .raw_args[.libpath_idx + 1]
  if (dir.exists(.user_lib)) { .libPaths(.user_lib); message(sprintf("Using R library path (--lib_path): %s", .user_lib)) }
} else if (!is.na(.repo_lib)) { .libPaths(.repo_lib); message(sprintf("Using repo-local R library: %s", .repo_lib)) }

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(data.table); library(TwoSampleMR); library(ieugwasr)
})

option_list <- list(
  make_option("--exp_info", type = "character", help = "Comma-separated info file(s) describing the exposure traits (required)."),
  make_option("--out_info", type = "character", help = "Comma-separated info file(s) describing the outcome traits (required)."),
  make_option("--data_dir", type = "character", default = NULL, help = "Base dir for relative file.name / ivs.file entries."),
  make_option("--exposures", type = "character", help = "Comma-separated Short names of the exposures to model JOINTLY (>=2; required)."),
  make_option("--outcomes", type = "character", help = "Comma-separated Short names of the outcome(s) (required)."),
  make_option("--out_prefix", type = "character", help = "Output path prefix, e.g. 'output/mvmr/' (required)."),
  make_option("--ld_ref", type = "character", help = "PLINK bfile prefix for clumping (required)."),
  make_option("--plink_bin", type = "character", default = NULL, help = "Path to PLINK binary."),
  make_option("--clump_p", type = "numeric", default = 5e-8, help = "P threshold for per-exposure instrument selection [default: %default]."),
  make_option("--clump_r2", type = "numeric", default = 0.01, help = "R2 for clumping / union re-clump [default: %default]."),
  make_option("--clump_kb", type = "numeric", default = 1000, help = "Clumping window in kb [default: %default]."),
  make_option("--mv_pval_threshold", type = "numeric", default = 1,
              help = "pval_threshold passed to TwoSampleMR::mv_multiple. Default 1 uses the full curated instrument set (selection already done upstream) [default: %default]."),
  make_option("--exclude_mhc", action = "store_true", default = FALSE, help = "Drop MHC instruments (flagged via each trait's build-specific MHC region)."),
  make_option("--mhc_region", type = "character", default = NULL, help = "Override MHC region CHR:START-END (default: per-trait from build)."),
  make_option("--f_stat", type = "numeric", default = 10, help = "Minimum per-SNP F for the clumped exposures [default: %default]."),
  make_option("--no_univariate", action = "store_true", default = FALSE, help = "Skip the univariate IVW overlay."),
  make_option("--egger", action = "store_true", default = FALSE, help = "Also run MVMR-Egger (needs MendelianRandomization)."),
  make_option("--no_plot", action = "store_true", default = FALSE, help = "Skip the univariable-vs-multivariable forest plot."),
  make_option("--dry_run", action = "store_true", default = FALSE, help = "Print the resolved traits + planned runs and exit."),
  make_option("--tmp_dir", type = "character", default = "./tmp_pipeline", help = "Temp dir [default: %default]."),
  make_option("--cache_dir", type = "character", default = "./mr_cache", help = "Cache dir [default: %default]."),
  make_option("--no_cache", action = "store_true", default = FALSE, help = "Disable on-disk caches."),
  make_option("--cache_max_gb", type = "numeric", default = 20, help = "Cache size cap in GB (0=unlimited) [default: %default]."),
  make_option("--lib_path", type = "character", default = NULL, help = "Optional custom R library.")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (isTRUE(opt$no_cache)) opt$cache_dir <- NULL
options(warn = 1)
opt$egger <- isTRUE(opt$egger)

# Source helpers (data layer + MVMR engine).
src <- function(f) {
  cand <- c(file.path(SCRIPT_DIR, f), f, file.path("pipeline_scripts", f))
  hit <- cand[file.exists(cand)][1]
  if (is.na(hit)) stop(sprintf("Cannot find %s in: %s", f, paste(cand, collapse = ", ")), call. = FALSE)
  source(hit)
}
src("pipeline_functions.R"); src("mvmr_functions.R")
prune_cache(opt$cache_dir, opt$cache_max_gb)

for (req in c("exp_info", "out_info", "exposures", "outcomes", "out_prefix", "ld_ref"))
  if (is.null(opt[[req]])) stop(sprintf("--%s is required.", req), call. = FALSE)
split_csv <- function(x) trimws(strsplit(x, ",")[[1]])

exp_names <- split_csv(opt$exposures); out_names <- split_csv(opt$outcomes)
if (length(exp_names) < 2) stop("MVMR needs at least 2 --exposures.", call. = FALSE)
exp_specs <- read_trait_manifest(split_csv(opt$exp_info), data_dir = opt$data_dir, select = exp_names, role = "exposure")
out_specs <- read_trait_manifest(split_csv(opt$out_info), data_dir = opt$data_dir, select = out_names, role = "outcome")
exp_specs <- exp_specs[intersect(gsub("[^A-Za-z0-9._+-]+", "_", exp_names), names(exp_specs))]  # preserve order
# MHC handling: exclude_mhc drops flagged SNPs at selection (per-trait region).
for (i in seq_along(exp_specs)) {
  if (!is.null(opt$mhc_region)) exp_specs[[i]]$mhc_region <- opt$mhc_region
}
set_label <- paste(vapply(exp_specs, `[[`, "", "name"), collapse = "_")

cat("\n=== MVMR model ===\n")
cat("exposures (joint):", paste(vapply(exp_specs, `[[`, "", "name"), collapse = " + "), "\n")
for (s in exp_specs) cat(sprintf("  %-12s type=%-10s ivs=%s\n    %s\n", s$name, s$type,
                                  ifelse(is.na(s$ivs), "(clump)", basename(s$ivs)), s$gwas))
cat("outcomes:", paste(vapply(out_specs, `[[`, "", "name"), collapse = ", "), "\n")
cat(sprintf("planned MVMR runs: %d (1 model x %d outcomes)\n\n", length(out_specs), length(out_specs)))
if (isTRUE(opt$dry_run)) { message("--dry_run set; not running MVMR."); quit(status = 0) }

out_par <- dirname(opt$out_prefix)
if (nzchar(out_par) && !dir.exists(out_par)) dir.create(out_par, recursive = TRUE, showWarnings = FALSE)

# ---- Pool + re-clump instruments once (shared across outcomes) --------------
pool <- pool_and_reclump(exp_specs, opt)
if (length(pool$pooled) < 2) stop("Fewer than 2 pooled instruments after re-clumping; cannot run MVMR.", call. = FALSE)

# Read + format each exposure at the pooled instrument set (once).
exp_formatted <- lapply(exp_specs, read_trait_at, keep_snps = pool$pooled, type = "exposure", opt = opt)
exposure_dat <- data.table::rbindlist(lapply(exp_formatted, data.table::as.data.table), fill = TRUE)

# Optional MHC exclusion: drop pooled SNPs in any exposure's build-specific MHC region.
if (isTRUE(opt$exclude_mhc)) {
  mhc_snps <- character(0)
  for (i in seq_along(exp_specs)) {
    fi <- data.table::as.data.table(exp_formatted[[i]])
    if (all(c("chr.exposure", "pos.exposure") %in% names(fi))) {
      reg <- parse_mhc_region(exp_specs[[i]]$mhc_region)
      hit <- fi[as.integer(chr.exposure) == reg$chr &
                as.numeric(pos.exposure) >= reg$start & as.numeric(pos.exposure) <= reg$end, SNP]
      mhc_snps <- union(mhc_snps, hit)
    } else {
      warning(sprintf("--exclude_mhc: %s has no chr/pos columns; its MHC SNPs cannot be flagged.", exp_specs[[i]]$name), call. = FALSE)
    }
  }
  if (length(mhc_snps)) {
    message(sprintf("[MVMR] --exclude_mhc: dropping %d MHC instrument(s).", length(mhc_snps)))
    pool$pooled  <- setdiff(pool$pooled, mhc_snps)
    exposure_dat <- exposure_dat[!(SNP %in% mhc_snps)]
    for (n in names(pool$per_exposure)) pool$per_exposure[[n]] <- setdiff(pool$per_exposure[[n]], mhc_snps)
  }
}

# ---- Run MVMR + univariate for each outcome ---------------------------------
all_res <- list()
for (o in out_specs) {
  message(sprintf("\n=== [%s] MVMR: %s -> %s ===", format(Sys.time(), "%T"), set_label, o$name))
  outcome_dat <- tryCatch(read_trait_at(o, pool$pooled, "outcome", opt),
                          error = function(e) { warning(sprintf("[%s] outcome read failed: %s", o$name, e$message), call. = FALSE); NULL })
  if (is.null(outcome_dat)) next
  mvdat <- tryCatch(TwoSampleMR::mv_harmonise_data(exposure_dat = as.data.frame(exposure_dat),
                                                   outcome_dat = as.data.frame(outcome_dat), harmonise_strictness = 2),
                    error = function(e) { warning(sprintf("[%s] mv_harmonise failed: %s", o$name, e$message), call. = FALSE); NULL })
  if (is.null(mvdat) || nrow(mvdat$exposure_beta) < 2) { warning(sprintf("[%s] <2 harmonised instruments; skipping.", o$name), call. = FALSE); next }

  mv <- tryCatch(run_mvmr(mvdat, o, opt),
                 error = function(e) { warning(sprintf("[%s] MVMR failed: %s", o$name, e$message), call. = FALSE); NULL })
  res_o <- list(mv)

  if (!isTRUE(opt$no_univariate)) {
    for (e in exp_specs) {
      uv <- tryCatch(univariate_ivw(e, o, pool$per_exposure[[e$name]], opt),
                     error = function(err) { warning(sprintf("[%s univariable %s] failed: %s", o$name, e$name, err$message), call. = FALSE); NULL })
      res_o[[length(res_o) + 1]] <- uv
    }
  }
  res_o <- data.table::rbindlist(Filter(Negate(is.null), res_o), fill = TRUE)
  if (nrow(res_o) == 0) next
  res_o[, exposure_set := set_label]
  data.table::fwrite(res_o, sprintf("%s%s_vs_%s_mvmr.csv", opt$out_prefix, set_label, o$name), na = "NA")
  all_res[[length(all_res) + 1]] <- res_o
}

if (length(all_res) == 0) { message("No MVMR results produced."); quit(status = 0) }
combined <- data.table::rbindlist(all_res, fill = TRUE)
data.table::fwrite(combined, sprintf("%s%s_all_mvmr_results.csv", opt$out_prefix, set_label), na = "NA")
message(sprintf("Wrote combined results: %s%s_all_mvmr_results.csv", opt$out_prefix, set_label))

# ---- Univariable-vs-multivariable forest, one panel per outcome -------------
if (!isTRUE(opt$no_plot) && requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  p <- combined[!is.na(b)]
  p[, is_or := !is.na(or)]
  p[, est := fifelse(is_or, or, b)]
  p[, lo  := fifelse(is_or, or_lci95, lo_ci)]
  p[, hi  := fifelse(is_or, or_uci95, up_ci)]
  p[, analysis := factor(analysis, levels = c("univariable", "multivariable"))]
  p[, row := paste(exposure, analysis)]
  gg <- ggplot(p, aes(est, exposure, colour = analysis)) +
    geom_vline(data = unique(p[, .(outcome, is_or)]),
               aes(xintercept = ifelse(is_or, 1, 0)), linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.3, position = position_dodge(width = 0.6)) +
    geom_point(size = 2.6, position = position_dodge(width = 0.6)) +
    facet_wrap(~ outcome, scales = "free_x", ncol = 1) +
    labs(title = sprintf("MVMR: %s", gsub("_", " + ", set_label)),
         x = "Estimate (OR for binary, β for continuous; 95% CI)", y = NULL, colour = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  ggsave(sprintf("%s%s_mvmr_forest.png", opt$out_prefix, set_label), gg,
         width = 8, height = 1.6 + 1.4 * length(out_specs) * length(exp_specs), dpi = 300)
  message(sprintf("Wrote forest: %s%s_mvmr_forest.png", opt$out_prefix, set_label))
}

message("MVMR complete.")
