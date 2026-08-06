#!/usr/bin/env Rscript
# ===========================================================================
# Per-SNP and diagnostic MR plots, from the harmonised data the pipeline saved
# (<prefix><exp>_vs_<out>_harmonized_data.rds). For each analysis it produces
# the four standard TwoSampleMR diagnostics:
#   * scatter        - exposure vs outcome effect, with method slopes (the main
#                      visual of the causal estimate + outliers)
#   * singlesnp      - per-SNP Wald ratios + meta estimate (which SNPs drive it)
#   * leave-one-out  - IVW dropping each SNP in turn (is one SNP responsible?)
#   * funnel         - single-SNP estimate vs precision (asymmetry = pleiotropy)
#
# Edit `out_dir` and `pattern`, then:  Rscript plot_mr_diagnostics.R
# ===========================================================================
suppressPackageStartupMessages({ library(TwoSampleMR); library(ggplot2) })

out_dir  <- "/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline/output"  # <-- edit
plot_dir <- file.path(out_dir, "plots"); dir.create(plot_dir, showWarnings = FALSE)

# Which harmonised files to plot. "" = all; or narrow it, e.g. "asthma_vs_POI".
# Which analyses to plot:
#   pheno_re : phenotypes to include (regex; matches the filename)
#   mhc_re   : "noMHC" = MHC-excluded runs only; set to "" for all MHC settings
pheno_re <- "POI_ofh|ANM_ofh"
mhc_re   <- "noMHC"

files <- list.files(out_dir, pattern = "_harmonized_data\\.rds$", full.names = TRUE)
files <- files[grepl(pheno_re, basename(files))]
if (nzchar(mhc_re)) files <- files[grepl(mhc_re, basename(files))]
stopifnot(length(files) > 0)
cat(sprintf("Plotting diagnostics for %d analyses:\n  %s\n",
            length(files), paste(basename(files), collapse = "\n  ")))

# Per-SNP forest & leave-one-out get one row per SNP, so they are unreadable
# (and exceed ggsave's 50-inch cap) above a few hundred SNPs. Skip those two
# for large instrument sets; scatter + funnel are fixed-size and fine for any N.
max_forest_snp <- 150

diagnostics <- function(rds) {
  stub <- file.path(plot_dir, sub("_harmonized_data\\.rds$", "", basename(rds)))
  dat  <- readRDS(rds)
  if ("mr_keep" %in% names(dat)) dat <- dat[dat$mr_keep == TRUE, , drop = FALSE]
  n <- nrow(dat)
  if (n < 3) { message(sprintf("Skipping %s (only %d SNPs).", basename(rds), n)); return(invisible()) }
  h <- min(0.13 * n + 2.5, 45)   # per-SNP-row height, capped under the 50-inch limit

  mr_res <- tryCatch(mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")),
                     error = function(e) NULL)
  ss  <- tryCatch(mr_singlesnp(dat),   error = function(e) NULL)
  loo <- tryCatch(mr_leaveoneout(dat), error = function(e) NULL)

  # Always: scatter + funnel (fixed size).
  if (!is.null(mr_res)) ggsave(paste0(stub, "_scatter.png"),
       mr_scatter_plot(mr_res, dat)[[1]], width = 7, height = 6, dpi = 300)
  if (!is.null(ss)) ggsave(paste0(stub, "_funnel.png"),
       mr_funnel_plot(ss)[[1]], width = 7, height = 6, dpi = 300)

  # Per-SNP forest + leave-one-out only when the instrument set is small enough.
  if (n <= max_forest_snp) {
    if (!is.null(ss))  ggsave(paste0(stub, "_singlesnp_forest.png"),
         mr_forest_plot(ss)[[1]], width = 7, height = h, dpi = 300, limitsize = FALSE)
    if (!is.null(loo)) ggsave(paste0(stub, "_leaveoneout.png"),
         mr_leaveoneout_plot(loo)[[1]], width = 7, height = h, dpi = 300, limitsize = FALSE)
  } else {
    message(sprintf("  %s: %d SNPs > %d - scatter+funnel only (per-SNP forest & leave-one-out skipped as unreadable).",
                    basename(rds), n, max_forest_snp))
  }
  message(sprintf("Done: %s (%d SNPs)", basename(rds), n))
}

invisible(lapply(files, diagnostics))
cat(sprintf("Diagnostics for %d analyses written to %s\n", length(files), plot_dir))
