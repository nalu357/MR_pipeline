#!/usr/bin/env Rscript
# ===========================================================================
# Forest plots for the MR grid.
# Reads every *_full_mr_results.csv in an output directory, tidies them
# (config from filename prefix; direction/phenotype from the exposure/outcome
# columns; OR for a binary outcome, beta for a continuous one) and draws a few
# forest plots. Edit `out_dir` and run:  Rscript plot_forest.R
# ===========================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

out_dir  <- "/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline/output"  # <-- edit
plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE)

# ---- 1. Read + tidy all result files ---------------------------------------
files <- list.files(out_dir, pattern = "_full_mr_results\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0)

read_one <- function(f) {
  d <- fread(f)
  b <- basename(f)
  d[, config := fifelse(grepl("noMHC", b) & grepl("rigid", b), "strict, no MHC",
                 fifelse(grepl("rigid", b), "strict, +MHC",
                 fifelse(grepl("noMHC", b), "lenient, no MHC", "lenient, +MHC")))]
  d
}
res <- rbindlist(lapply(files, read_one), fill = TRUE)

# Ensure OR/CI columns exist even if all files were continuous (or vice-versa).
for (col in c("or", "or_lci95", "or_uci95", "lo_ci", "up_ci", "b", "se"))
  if (!col %in% names(res)) res[, (col) := NA_real_]

# Unified estimate/CI: OR for binary outcomes, beta for continuous.
res[, is_or := !is.na(or)]
res[, `:=`(est   = fifelse(is_or, or,       b),
           lo    = fifelse(is_or, or_lci95, lo_ci),
           hi    = fifelse(is_or, or_uci95, up_ci),
           scale = fifelse(is_or, "OR", "beta"))]
res[, direction := paste(exposure, "→", outcome)]

# Nice method labels + display order (top -> bottom).
method_map <- c("Inverse variance weighted" = "IVW",
                "Weighted median"           = "Weighted median",
                "MR Egger"                  = "MR-Egger",
                "mr_ivw_steiger"            = "IVW-Steiger",
                "mr_presso_raw"             = "MR-PRESSO (raw)",
                "mr_presso_corrected"       = "MR-PRESSO (corrected)")
res <- res[method %in% names(method_map)]
res[, method := factor(method_map[method], levels = rev(unname(method_map)))]
res[, config := factor(config, levels = c("lenient, +MHC", "strict, +MHC",
                                          "lenient, no MHC", "strict, no MHC"))]

# ---- 2. A reusable forest-plot helper --------------------------------------
forest <- function(d, yvar, title, xlab, null_x, logx = FALSE) {
  p <- ggplot(d, aes(x = est, y = .data[[yvar]], colour = config)) +
    geom_vline(xintercept = null_x, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25,
                   position = position_dodge(width = 0.65)) +
    geom_point(position = position_dodge(width = 0.65), size = 2.2) +
    labs(title = title, x = xlab, y = NULL, colour = "Instruments") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  if (logx) p <- p + scale_x_log10()
  p
}

# ---- 3. Plots ---------------------------------------------------------------

# (a) The headline: asthma -> POI, all methods x config (OR, log scale)
d_poi <- res[exposure == "asthma" & grepl("POI", outcome)]
if (nrow(d_poi)) {
  ggsave(file.path(plot_dir, "forest_asthma_to_POI.png"),
         forest(d_poi, "method", "Asthma → POI",
                "OR of POI (95% CI) per unit asthma liability", null_x = 1, logx = TRUE),
         width = 8, height = 5, dpi = 300)
}

# (b) Menopause traits -> asthma, IVW only (OR); one row per exposure x config
d_rev <- res[outcome == "asthma" & method == "IVW"]
if (nrow(d_rev)) {
  ggsave(file.path(plot_dir, "forest_meno_to_asthma_IVW.png"),
         forest(d_rev, "exposure", "Menopause traits → asthma (IVW)",
                "OR of asthma (95% CI)", null_x = 1, logx = TRUE),
         width = 8, height = 4.5, dpi = 300)
}

# (c) asthma -> continuous menopause (ANM/ReproGen), IVW only (beta)
d_fwd_cont <- res[exposure == "asthma" & scale == "beta" & method == "IVW"]
if (nrow(d_fwd_cont)) {
  ggsave(file.path(plot_dir, "forest_asthma_to_menopause_beta_IVW.png"),
         forest(d_fwd_cont, "outcome", "Asthma → age at menopause (IVW)",
                "beta (95% CI), years per unit asthma liability", null_x = 0),
         width = 8, height = 4, dpi = 300)
}

# (d) Everything, faceted: one panel per direction, methods on y, colour = config
p_all <- forest(res, "method", "All MR analyses", "Estimate (OR or beta; 95% CI)",
                null_x = NA) +
  geom_vline(data = unique(res[, .(direction, scale)]),
             aes(xintercept = ifelse(scale == "OR", 1, 0)),
             linetype = "dashed", colour = "grey50") +
  facet_wrap(~ direction, scales = "free_x", ncol = 2)
ggsave(file.path(plot_dir, "forest_all_analyses.png"), p_all,
       width = 12, height = 0.6 * length(unique(res$direction)) + 4, dpi = 300)

cat(sprintf("Read %d result files. Plots written to %s\n", length(files), plot_dir))
