#!/usr/bin/env Rscript
# ===========================================================================
# Forest plots for the MR grid.
# Reads *_full_mr_results.csv, tags each with its config (from the filename
# prefix) and phenotype, PRINTS what it found (so you can see stale/duplicate
# files), then draws one clean 2-panel forest per phenotype with the four
# configs colour-coded and separated. Edit `out_dir` and `keep`, then:
#   Rscript plot_forest.R
# ===========================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

out_dir  <- "/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Ana/MR_pipeline/output/eos"  # <-- edit
plot_dir <- file.path(out_dir, "plots"); dir.create(plot_dir, showWarnings = FALSE)

# `anchor` is the trait that is common to every pair in this grid (the shared
# exposure/outcome). One forest is drawn per OTHER trait, so a grid of
# EO vs {asthma, POI, ANM} anchored on "EO" gives three 2-panel forests.
# For the old menopause<->asthma grid, set anchor <- "asthma".
anchor <- "EO"

# Only plot files whose name matches this (drops stale earlier runs). Default
# keeps every file that mentions the anchor trait.
keep <- anchor

# ---- Read + tidy -----------------------------------------------------------
files <- list.files(out_dir, pattern = "_full_mr_results\\.csv$", full.names = TRUE)
if (nzchar(keep)) files <- files[grepl(keep, basename(files))]
stopifnot(length(files) > 0)

read_one <- function(f) {
  # NB: use `fn` not `b` for the filename - the results table has a column `b`
  # (beta), which would shadow a local `b` inside the data.table j-expression.
  d <- fread(f); fn <- basename(f)
  d[, config := fifelse(grepl("noMHC", fn) & grepl("rigid", fn), "strict, no MHC",
                 fifelse(grepl("rigid", fn), "strict, +MHC",
                 fifelse(grepl("noMHC", fn), "lenient, no MHC", "lenient, +MHC")))]
  d[, file := fn]
  d
}
res <- rbindlist(lapply(files, read_one), fill = TRUE)
for (col in c("or","or_lci95","or_uci95","lo_ci","up_ci","b","se"))
  if (!col %in% names(res)) res[, (col) := NA_real_]

# unified estimate/CI + phenotype/direction
res[, is_or := !is.na(or)]
res[, `:=`(est = fifelse(is_or, or, b), lo = fifelse(is_or, or_lci95, lo_ci),
           hi = fifelse(is_or, or_uci95, up_ci), scale = fifelse(is_or, "OR", "beta"))]
# Anchor on `anchor`: the "phenotype" is whichever trait is NOT the anchor,
# and the direction is written relative to the anchor. Strip the "_ofh_meta"
# suffix for clean labels (POI_ofh_meta -> POI, ANM_ofh_meta -> ANM).
clean <- function(x) gsub("_ofh_meta$", "", x)
res <- res[exposure == anchor | outcome == anchor]   # keep only pairs involving the anchor
res[, phenotype := clean(fifelse(exposure == anchor, outcome, exposure))]
res[, direction := fifelse(exposure == anchor,
                           paste(anchor, "→", phenotype),
                           paste(phenotype, "→", anchor))]

method_map <- c("Inverse variance weighted"="IVW","Weighted median"="Weighted median",
                "MR Egger"="MR-Egger","mr_ivw_steiger"="IVW-Steiger",
                "mr_presso_raw"="MR-PRESSO (raw)","mr_presso_corrected"="MR-PRESSO (corrected)")
res <- res[method %in% names(method_map)]
res[, method := factor(method_map[method], levels = rev(unname(method_map)))]
res[, config := factor(config, levels = c("lenient, +MHC","strict, +MHC",
                                          "lenient, no MHC","strict, no MHC"))]

# ---- Transparency: show exactly what is being plotted ----------------------
cat("\nFiles being plotted (should be ONE per phenotype x direction x config):\n")
print(unique(res[, .(file, phenotype, exposure, outcome, config)]), nrow = Inf)
dupes <- res[method == "IVW", .N, by = .(phenotype, direction, config)][N > 1]
if (nrow(dupes)) { cat("\nWARNING: duplicate config/direction combos (stale files?):\n"); print(dupes) }

# ---- One 2-panel forest per phenotype --------------------------------------
forest_pheno <- function(d, ph) {
  ggplot(d, aes(x = est, y = method, colour = config)) +
    geom_vline(data = unique(d[, .(direction, scale)]),
               aes(xintercept = ifelse(scale == "OR", 1, 0)),
               linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.3,
                   position = position_dodge(width = 0.7)) +
    geom_point(position = position_dodge(width = 0.7), size = 2.2) +
    facet_wrap(~ direction, scales = "free_x") +
    labs(title = paste("MR:", anchor, "vs", ph), x = "Estimate (OR or beta; 95% CI)",
         y = NULL, colour = "Instruments") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

for (ph in unique(res$phenotype)) {
  d <- res[phenotype == ph]
  ggsave(file.path(plot_dir, sprintf("forest_%s.png", ph)),
         forest_pheno(d, ph), width = 11, height = 5, dpi = 300)
}

cat(sprintf("\nWrote %d per-phenotype forest plots to %s\n", length(unique(res$phenotype)), plot_dir))
