#!/usr/bin/env Rscript

# ===========================================================================
# MR_pipeline dependency installer
# ---------------------------------------------------------------------------
# Installs all packages the pipeline needs into a repo-local `Rpackages/`
# library. Run this ONCE on the machine where you will run the pipeline
# (e.g. the cluster), with the R version you intend to use:
#
#     Rscript install_dependencies.R
#
# mr_pipeline.R auto-detects `Rpackages/` next to the repo, so after this
# finishes you can run the pipeline with no --lib_path needed.
#
# Reproducibility: fill the PINS block below with the exact versions/SHAs you
# validated against (from the version-dump command in the README). Any entry
# left as "" installs the latest available version.
# ===========================================================================

# ---- Locate the repo-local library (next to this script) ------------------
get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) == 1) return(dirname(normalizePath(sub("^--file=", "", f))))
  getwd()
}
script_dir <- get_script_dir()
lib <- file.path(script_dir, "Rpackages")
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))

repos <- "https://cloud.r-project.org"

cat("=====================================================\n")
cat("MR_pipeline dependency installer\n")
cat("  Library : ", lib, "\n")
cat("  R       : ", R.version.string, "\n")
cat("=====================================================\n\n")

# ===========================================================================
# PINS  ---  fill with your validated versions/SHAs, or leave "" for latest.
# ===========================================================================
# CRAN packages: pin an exact version string, e.g. "1.15.4".
cran_pins <- c(
  remotes    = "",   # bootstrap installer (used below)
  optparse   = "",
  dplyr      = "",
  data.table = ""
)
# GitHub packages: value is the git ref (SHA or tag), e.g. "0d1e2f3" or "v0.6.3".
github_pins <- c(
  "MRCIEU/TwoSampleMR"        = "",
  "MRCIEU/ieugwasr"           = "",
  "MRCIEU/genetics.binaRies"  = "",  # optional PLINK-binary helper
  "rondolab/MR-PRESSO"        = ""
)
# ===========================================================================

installed_in_lib <- function(pkg) pkg %in% rownames(installed.packages(lib.loc = lib))

install_cran <- function(pkg, version) {
  if (installed_in_lib(pkg) &&
      (version == "" || as.character(packageVersion(pkg)) == version)) {
    cat(sprintf("[skip ] %-14s already in Rpackages/ (%s)\n", pkg, packageVersion(pkg)))
    return(invisible(TRUE))
  }
  cat(sprintf("[cran ] %-14s %s\n", pkg, if (version == "") "(latest)" else version))
  if (version == "") {
    install.packages(pkg, lib = lib, repos = repos)
  } else {
    remotes::install_version(pkg, version = version, lib = lib, repos = repos, upgrade = "never")
  }
}

install_gh <- function(repo, ref) {
  pkg <- sub("^.*/", "", repo)
  pkg <- gsub("MR-PRESSO", "MRPRESSO", pkg)  # repo name != package name
  spec <- if (ref == "") repo else paste0(repo, "@", ref)
  if (installed_in_lib(pkg) && ref == "") {
    cat(sprintf("[skip ] %-14s already in Rpackages/ (%s)\n", pkg, packageVersion(pkg)))
    return(invisible(TRUE))
  }
  cat(sprintf("[github] %-13s %s\n", pkg, spec))
  remotes::install_github(spec, lib = lib, repos = repos, upgrade = "never")
}

# ---- Bootstrap: remotes must exist before install_version/install_github --
if (!requireNamespace("remotes", quietly = TRUE)) {
  cat("[cran ] remotes        (bootstrap)\n")
  install.packages("remotes", lib = lib, repos = repos)
}

# ---- Install CRAN, then GitHub packages -----------------------------------
for (p in names(cran_pins)) install_cran(p, cran_pins[[p]])
for (r in names(github_pins)) install_gh(r, github_pins[[r]])

# ---- Verify ---------------------------------------------------------------
cat("\n----- Verifying installation -----\n")
required <- c("optparse", "dplyr", "data.table", "TwoSampleMR", "ieugwasr", "MRPRESSO")
optional <- c("genetics.binaRies")
ok <- TRUE
for (p in c(required, optional)) {
  present <- requireNamespace(p, quietly = TRUE, lib.loc = lib)
  if (present) {
    cat(sprintf("  [ok ] %-16s %s\n", p, packageVersion(p)))
  } else {
    tag <- if (p %in% optional) "[opt]" else "[MISS]"
    cat(sprintf("  %s %-16s not installed\n", tag, p))
    if (!(p %in% optional)) ok <- FALSE
  }
}
cat("-----------------------------------\n")
if (!ok) {
  stop("One or more REQUIRED packages failed to install. See messages above.", call. = FALSE)
}
cat("All required packages installed into Rpackages/. You can now run mr_pipeline.R.\n")
