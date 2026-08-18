# ===========================================================================
# Multivariable MR (MVMR) engine.
#
# Reuses the univariate pipeline's data layer (read_trait_manifest, read_gwas,
# select_instruments, read_iv_list, format_gwas, resolve_plink from
# pipeline_functions.R) and adds the MVMR-specific steps:
#   * pool_and_reclump()  - union the exposures' instruments, then re-clump the
#     union to independence (so the COMBINED instrument set is mutually
#     independent, which MVMR-IVW assumes).
#   * run_mvmr()          - MVMR-IVW (TwoSampleMR::mv_multiple) + conditional
#     F-statistics / Q-pleiotropy / qhet (MVMR package) + optional MVMR-Egger
#     (MendelianRandomization).
#   * univariate_ivw()    - the plain IVW estimate for one exposure->outcome, for
#     the univariable-vs-multivariable comparison.
# Requires pipeline_functions.R to be sourced first.
# ===========================================================================

#' Instrument SNP IDs for one exposure: its ivs.file (authoritative) if present,
#' else clump the GWAS at clump_p / clump_r2 / clump_kb (F>=f_stat).
mvmr_instrument_snps <- function(exp_spec, opt) {
  if (!is.na(exp_spec$ivs) && nzchar(exp_spec$ivs) && file.exists(exp_spec$ivs)) {
    ivs <- read_iv_list(exp_spec$ivs)
    message(sprintf("[MVMR] %s: %d instruments from %s.", exp_spec$name, length(ivs), basename(exp_spec$ivs)))
    return(ivs)
  }
  raw <- read_gwas(exp_spec$gwas, type = "exposure", col_args = exp_spec$col_args,
                   trait_name = exp_spec$name, tmp_dir = opt$tmp_dir, pval_thresh = opt$clump_p,
                   n_total = exp_spec$n_total, ncase_total = exp_spec$ncase_total, cache_dir = opt$cache_dir)
  sel <- select_instruments(raw, exp_spec$name, clump_p = opt$clump_p, clump_kb = opt$clump_kb,
                            clump_r2 = opt$clump_r2, ld_ref = opt$ld_ref, plink_bin = opt$plink_bin,
                            min_f_stat = opt$f_stat, skip_clump = FALSE, iv_list = NULL,
                            mhc_region = exp_spec$mhc_region, exclude_mhc = FALSE)
  message(sprintf("[MVMR] %s: %d instruments (clumped r2<%g).", exp_spec$name, nrow(sel), opt$clump_r2))
  sel$SNP
}

#' Pool the exposures' instruments and re-clump the union to independence.
#' Each union SNP is tagged with its MINIMUM p across the exposures so the
#' re-clump keeps the most-significant lead per LD block (deterministic).
#' @return list(pooled = final SNP vector, per_exposure = named list of SNP vecs).
pool_and_reclump <- function(exp_specs, opt) {
  per <- lapply(exp_specs, mvmr_instrument_snps, opt = opt)
  names(per) <- vapply(exp_specs, `[[`, "", "name")
  union_snps <- unique(unlist(per, use.names = FALSE))
  if (length(union_snps) == 0) stop("MVMR: no instruments for any exposure.", call. = FALSE)

  # minimum p across exposures at each union SNP (read each exposure once; cached)
  minp <- setNames(rep(Inf, length(union_snps)), union_snps)
  for (e in exp_specs) {
    raw <- data.table::as.data.table(read_gwas(
      e$gwas, type = "exposure", col_args = e$col_args, trait_name = e$name, tmp_dir = opt$tmp_dir,
      keep_snps = union_snps, n_total = e$n_total, ncase_total = e$ncase_total, cache_dir = opt$cache_dir))
    pv <- suppressWarnings(as.numeric(raw$pval)); names(pv) <- raw$SNP
    common <- intersect(names(pv), union_snps)
    minp[common] <- pmin(minp[common], pv[common], na.rm = TRUE)
  }
  minp[!is.finite(minp)] <- 1

  plink <- resolve_plink(opt$plink_bin)
  clumped <- ieugwasr::ld_clump(
    dplyr::tibble(rsid = union_snps, pval = as.numeric(minp[union_snps])),
    plink_bin = plink, bfile = opt$ld_ref,
    clump_r2 = opt$clump_r2, clump_kb = opt$clump_kb, clump_p = 1)
  pooled <- clumped$rsid
  message(sprintf("[MVMR] pooled instruments: %d union -> %d after re-clump (r2<%g, %dkb).",
                  length(union_snps), length(pooled), opt$clump_r2, opt$clump_kb))
  list(pooled = pooled, per_exposure = per)
}

#' Read + format one exposure (or outcome) restricted to `keep_snps`.
read_trait_at <- function(spec, keep_snps, type, opt) {
  raw <- read_gwas(spec$gwas, type = type, col_args = spec$col_args, trait_name = spec$name,
                   tmp_dir = opt$tmp_dir, keep_snps = keep_snps,
                   n_total = spec$n_total, ncase_total = spec$ncase_total, cache_dir = opt$cache_dir)
  format_gwas(raw, type = type, trait_name = spec$name)
}

#' Multivariable MR for one harmonised mvdat against one outcome.
#' @return data.table, one row per exposure (conditional effect + F/Q/qhet [+ egger]).
run_mvmr <- function(mvdat, out_spec, opt) {
  outcome_is_binary <- !identical(tolower(trimws(as.character(out_spec$type))), "continuous")
  # exposure names in the column order of the harmonised beta matrix
  exp_order <- colnames(mvdat$exposure_beta)
  exp_names <- mvdat$expname$exposure[match(exp_order, mvdat$expname$id.exposure)]

  res <- TwoSampleMR::mv_multiple(mvdat, pval_threshold = opt$mv_pval_threshold, plots = FALSE)
  dt <- data.table::as.data.table(res$result)[, .(exposure, nsnp, b, se, pval)]
  dt[, outcome := out_spec$name]

  # MVMR package: conditional F (instrument strength), Q (pleiotropy), qhet (robust estimate)
  if (requireNamespace("MVMR", quietly = TRUE)) {
    Fdat <- MVMR::format_mvmr(BXGs = as.data.frame(mvdat$exposure_beta),
                              BYG  = mvdat$outcome_beta,
                              seBXGs = as.data.frame(mvdat$exposure_se),
                              seBYG  = mvdat$outcome_se,
                              RSID = rownames(mvdat$exposure_beta))
    sres <- tryCatch(MVMR::strength_mvmr(Fdat, gencov = 0), error = function(e) NULL)
    pres <- tryCatch(MVMR::pleiotropy_mvmr(Fdat, gencov = 0), error = function(e) NULL)
    qh   <- tryCatch(MVMR::qhet_mvmr(Fdat, diag(length(exp_names)), CI = FALSE, iterations = 100),
                     error = function(e) NULL)
    if (!is.null(sres)) { f <- setNames(as.numeric(sres[1, ]), exp_names); dt[, cond_Fstat := f[exposure]] }
    if (!is.null(pres)) { dt[, Qstat := as.numeric(pres$Qstat)]; dt[, Qpval := as.numeric(pres$Qpval)] }
    if (!is.null(qh))   { qcol <- if ("Effect Estimates" %in% names(qh)) qh[["Effect Estimates"]] else qh[[1]]
                          qv <- as.numeric(qcol)
                          if (length(qv) == length(exp_names)) { qn <- setNames(qv, exp_names); dt[, b_qhet := qn[exposure]] } }
  } else {
    warning("Package 'MVMR' not installed; conditional F / Q / qhet omitted.", call. = FALSE)
  }

  # Optional MVMR-Egger (intercept = directional pleiotropy across exposures)
  if (isTRUE(opt$egger) && requireNamespace("MendelianRandomization", quietly = TRUE)) {
    eg <- tryCatch({
      mvin <- MendelianRandomization::mr_mvinput(
        bx = as.matrix(mvdat$exposure_beta), bxse = as.matrix(mvdat$exposure_se),
        by = as.numeric(mvdat$outcome_beta), byse = as.numeric(mvdat$outcome_se),
        snps = rownames(mvdat$exposure_beta))
      MendelianRandomization::mr_mvegger(mvin, orientate = 1)
    }, error = function(e) NULL)
    if (!is.null(eg)) {
      ebeta <- setNames(as.numeric(eg$Estimate), exp_names)
      dt[, egger_b := ebeta[exposure]]
      dt[, egger_intercept := as.numeric(eg$Intercept)]
      dt[, egger_intercept_pval := as.numeric(eg$Pvalue.Int)]
    }
  }

  if (outcome_is_binary) {
    dt[, `:=`(or = exp(b), or_lci95 = exp(b - 1.96 * se), or_uci95 = exp(b + 1.96 * se))]
  } else {
    dt[, `:=`(lo_ci = b - 1.96 * se, up_ci = b + 1.96 * se)]
  }
  dt[, analysis := "multivariable"]
  dt[]
}

#' Univariate IVW for one exposure->outcome, using the exposure's OWN instruments.
univariate_ivw <- function(exp_spec, out_spec, ivs, opt) {
  ivs <- unique(ivs); if (length(ivs) == 0) return(NULL)
  exp_f <- tryCatch(read_trait_at(exp_spec, ivs, "exposure", opt), error = function(e) NULL)
  out_f <- tryCatch(read_trait_at(out_spec, ivs, "outcome",  opt), error = function(e) NULL)
  if (is.null(exp_f) || is.null(out_f)) return(NULL)
  h <- tryCatch(TwoSampleMR::harmonise_data(exp_f, out_f, action = 2), error = function(e) NULL)
  if (is.null(h)) return(NULL)
  if ("mr_keep" %in% names(h)) h <- h[h$mr_keep, , drop = FALSE]
  if (nrow(h) < 2) return(NULL)
  r <- tryCatch(TwoSampleMR::mr(h, method_list = "mr_ivw"), error = function(e) NULL)
  if (is.null(r) || nrow(r) == 0) return(NULL)
  outcome_is_binary <- !identical(tolower(trimws(as.character(out_spec$type))), "continuous")
  dt <- data.table::data.table(exposure = exp_spec$name, outcome = out_spec$name,
                               nsnp = r$nsnp[1], b = r$b[1], se = r$se[1], pval = r$pval[1],
                               analysis = "univariable")
  if (outcome_is_binary) dt[, `:=`(or = exp(b), or_lci95 = exp(b - 1.96 * se), or_uci95 = exp(b + 1.96 * se))]
  else                   dt[, `:=`(lo_ci = b - 1.96 * se, up_ci = b + 1.96 * se)]
  dt
}
