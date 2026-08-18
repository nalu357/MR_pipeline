# Mendelian Randomization Pipeline

## Purpose

This R-based pipeline performs a standard two-sample Mendelian Randomization (MR) analysis. It takes GWAS summary statistics for an exposure and an outcome, selects genetic instruments for the exposure using LD clumping, harmonizes the data, runs various MR methods, performs sensitivity analyses, and outputs the results. It is designed to be flexible and efficient, supporting the following scenarios:

- **Single exposure, single outcome:** Analyze one exposure-outcome pair.
- **Multiple exposures, multiple outcomes:** Use a shell loop to analyze all exposures in a directory against all outcomes in another directory.
## Dependencies

### R Packages

The pipeline needs these R packages:

*   `optparse` (CRAN): command-line argument parsing.
*   `dplyr` (CRAN): data manipulation.
*   `data.table` (CRAN): fast file reading (`fread`) and manipulation.
*   `TwoSampleMR` (GitHub, MRCIEU): core MR analysis functions.
*   `ieugwasr` (GitHub, MRCIEU): LD clumping via PLINK.
*   `genetics.binaRies` (GitHub, MRCIEU): optional helper to locate PLINK binaries if PLINK is not on PATH.
*   `MRPRESSO` (GitHub, rondolab): required only if `--presso` is used.

#### Recommended: one-step install into a repo-local library

Run the installer once on the machine where you will run the pipeline (use the R version you intend to run with):

```bash
Rscript install_dependencies.R
```

This installs everything into a repo-local `Rpackages/` directory (git-ignored). `mr_pipeline.R` **auto-detects `Rpackages/`**, so no `--lib_path` is needed afterwards. Precedence for the library path is: `--lib_path` (explicit) > repo-local `Rpackages/` > R's default.

#### Pinning exact versions (reproducibility)

By default `install_dependencies.R` installs the latest versions. To lock the versions you validated against, first dump them on the target machine:

```r
pkgs <- c("optparse","dplyr","data.table","TwoSampleMR","ieugwasr","genetics.binaRies","MRPRESSO")
cat("R:", R.version.string, "\n")
for (p in pkgs) {
  if (requireNamespace(p, quietly=TRUE)) {
    d <- packageDescription(p); sha <- d$RemoteSha; if (is.null(sha)) sha <- d$GithubSHA1
    cat(sprintf("%-18s v=%-10s sha=%s\n", p, d$Version, if (is.null(sha)) "CRAN" else sha))
  } else cat(sprintf("%-18s NOT INSTALLED\n", p))
}
```

Then fill the `cran_pins` (exact versions) and `github_pins` (git SHAs/tags) blocks near the top of `install_dependencies.R` and re-run it.

### External Software

*   **PLINK (v1.9 or higher):** Required for local LD clumping. The pipeline uses `ieugwasr::ld_clump` which calls PLINK. Ensure PLINK is installed and either:
    *   Its location is in your system's PATH environment variable, OR
    *   You provide the path to the PLINK executable using the `--plink_bin` argument.

## Input GWAS File Requirements

The pipeline requires two input GWAS summary statistics files: one for the exposure and one for the outcome.

*   **Format:** Files should be plain text (e.g., `.tsv`, `.csv`, `.txt`) and can be compressed with gzip (`.gz`). They must be readable by `data.table::fread`. Ensure the file uses a consistent delimiter (tab, comma, space).
*   **Headers:** The script expects the files to have header rows. The exact names of the columns can vary, as you will map them using command-line arguments.
*   **Essential Columns:** The following information *must* be present in both files for the core analysis to run:
    *   **SNP Identifier:** Usually an rsID. Map using `--exp_snp` / `--out_snp`.
    *   **Effect Allele (EA):** The allele associated with the effect size. Must be A, T, C, or G. Map using `--exp_ea` / `--out_ea`.
    *   **Non-Effect Allele (NEA):** The other allele. Must be A, T, C, or G. Map using `--exp_nea` / `--out_nea`.
    *   **Effect Size (Beta):** The regression coefficient for the EA. Map using `--exp_beta` / `--out_beta`.
    *   **Standard Error (SE):** The standard error of the Beta. Map using `--exp_se` / `--out_se`.
    *   **P-value (P):** The p-value for the association. Map using `--exp_p` / `--out_p`.
*   **Recommended Columns:** These columns are not strictly essential for basic MR but are needed for specific analyses or improve harmonization:
    *   **Effect Allele Frequency (EAF):** Highly recommended for harmonization (allele checks, handling ambiguous SNPs). Map using `--exp_eaf` / `--out_eaf`.
    *   **Sample Size (N):** Used for Steiger filtering (`--steiger`). Map a per-SNP column with `--exp_n` / `--out_n`, or, if the file has no N column, supply the study's total sample size with `--exp_n_total` / `--out_n_total` (applied as a constant N to every SNP). (Note: the F-statistic is computed as beta²/se² and does not need N.)
    *   **Chromosome (CHR):** Useful for some functions and debugging. Map using `--exp_chr` / `--out_chr`.
    *   **Position (POS):** Useful for some functions and debugging. Map using `--exp_pos` / `--out_pos`.
    *   **Number of Cases (NCASE):** Required for binary outcomes (`--out_type binary`) if Steiger filtering (`--steiger`) is enabled. Map using `--out_ncase`.
*   **Data Cleaning:** The script performs basic cleaning:
    *   Removes rows with missing values in essential columns (SNP, Beta, SE, EA, NEA, P).
    *   Ensures Effect/Non-Effect Alleles are A, T, C, or G (uppercase).
    *   Converts relevant columns (Beta, SE, P, EAF, N, etc.) to numeric types. Errors will occur if these columns contain non-numeric values.

## LD Reference Panel

LD clumping requires a reference panel in PLINK binary format (`.bed`, `.bim`, `.fam`).
*   Provide the **full path prefix** (without the `.bed`, `.bim`, `.fam` extension) to the reference panel using the `--ld_ref` argument. Example: `--ld_ref /path/to/ld_reference/1kg_eur`.
*   The reference panel should ideally match the ancestry of your GWAS data.

## Instrument Independence (important for MR)

MR methods (IVW, weighted median, MR-Egger) assume the instruments are **statistically independent**. Correlated instruments understate the IVW standard error and inflate false positives. The MR convention is therefore strict LD clumping at **r² < 0.001** (this pipeline's default `--clump_r2`).

> **Do not confuse gene-mapping "independent signals" with MR instruments.** Signal-selection / fine-mapping pipelines (e.g. GCTA-COJO stepwise selection) deliberately keep *conditionally* independent signals at a much looser threshold (often r² < 0.05), and report **marginal** effect sizes. Two such signals can sit a few kb apart at the same locus. Feeding those directly into MR as independent instruments is **not valid**.

**Recommended handling of a pre-selected signal list:** pass it as the exposure and let the pipeline LD-clump it (i.e. **do not** use `--skip_clump`). `ieugwasr::ld_clump` only needs the SNP list plus an LD reference panel, so it will re-prune your signals to r² < 0.001 correctly. Use `--skip_clump` **only** for instruments already independent at MR standard (e.g. cis-QTL instruments pre-clumped at r² < 0.001).

**MHC region:** by default the pipeline **flags** (does not drop) instruments in the MHC (`--mhc_region`, default `6:25000000-34000000`, GRCh37 extended MHC) via an `mhc` column in the instrument table, so you can run a sensitivity analysis with/without them. Set `--mhc_region` to match your GWAS build, or pass `--exclude_mhc` to drop them outright.

**Non-rsID instruments:** IDs that are not rsIDs (e.g. indels named `13:60994514_GT_G`) usually fail rsID-based LD matching and harmonisation and are dropped; the pipeline warns you when it sees them.

### Key independence / QC options

| Option | Default | Purpose |
|---|---|---|
| `--clump_r2` | `0.001` | LD r² threshold for clumping (MR standard). |
| `--clump_kb` | `10000` | Clumping window (kb). |
| `--clump_p` | `5e-8` | P-value threshold for instrument selection. |
| `--skip_clump` | off | Skip clumping. Only for already-independent (r²<0.001) instruments. |
| `--exp_ivs` | none | File of pre-defined instrument SNP IDs; instruments are restricted to this list (no clumping) while `--exp_gwas` provides full summary stats so proxies can be found. |
| `--exp_ncase` / `--out_ncase` | / `ncases` | Case-count **column**; for a **binary** trait (with `--steiger_logodds`) this makes Steiger use the log-odds r². |
| `--exp_ncase_total` / `--out_ncase_total` | none | Total number of cases as a **constant** (when the file has no per-SNP ncase column), analogous to `--*_n_total`. |
| `--exp_prevalence` / `--out_prevalence` | none | Population prevalence for a binary trait's Steiger r² (TwoSampleMR assumes 0.1 if unset). |
| `--mhc_region` | `6:25000000-34000000` | MHC region (CHR:START-END) to flag; set to your build. |
| `--exclude_mhc` | off | Drop MHC instruments before harmonisation. |
| `--analysis_exclude_mhc` | off | Drop flagged MHC SNPs after harmonisation; use for the no-MHC branch with the same inputs to reuse caches. |
| `--proxies` | off | Substitute an LD proxy for any instrument missing from the outcome (uses `--ld_ref`). |
| `--proxy_r2` | `0.8` | Minimum r² for an LD proxy. |
| `--proxy_kb` | `1000` | Window (kb) to search for LD proxies. |
| `--proxy_mem` | `30000` | Memory (MB) cap for the PLINK proxy search. |
| `--steiger` / `--no_steiger` | **on** | Steiger directionality filtering; `--no_steiger` disables it. |
| `--presso` / `--no_presso` | **on** | MR-PRESSO outlier test, run only when the main IVW p-value is <0.05; `--no_presso` disables it. |
| `--presso_nbdist` | `1000` | MR-PRESSO simulations; raise (e.g. 10000) for large instrument sets. |
| `--out_type` | `binary` | `binary` reports odds ratios; `continuous` reports the beta (no OR). |
| `--f_stat` | `10` | Minimum per-SNP F-statistic. |
| `--lib_path` | none | Optional custom R library path (no need to edit the script). |

> **Steiger runs by default. MR-PRESSO is enabled by default but only runs for pairs with IVW p<0.05.** Turn them off with `--no_steiger` / `--no_presso`; MR-PRESSO cost scales ~n_snps² × NbDistribution.

**LD proxies:** with `--proxies`, any instrument that is absent from the outcome GWAS is looked up against the LD reference (`--ld_ref`, via PLINK `--r2`); the highest-r² proxy that is present in **both** the exposure and the outcome (r² ≥ `--proxy_r2`, within `--proxy_kb`) replaces it, using the proxy SNP's own effects in both datasets. A proxy must also meet the exposure association threshold (`--clump_p`, 5e-8 by default) and the F-statistic threshold (`--f_stat`, 10 by default). Substitutions are written to `<prefix>_proxies.tsv` (`instrument`, `proxy`, `r2`). Requires `--ld_ref`; needs the full exposure summary stats (so the proxy has an exposure effect), so it's most useful when the exposure is a full GWAS rather than a pre-selected signal list.

## Performance for GWAS-by-GWAS screens

The pipeline caches cleaned, SNP-subsetted GWAS tables, LD-proxy substitutions, and harmonised data in `./mr_cache` by default. Cache keys include the input file path, size, modification time, column mappings, selected SNPs, and relevant thresholds, so a changed input automatically invalidates its old cache. Put `--cache_dir` on fast shared cluster storage when using multiple jobs; use `--no_cache` for a one-off clean run. Each pair also writes `<prefix>_runtime.tsv`, giving cumulative stage timings for profiling.

For an MHC sensitivity analysis, run the MHC-included configuration first, then repeat the same command and replace the MHC setting with `--analysis_exclude_mhc` (and change `--out_prefix`). This excludes MHC after harmonisation, allowing the cleaned GWAS and harmonised cache entries to be reused. The supplied local and Slurm grid launchers use this mode for their no-MHC configurations.

Independent pairs/configurations can be parallelised as a Slurm job array with [run_mr_slurm_array.sh](run_mr_slurm_array.sh). Create a trusted manifest containing one fully quoted `Rscript pipeline_scripts/mr_pipeline.R ... --cache_dir /shared/mr_cache` command per line, then submit `sbatch --array=1-N%K run_mr_slurm_array.sh manifest.txt`, choosing `K` for the available memory and reference-panel I/O capacity.

## Output Files

The script will generate files in the directory specified by the `--out_prefix` argument (e.g., `output/exposure_vs_outcome_*`):

*   `_harmonized_data.rds`: An R Data Serializable file containing the harmonized data frame used for the MR analysis.
*   `_full_mr_results.csv`: A comma-separated file with the results from all MR methods run, including sensitivity analyses statistics (heterogeneity Q-stat, MR-Egger intercept, MR-PRESSO results if run). Effect columns depend on `--out_type`: for a **binary** outcome the file includes odds ratios (`or`, `or_lci95`, `or_uci95` = `exp(beta)`); for a **continuous** outcome the OR columns are omitted and the effect is the beta with its CI (`b`, `lo_ci`, `up_ci`) — since `exp(beta)` is not an odds ratio for a continuous outcome.
*   `all_processed_mr_results.csv`: A summary file of all processed results and flags based on sensitivity tests.
*   `_exposure_ivs.tsv`: A table of the genetic variants selected as instruments for the exposure after clumping and F-statistic filtering. Includes chr/pos, alleles, EAF, beta/se/p, F-statistic, and an `mhc` flag column.

## Usage Example Single Exposure, Single Outcome

```bash
Rscript pipeline_scripts/mr_pipeline.R \
 --exp_gwas /path/to/exposure_gwas.tsv.gz \
 --out_gwas /path/to/outcome_gwas.txt.gz \
 --ld_ref /path/to/ld_reference/1kg_eur \
 --out_prefix output/exposure_vs_outcome \
 \
 --exp_name "ExposureTrait" \
 --exp_snp "variant_id" \
 --exp_beta "beta_exposure" \
 --exp_se "se_exposure" \
 --exp_ea "effect_allele" \
 --exp_nea "other_allele" \
 --exp_p "p_value" \
 --exp_eaf "eaf" \
 --exp_n "sample_size" \
 \
 --out_name "OutcomeTrait" \
 --out_snp "RSID" \
 --out_beta "BETA" \
 --out_se "SE" \
 --out_ea "ALLELE1" \
 --out_nea "ALLELE0" \
 --out_p "PVAL" \
 --out_eaf "EAF" \
 --out_n "N" \
 --out_ncase "N_cases" \
 --out_type "binary" \
 \
 --clump_p 5e-8 \
 --clump_r2 0.001 \
 --clump_kb 10000 \
 --f_stat 10 \
 --steiger \
 --presso
```

*Adjust file paths, column names (`--exp_*`, `--out_*`), and parameters as needed.*

## Usage Example Multiple Exposures, Multiple Outcomes

Place all your outcome GWAS files (e.g., `.tsv` files) in a single directory (e.g., `outcomes/`).  Place all your exposure GWAS files (e.g., `.tsv` files) in a single directory (e.g., `exposures/`).  
Run:

```bash
Rscript pipeline_scripts/mr_pipeline.R \
  --exposure_dir /path/to/exposures/ \
  --outcome_dir /path/to/outcomes/ \
  --ld_ref /path/to/ld_reference/1kg_eur \
  --out_prefix output/ \
  [other arguments as needed]
```
- The script will process each outcome file in the directory.
- Output files will be named like `output/exposure1_vs_outcome1_full_mr_results.tsv`, `output/exposure1_vs_outcome2_full_mr_results.tsv`, etc.

### Important Note

When running the pipeline with multiple exposures or multiple outcomes (i.e., using `--outcome_dir` or looping over exposures), **all exposure files must have the same column headers, and all outcome files must have the same column headers**.

This ensures that the column mapping arguments (`--exp_*` and `--out_*`) work consistently across all files in the batch.

If your files have different header names, please standardize them or adjust the column mapping arguments accordingly before running the pipeline.

## Trait info files (manifests) — `mr_grid.R`

`mr_pipeline.R` runs **one exposure against one outcome**, with every column name,
sample size and build passed as `--exp_*` / `--out_*` flags — so a batch requires all
files to share headers (see the note above). `mr_grid.R` removes that constraint: you
describe each trait **once** in an **info file** (a CSV manifest), and the driver runs the
whole exposure × outcome grid — both directions and, optionally, the 4-config sensitivity
grid — reading every column/N/type/build/instrument detail per trait from the manifest.
Traits can have completely different headers.

### Manifest columns

One row per trait. Header names are case-insensitive and `.`/space/`_` are interchangeable
(so `chr.col`, `Chr Col`, `chr_col` all match). Required columns are marked ✔.

| Column | Meaning |
|---|---|
| `Short name` ✔ | Trait ID — used in result labels and output filenames, and as the key for `--exposures`/`--outcomes`. |
| `file.name` ✔ | GWAS summary-stats file. Absolute path used as-is; otherwise resolved against `--data_dir`. |
| `snp.col`,`ea.col`,`nea.col`,`beta.col`,`se.col`,`pval.col` ✔ | Column names in that GWAS. |
| `chr.col`,`pos.col`,`eaf.col` | Optional columns; leave blank if absent (EAF-less traits are handled). |
| `N` | Total sample size → constant per-SNP N (enables Steiger). |
| `Ncases` | Number of cases. **Populated ⇒ the trait is treated as BINARY** (odds ratios, case-aware Steiger); **blank ⇒ CONTINUOUS** (beta). |
| `Population prevalence` | Used by Steiger's log-odds r² for binary traits (with `--steiger_logodds`). |
| `Genome build` | `b37`/`b38` — selects the MHC region used to flag/exclude that trait's instruments when it is the exposure. |
| `ivs.file` | Path to the trait's **pre-defined instrument list** (independent signals from the original publication). Drives the *lenient* config (see below). Blank ⇒ the lenient config clumps at r²=0.2 instead. |
| `Downloaded` | `no` (or a blank `file.name`) ⇒ the row is **skipped** with a warning. |

Numbers may contain thousands separators (`"408,112"` is fine).

### How the grid runs

- **Selection.** `--exp_info a.csv[,b.csv]` defines the exposure pool, `--out_info c.csv`
  the outcome pool (the same file may serve both). By default **all** traits are used;
  `--exposures name1,name2` / `--outcomes name1,name2` restrict to a subset of Short names.
- **Directions.** `--bidirectional` also runs each outcome-as-exposure.
- **Sensitivity grid.** `--sensitivity_grid` runs 4 configs per pair, distinguished by
  instrument selection and MHC, with output-filename prefixes matching the rest of the
  pipeline (so `plot_forest.R` works unchanged):
  - `` (lenient) — use the exposure's `ivs.file` (pre-defined signals), keep MHC.
  - `rigid_` (strict) — re-clump the full sumstats at r²=0.001, keep MHC.
  - `noMHC_` — lenient + exclude MHC at analysis.
  - `noMHC_rigid_` — strict + exclude MHC.

  Without `--sensitivity_grid` a single config runs, using `--clump_r2`/`--exclude_mhc`
  (the exposure's `ivs.file` is still used as instruments when present).
- Output files are `<out_prefix><config>_<exposure>_vs_<outcome>_full_mr_results.csv`,
  plus a processed summary `<out_prefix><config>all_processed_mr_results.csv` per config.

### Example

```bash
Rscript pipeline_scripts/mr_grid.R \
  --exp_info read_me_autoimmune_traits.csv \
  --out_info read_me_menopause_traits.csv \
  --data_dir /path/to/gwas \
  --ld_ref   /path/to/ld_reference/Phase3_eur \
  --bidirectional --sensitivity_grid \
  --out_prefix output/autoimmune/ \
  --proxies --proxy_r2 0.8 --proxy_kb 1000
# add --dry_run first to print the resolved traits and the exact grid, then exit
# add --exposures asthma,EO / --outcomes POI,ANM to restrict the sets
```

See `run_mr_grid_slurm.sh` for a ready-to-edit SLURM wrapper.

## Multivariable MR — `mr_mvmr.R`

`mr_mvmr.R` estimates the **joint (conditional) effect** of two or more exposures on
an outcome, using the same manifests. Each invocation is ONE model: every trait in
`--exposures` goes in together, run against each `--outcomes` trait.

How it works:
1. **Per-exposure instruments** — each exposure's `ivs.file` (authoritative) or, if none,
   a clump of its GWAS.
2. **Pool + re-clump** — union the exposures' instruments, tag each SNP with its minimum
   p across exposures, and **re-clump the union** (`--clump_r2`, `--clump_kb`) so the
   combined instrument set is mutually independent (MVMR-IVW assumes this).
3. **Harmonise** all exposures + the outcome at that set (`mv_harmonise_data`).
4. **MVMR-IVW** (`TwoSampleMR::mv_multiple`) for the conditional effects, plus **conditional
   F-statistics** (instrument strength per exposure), **Q / Q-pval** (horizontal pleiotropy),
   and **qhet** (weak-instrument-robust estimate) from the `MVMR` package; optional
   **MVMR-Egger** intercept (`--egger`, needs `MendelianRandomization`).
5. Also emits the **univariate IVW** for each exposure→outcome and a
   **univariable-vs-multivariable forest**.

```bash
# asthma + smoking, modelled jointly, on POI and ANM:
Rscript pipeline_scripts/mr_mvmr.R \
  --exp_info manifests/read_me_respiratory_traits.xlsx \
  --out_info manifests/read_me_reproductive_traits.xlsx \
  --exposures asthma,smoking --outcomes poi,anm \
  --ld_ref /path/to/Phase3_eur --out_prefix output/mvmr/ \
  --clump_r2 0.01 --clump_kb 1000 --egger
# exposures spanning two manifests: pass both to --exp_info (comma-separated)
#   --exp_info respiratory.xlsx,reproductive.xlsx --exposures asthma,bweight
# --exclude_mhc drops MHC instruments; --dry_run prints the model and exits
```

Outputs per model×outcome: `<prefix><expA_expB>_vs_<outcome>_mvmr.csv` (one row per
exposure: conditional `b/se/pval`, `or`/CI for a binary outcome, `cond_Fstat`, `Qstat`,
`Qpval`, `b_qhet`, and — with `--egger` — `egger_*`), a combined `_all_mvmr_results.csv`,
and a `_mvmr_forest.png`. See `run_mr_mvmr_slurm.sh`. Interpretation: a **conditional
F ≥ 10** per exposure indicates adequate instrument strength; if an exposure's univariate
effect shrinks toward the null once the other is included, its marginal effect was
(partly) mediated/confounded by the co-exposure.

### Limitation

Each manifest row has a single `pval.col`, so a trait that needs a *different* p-value
column depending on whether it is the exposure or the outcome (e.g. the ANM
`P_BOLT_LMM` vs `P_BOLT_LMM_INF` case) can't be expressed in one row — keep those on the
per-trait bash runner (`run_mr_slurm_eos.sh`). Standard full-summary-statistic traits
(one p-value column) are the target here.

## Special Case: QTLs as Exposures

If you are using QTLs (e.g., eQTLs, pQTLs, mQTLs) as exposures, you can provide a directory of QTL summary statistics files (e.g., qtl_exposures/) with one file per molecular trait (for example, one file per gene, protein, or metabolite).
Each file should contain the summary statistics for the independent QTLs (IVs) for that molecular trait.

Use the --exposure_dir option to specify this directory.
Use the --skip_clump option **only if** the IVs are already independent at MR standard (r² < 0.001). If your QTL instruments were selected at a looser threshold, omit --skip_clump so the pipeline LD-clumps them (see [Instrument Independence](#instrument-independence-important-for-mr)).

Run :

```bash
Rscript pipeline_scripts/mr_pipeline.R \
  --exposure_dir qtl_exposures/ \
  --outcome_dir outcomes/ \
  --out_prefix output/ \
  --skip_clump \
  [other arguments as needed]
```
- All QTL exposure files must have the same column headers.
- All outcome files must have the same column headers.
- The script will run MR for every QTL trait (file) against every outcome file.
