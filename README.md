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
    *   **Sample Size (N):** Required for calculating instrument strength (F-statistic) and for Steiger filtering (`--steiger`). Map using `--exp_n` / `--out_n`.
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
| `--mhc_region` | `6:25000000-34000000` | MHC region (CHR:START-END) to flag; set to your build. |
| `--exclude_mhc` | off | Drop MHC instruments instead of flagging them. |
| `--lib_path` | none | Optional custom R library path (no need to edit the script). |

## Output Files

The script will generate files in the directory specified by the `--out_prefix` argument (e.g., `output/exposure_vs_outcome_*`):

*   `_harmonized_data.rds`: An R Data Serializable file containing the harmonized data frame used for the MR analysis.
*   `_full_mr_results.csv`: A comma-separated file with the results from all MR methods run, including sensitivity analyses statistics (heterogeneity Q-stat, MR-Egger intercept, MR-PRESSO results if run).
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
