# Mendelian Randomization Pipeline

## Purpose

This R-based pipeline performs a standard two-sample Mendelian Randomization (MR) analysis. It takes GWAS summary statistics for an exposure and an outcome, selects genetic instruments for the exposure using LD clumping, harmonizes the data, runs various MR methods, performs sensitivity analyses, and outputs the results. It is designed to be flexible and efficient, supporting the following scenarios:

- **Single exposure, single outcome:** Analyze one exposure-outcome pair.
- **Multiple exposures, multiple outcomes:** Use a shell loop to analyze all exposures in a directory against all outcomes in another directory.
## Dependencies

### R Packages

You need to have the following R packages installed:

*   `optparse`: For parsing command-line arguments.
*   `dplyr`: For data manipulation.
*   `data.table`: For fast file reading (`fread`) and data manipulation.
*   `TwoSampleMR`: The core package for MR analysis functions.
*   `ieugwasr`: Required for LD clumping using the 1000 Genomes reference panel.
*   `genetics.binaRies`: Helper for `ieugwasr` to find PLINK binaries (optional, if PLINK is not in PATH).
*   `MRPRESSO`: Required if the `--presso` option is used for outlier detection (install from GitHub: `devtools::install_github("rondolab/MR-PRESSO")`).

You can install most of these from CRAN:
`install.packages(c("tidyselect", "generics", "optparse", "dplyr", "data.table", "remotes", "ieugwasr"))`

You can install most of these from GitHub:
`remotes::install_github("MRCIEU/TwoSampleMR")`
`remotes::install_github("MRCIEU/genetics.binaRies")`

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

## Output Files

The script will generate files in the directory specified by the `--out_prefix` argument (e.g., `output/exposure_vs_outcome_*`):

*   `_harmonized_data.rds`: An R Data Serializable file containing the harmonized data frame used for the MR analysis.
*   `_full_mr_results.tsv`: A tab-separated file with the results from all MR methods run, including sensitivity analyses statistics (heterogeneity Q-stat, MR-Egger intercept, MR-PRESSO results if run).
*   `_processed_mr_results.tsv`: (To be implemented) A summary file potentially highlighting key results and flags based on sensitivity tests.
*   `_exposure_ivs.tsv`: (To be implemented) A file listing the genetic variants selected as instruments for the exposure after clumping and F-statistic filtering.

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
Use the --skip_clump option, since the IVs are already independent.

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
