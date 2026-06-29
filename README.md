# CeRberus

CeRberus is an R package for aggregating guide-level CRISPR genetic interaction (GI) scores to gene-pair-level scores using limma-based models.

CeRberus is intended to serve two roles:

- as the R analysis component of the broader Cerberus Python package workflow, and
- as a standalone R package that can be used independently in other GI scoring pipelines.

## What CeRberus does

Starting from guide-level GI scores in long-table format, CeRberus:

1. imports the score table,
2. infers the screen structure,
3. generates alternative replicate/blocking configurations,
4. estimates duplicate correlations,
5. fits limma-based models,
6. aggregates guide-level scores to gene-pair-level GI scores, and
7. optionally writes result tables and intermediate objects to disk.

CeRberus supports both:

- **fixed-pair screens**, where query and library gene sets differ, and
- **multiplex screens**, where query and library gene sets are identical.

## Installation

CeRberus depends on CRAN packages and on [`limma`](https://bioconductor.org/packages/release/bioc/html/limma.html) from Bioconductor.

### 1. Install BiocManager if needed

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
```

### 2. Install limma

```r
BiocManager::install("limma")
```

### 3. Install CeRberus from GitHub

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("thorohde/CeRberus")
```

### 4. Load the package

```r
library(CeRberus)
```

## Two ways to use CeRberus

CeRberus can be used in two main ways.

### 1. Full pipeline execution with `full_run()`

Use `full_run()` when you want a YAML-driven analysis workflow that reads input data, evaluates configurations, computes GI scores, and optionally writes output files.

```r
full_run("input_yaml_file.yaml")
```

### 2. Direct object construction with `GIScores()`

Use `GIScores()` when you want finer control inside an R workflow.

```r
gi_obj <- GIScores(input = my_scores)
```

This creates a CeRberus screen object that can then be processed with downstream package methods.

## Quick start: YAML-based workflow

The typical standalone workflow is:

1. store guide-level GI scores in a `.csv` or `.rds` file,
2. create a YAML instruction file,
3. run `full_run()`, and
4. inspect the returned objects and/or written output files.

### Minimal example

```r
library(CeRberus)

result <- full_run("input_yaml_file.yaml")
```

By default, `full_run()` returns the selected CeRberus screen object configuration as a named list. If `return_output = FALSE`, it returns `NULL` after running the pipeline.

## YAML instruction file

The YAML instruction file controls the pipeline run.

### Required fields

| parameter | description |
|---|---|
| `scores_file` | Path to a guide-level GI score file in `.csv` or `.rds` format |
| `output_directory` | Directory where CeRberus writes outputs |

### Optional fields

| parameter | description | default |
|---|---|---|
| `FDR` | Multiple-testing correction method. Currently supported: `BH`, `bonferroni` | `BH` |
| `overwrite_output` | Whether to write output files to `output_directory` | `TRUE` |
| `make_symmetric` | Aggregate symmetric multiplex screens into one score per unordered gene pair | `FALSE` |
| `keep_all_configurations` | Keep all tested replicate/blocking configurations instead of only the selected one | `FALSE` |
| `verbose` | Print additional progress information and screen summaries | `FALSE` |

### Example YAML file

```yaml
scores_file: "path/to/guide_scores.csv"
output_directory: "path/to/output"
FDR: "BH"
overwrite_output: true
make_symmetric: false
keep_all_configurations: false
verbose: false
```

## Input data format

The input score table must be in **long format**, with one guide-pair observation per row.

CeRberus expects the following columns by default:

- `query_gene`
- `library_gene`
- `bio_rep`
- `tech_rep`
- `guide_pair`
- `GI`

Example input:

| bio_rep | tech_rep | guide_pair | query_gene | library_gene | GI |
|---|---|---|---|---|---|
| b1 | t1 | g1 | IGF2 | ADRB2 | 1.02869961 |
| b1 | t1 | g2 | IGF2 | ADRB2 | -0.46078453 |
| b1 | t1 | g3 | IGF2 | ADRB2 | 1.01662498 |
| b1 | t2 | g1 | IGF2 | ADRB2 | -0.94890821 |
| b1 | t2 | g2 | IGF2 | ADRB2 | -0.44307736 |
| b1 | t2 | g3 | IGF2 | ADRB2 | -0.51468140 |
| b2 | t1 | g1 | MTOR | IMPDH1 | -1.14823650 |
| b2 | t1 | g2 | MTOR | IMPDH1 | -2.68735930 |
| b2 | t1 | g3 | MTOR | IMPDH1 | -0.67141920 |
| b2 | t2 | g1 | MTOR | IMPDH1 | -0.79467931 |
| b2 | t2 | g2 | MTOR | IMPDH1 | -2.05650114 |
| b2 | t2 | g3 | MTOR | IMPDH1 | -0.32426398 |

If needed, column names can be customized when using `GIScores()` directly.

## Analysis procedure

CeRberus follows this general procedure:

1. **Infer screen type**  
   Based on the query and library gene sets, CeRberus determines whether the data represent a fixed-pair or multiplex screen.

2. **Generate alternative replicate/blocking configurations**  
   CeRberus tests different ways of using or collapsing `tech_rep`, `bio_rep`, and `guide_pair` layers.

3. **Estimate duplicate correlations**  
   CeRberus evaluates which configuration produces the most suitable duplicate-correlation structure for downstream modelling.

4. **Fit limma models and compute gene-level GI scores**  
   Guide-level scores are aggregated to gene-pair-level scores, p-values are computed, and multiple-testing correction is applied.

5. **Return and optionally export results**  
   Depending on settings, CeRberus returns the selected configuration and may also write output files to disk.

## Output

When `overwrite_output = TRUE`, CeRberus writes output files to `output_directory`. These can include:

- `all_GI_objects.rds` — intermediate CeRberus screen objects before final selection
- `duplicateCorrelationPlot.png` — duplicate-correlation summary plot
- `duplicate_correlation.csv` — duplicate-correlation summary table
- `GI_scores_<configuration>.csv` — gene-level GI scores for retained configurations

When `overwrite_output = FALSE`, CeRberus still runs the analysis but does not write these files.

## Standalone use in R pipelines

Although CeRberus is part of the larger Cerberus Python package ecosystem, it is fully usable as a standalone R package.

This is useful when you:

- already have guide-level GI scores available in R,
- want to integrate CeRberus into an R-based analysis workflow, or
- want more direct control over object construction and downstream processing.

## Development status

CeRberus is under active development. If you use it as a standalone package, it is a good idea to pin a specific GitHub commit or release in reproducible analysis workflows.
