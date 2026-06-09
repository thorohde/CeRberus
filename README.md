# CeRberus

The CeRberus module is part of the Cerberus package and performs the final steps of aggregating guide-level genetic interaction (GI) scores to gene-pair level scores. It can also be used independently as part of other `R` GI scoring pilelines. 

# Quick start

1 Install `devtools` from CRAN:

```{r}
if (!require("devtools", quietly = TRUE)) {
	install.packages("devtools")
	}
```

2. Install the `CeRberus` package:

```{r}
devtools::install_github("thorohde/CeRberus")
```


3. Follow the [Bioconductor](https://www.bioconductor.org/) installation routine to install the [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) package:

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")
```

4. Load the `CeRberus` package into the workspace:

```{r}
library(CeRberus)
```

5. After storing guide GI scores in a csv file, call the modules `full_run()` routine with an instruction yaml file: 

```{r}
full_run(yaml_fpath = "input_yaml_file.yaml")
```



# Input parameters

| parameter | description | default | required |
|---|---|---|---|
| scores_file | path to a `csv` file of guide GI scores |  | * |
| output_directory | path where the package stores the computed scores |  | * |
| FDR | method used for to estimate the False Discovery rate. | `BH` |  |
| skip_update | Can be used to suppress package updates. | `False` |  |
| make_symmetric | Can be used to aggregate symmetric multiplex screens to produce one GI and FDR per gene pair. Useful if all pairs exist twice, e.g. Cas12/Cas9 & Cas9/ Cas12. | `False` |  |
| verbose | should the pipeline print a feedback while running? | `True` |  | 
| keep_all_configurations | Skip discarding configurations with negative duplicateCorrelation. |  |  |

# Input file format

The input table of guide GI scores is required to be in long format, one observation per row. The columns `bio_rep`, `tech_rep` and `guide_pair` describe the experimental setup and will be used as replication layers by `CeRberus`. 


| bio_rep | tech_rep | guide_pair | query_gene | library_gene | GI |
|---|---|---|---|---|---|
| b1 | t1 | g1 | IGF2 | ADRB2 | 1.02869961 |
| b1 | t1 | g2 | IGF2 | ADRB2 | -0.46078453 |
| b1 | t1 | g3 | IGF2 | ADRB2 | 1.01662498 |
| b1 | t2 | g1 | IGF2 | ADRB2 | -0.94890821 |
| b1 | t2 | g2 | IGF2 | ADRB2 | -0.44307736 |
| b1 | t2 | g3 | IGF2 | ADRB2 | -0.5146814 |
| b2 | t1 | g1 | MTOR | IMPDH1 | -1.1482365 |
| b2 | t1 | g2 | MTOR | IMPDH1 | -2.6873593 |
| b2 | t1 | g3 | MTOR | IMPDH1 | -0.6714192 |
| b2 | t2 | g1 | MTOR | IMPDH1 | -0.79467931 |
| b2 | t2 | g2 | MTOR | IMPDH1 | -2.05650114 |
| b2 | t2 | g3 | MTOR | IMPDH1 | -0.32426398 |



# Procedure

1. After importing the guide GI scores, CeRberus infers the experimental structure based on the query and library gene set. 

- Fixed Pair Screen: Query genes and library genes are different.
- Multiplex Screen: Query and library gene sets are identical.


2.  CeRberus tests which configuration produces desirable duplicateCorrelation values. A positive correlation between either the biological or technical replicates is expected. 

3. For the configurations that pass these filters, genetic interactions are computed, using linear model fits. P-values are computed using Empirical Bayes Statistics, and corrected for multiple testing by a method of choice.

4. The computed scores are exported to the given output directory. 


