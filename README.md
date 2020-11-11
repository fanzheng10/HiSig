# HiSig
Enrichment analysis for nested and overlapping gene sets

# News

HiSig now available as a self-contained R package!  

Install beta version via `devtools::install_github("fanzheng10/HiSig", subdir="HiSig")`. See basic usage [here](https://github.com/fanzheng10/HiSig/blob/master/HiSig/examples/main.R). More documents and vignettes are under construction. 

# Introduction
HiSig is a program to identify a parsimonious list of gene sets (systems) to explain the signals associated with individual genes. It satisfies the purpose of Gene Set Enrichment Analysis, but much reduces the redundancy in results when provided gene sets that are overlapping or nested (e.g. the gene sets in the Gene Ontology, or other hierarchically structured gene sets). HiSig applies Lasso at multiple stringency of L1-regularization.

We use HiSig to study the gene sets which significantly aggregate somatic mutations in cancer. A manuscript is under preparation. 

# Dependencies
* A working installation of R. We have tested on R 3.4.  Require libraries `glmnet`, `Matrix`, and `parallel`
* Python package `statsmodels (0.9)`, `numpy` and `scipy` (any version should work)
* For efficient permutation test, need multiple CPU cores.


# Usage

Sample files are in the `/examples` directory.

## 1. prepare the input
One can start with a file describing a hierarchical model, which is a 3-column text file  defined in the `DDOT` package (see `sample.ont`),  and another 2-column text file with signals on leaf-nodes in the hierarchical model (see `sample_genescore.tsv`). In our use case, leaf-nodes are interpreted as genes and the signals on leaf-nodes are interpreted as the number of observed mutations (which could have been appropriately adjusted or transformed) of each gene in a cancer cohort.  

Example usage:
`python prepare_input.py --ont sample.ont --sig sample_genescore.tsv --out sample_signal.txt`


After running `prepare_input.py`, one should get the following files:   
(1) A sparse binary matrix defining gene-to-system (genes as row; systems as columns) membership (in TXT format, see `sample_conn.txt`). It corresponds to in the `[I, H]` in the following figure. Here, `I` is an identity matrix since we treat individual genes as systems as well, in order to prevent systems where most signals can be explained by a single gene.  
(2) A text file with real values (see `sample_signals.txt`), genes in the input file (`sample_genescore.tsv`) but not in the hierarchy (`sample.ont`) will be omitted. It corresponds to the `y` vector in the following figure.  
(3) Two files called `genes.txt` and `terms.txt`, they will be needed in the step 3.


<p align="center">
  <img src="fig1.png" width="600" align="center">
</p>



## 2. running Lasso regression

`R -f R/glmnet.R --args sample_conn.txt sample_signals.txt sample_ms_impact 10`

The first two arguments of this script are the 2 outputs of the previous step; the 3rd argument (`"sample_ms_impact"`) defines the file name of the R script output; the 4th argument (`"10"`) is for batch size of permutation. The batch number is 10 to enable parallelization. The number of total permutation is `batch_number * batch_size`, and thus it is 100 in the demo. By default, batch size is set to 1000, so it performs 10000 permutations.

**By default the script use 7 CPU cores; to change it, edit the `max_cores` in `glmnet.R` script**  
**P-values are calculated in an empirical distribution. So the minimum possible p-value depends on the number of permutations**

This step generates two outputs: `sample_ms_impact.coef` and `sample_ms_impact.impact-w-rand.tsv` (big files, not included in `/examples`)

## 3. parse the results

Use `parse.py` to parse the results. Example of usage:

`python parse.py --ont_conn sample_conn.txt --rout sample_ms_impact.impact-w-rand.tsv --terms terms.txt --genes genes.txt --signal sample_signals.txt --out sample_ms_impact_summary.tsv`

The final result is `sample_ms_impact_summary.tsv`, in which each row is a gene set (system); gene sets are ordered by their q-value (Benjamini-Hochberg FDR). The columns `Mutation model input` and `Rank of model` represent genes' signals in the input and their ranks among all genes, to help understand the results.   
