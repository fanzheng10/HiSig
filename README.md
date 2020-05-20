# HiSig

![Figure](figs3.png)


Given a hierarhical model of inter-connected systems, HiSig is a program that searches for a parsimonious set of systems as features explaining the signals observed on the leaves (bottom nodes) of such hierarchy, at multiple resolutions. 

# Dependencies
* The DDOT (`Data-driven Ontology Toolkit`) package (url), ensure all the Python dependencies specified there.
* A working installation of R. We have tested on R 3.4.  Require libraries `glmnet`, `Matrix`, and `parallel`
* Python package `statsmodels`
* for efficient permutation test, need multiple CPU cores.


# Usage

The Jupyter notebook `examples/demo.ipynb` illustrates the usage of the package step by step. 

## prepare the input
One should start with a file describing a hierarchical model, which is a 3-column text file  defined in the `DDOT` package (see `sample.ont`),  and another 2-column text file with signals on leaves nodes (see `sample_genescore.tsv`). In our use case, leaves nodes are interpreted as genes and the signals on leaves are interpreted as the (transformed) number of observed mutations of each gene.   
**TODO: sample_genescore is not here yet**

Example usage:
`python prepare_input.py --ont sample.ont --sig sample_genescore.tsv --out sample_signal.txt`


After running `prepare_input.py`, one should get two files: (1) a sparse matrix defining gene-to-system membership (in text format, see `sample_conn.txt`); (2) a text file with real values (see `sample_signals.txt`), genes in the input file (`sample_genescore.tsv`) but not in the hierarchy (`sample.ont`) will be omitted. 


## running Lasso regression

`R -f R/glmnet.R --args sample_conn.txt sample_signals.txt sample_ms_impact 10`

The first two arguments of this script are the 2 outputs of the previous step; the 3rd argument defines the file name of the R script output; the 4th argument (optional) is for batch size of permutation. The batch number is 10 to enable parallelization. The number of total permutation is `batch_number * batch_size`, and thus it is 100 in the demo. By default, batch size is set to 1000, so it performs 10000 permutations.

**By default the script use 7 CPU cores; to change it, edit the `max_cores` in `glmnet.R` script**

This step generates two outputs: `sample_ms_impact.coef` and `sample_ms_impact.impact-w-rand.tsv`.

*TODO: the second file is too big and thus not included; create small examples later*

## parse the results

Use `parse.py` to parse the results. Example of usage:

`python parse.py --ont sample.ont --rout sample_ms_impact.impact-w-rand.tsv --signal sample_signals.txt --out sample_ms_impact_summary.tsv`

The final result is `sample_ms_impact_summary.tsv`, in which each row is a gene set (system); gene sets are ordered by their q-value (Benjamini-Hochberg FDR). The columns `Mutation model input` and `Rank of model` represent genes' signals in the input and their ranks among all genes, to help understand the results.   