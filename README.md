
Given a hierarhical model of inter-connected systems, HiSig is a program that searches for a parsimonious set of systems as features explaining the signals observed on the leaves (bottom nodes) of such hierarchy, at multiple resolutions. 

# Dependencies
* The DDOT (`Data-driven Ontology Toolkit`) package (url), ensure all the Python dependencies specified there.   (**TODO: make it unnecessary to install DDOT**)
* A working installation of R. We have tested on R 3.4.  Require libraries `glmnet`, `Matrix`, and `parallel`
* Python packages `statsmodels`


# Usage

The Jupyter notebook `examples/demo.ipynb` illustrates the usage of the package step by step. 

## prepare the input
One should start with a file describing a hierarchical model, which is a 3-column text file  defined in `DDOT` (see `examples/sample.ont`),  and another 2-column text file with signals on leaves nodes (see `examples/`). In our use case, leaves nodes are interpreted as genes and the signals on leaves are interpreted as the (transformed) number of observed mutations of each gene. 

After running `prepare_input.py`, one should get two files: (1) a sparse matrix defining gene-to-system membership (in text format, see `examples/sample_conn.txt`); (2) a text file with real values (see `examples/gene_signals.txt`), genes in `blca_signal.txt` that are not in the input hierarchy (`examples/sample.ont`) will be omitted. 

## running Lasso regression

`R -f R/glmnet.R --args examples/sample_conn.txt example/gene_signals.txt example/sample_ms_impact`

## parse the results

