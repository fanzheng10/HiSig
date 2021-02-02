#!/usr/bin/env bash

R -f ../R/glmnet.R --args sim_conn.txt genescore.tsv sim_ms_impact

python ../parse.py --ont_conn sim_conn.txt --rout sim_ms_impact.impact-w-rand.tsv --terms terms.txt --genes genes.txt --signal genescore.tsv --out sim_ms_impact_summary.tsv

grep -Ff selected_terms.txt sim_ms_impact_summary.tsv |awk '{print $3,$6}' > selected_terms_p.txt