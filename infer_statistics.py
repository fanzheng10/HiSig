import pandas as pd
import argparse
from utils import *
from ddot import *

def infer_stats_for_systems(ont, conn, coef_file, term_name_list, signal, onp, outf=None):
    coef_adjust = redistribute_gene_score(coef_file, conn, signal, exponential=True)

    idy, best_lam = feature_best_lambda(coef_file)
    dict_lam = {idy[i]:best_lam[i] for i in range(len(idy))}

    df_onp = pd.read_table(onp, sep='\t', index_col=0)
    term_names = [l.strip() for l in open(term_name_list).readlines()]
    system_mut_count_pertumor = estimate_nsamples_per_term(ont, coef_adjust, term_names, dict_lam, df_onp, outf=outf, save_per_patient=True)

    return system_mut_count_pertumor


if __name__ == "__main__":
    par = argparse.ArgumentParser()
    par.add_argument('--ont', required=True, help='the ontology file')
    par.add_argument('--conn', required=True, help='sparse matrix, for gene-system membership')
    par.add_argument('--coef', required=True, help='the output of R; if there are multiple files concatenate them')
    par.add_argument('--term_list', required=True, help='a file which contains a list of all system names')
    par.add_argument('--signal', required=True, help='per gene mutation signal - the actual input of the Lasso regression (transformation of mutation);')
    par.add_argument('--onp', required=True, help='tumor profile')
    par.add_argument('--outf', required=True, help='output file')
    args = par.parse_args()

    ont = Ontology.from_table(args.ont, clixo_format=True, is_mapping = lambda x:x[2]=='gene')
    ont.propagate('forward', inplace=True) # important!
    infer_stats_for_systems(ont, args.conn, args.coef, args.term_list, args.signal, args.onp, outf = args.outf)






