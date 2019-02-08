import pandas as pd
import argparse
from utils import *

def infer_stats_for_systems(ont, conn, coef_file, sysinfo, signal, onp, qcut=0.3, outf=None):
    coef_adjust = redistribute_gene_score(coef_file, conn, signal, exponential=True)

    df_terms = pd.read_table(sysinfo, sep='\t', index_col=0)
    df_terms_sub = df_terms.loc[df_terms['q'] <=qcut, :]
    idy, best_lam = feature_best_lambda(coef_file)
    dict_lam = {i:best_lam[i] for i in range(len(idy))}

    df_onp = pd.read_table(onp, sep='\t', index_col=0)
    system_mut_count_pertumor = estimate_nsamples_per_term(ont, coef_adjust, df_terms_sub['System clixo/louvain name'],
                                                           dict_lam, df_onp)
    system_mut_count_inferred = np.sum(system_mut_count_pertumor, axis=1)
    df_terms_sub['inferred_system_mutation_count'] = np.array(system_mut_count_inferred)
    df_terms_sub = df_terms_sub.round(6)
    if outf != None:
        df_terms_sub.to_csv(outf, sep='\t', index=False)
    return system_mut_count_pertumor


if __name__ == "__main__":
    par = argparse.ArgumentParser()
    par.add_argument('--ont', required=True, help='the ontology file')
    par.add_argument('--conn', required=True, help='sparse matrix, for gene-system membership')
    par.add_argument('--coef', required=True, help='the output of R; if there are multiple files concatenate them')
    par.add_argument('--sysinfo', required=True, help='an existing data frame of system information')
    par.add_argument('--signal', required=True, help='per gene mutation signal - the actual input of the Lasso regression (transformation of mutation);')
    par.add_argument('--onp', required=True, help='tumor profile')
    par.add_argument('--outf', required=True, help='output file')
    args = par.parse_args()

    infer_stats_for_systems(args.ont, args.conn, args.coef, args.sysinfo, args.signal, args.onp, args.outf)




