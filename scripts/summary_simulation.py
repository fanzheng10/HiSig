import argparse, os
import pandas as pd
import numpy as np
from sklearn import metrics
# import re

par = argparse.ArgumentParser()
par.add_argument('--out', required=True, help='output file')
args = par.parse_args()

dirs = sorted([d for d in os.listdir('.') if os.path.isdir(d)])

odir = os.getcwd()

records = []

for d in dirs:
    os.chdir(odir)
    os.chdir(d)

    # re.findall(r"[-+]?\d*\.?\d+|\d+", x)

    N, m, p = [''.join(x[1:]) for x in d.split('_')[:3]]
    k = d.split('_')[-1]

    if not os.path.isfile('selected_terms.txt'):
        continue

    terms = [l.strip() for l in open('terms.txt')]
    sel_terms = [l.strip() for l in open('selected_terms.txt')]
    df_all = pd.DataFrame(index = terms)
    df_all['selected'] = 0
    df_all.loc[sel_terms, 'selected'] = 1

    # Hisig
    hisig_f = 'sim_ms_impact_summary.tsv'
    if os.path.isfile(hisig_f):
        df_hisig = pd.read_csv(hisig_f, sep='\t', index_col=0)
        df_hisig = df_hisig[['System_name', 'p']]
        df_hisig.set_index('System_name', inplace=True, drop=True)
        df_hisig.rename(columns={'p':'p_hisig'}, inplace=True)
        df_all = df_all.merge(df_hisig, how='left', left_index=True, right_index=True)
        df_all['p_hisig'].fillna(1.0, inplace=True)
        fpr, tpr, _ = metrics.roc_curve(np.array(df_all['selected']),
                                          -np.array(df_all['p_hisig']))
        auc = metrics.auc(fpr, tpr)
        records.append([N, m, p, auc, 'hisig', k])

    # GSEA
    gsea_f = 'gsea_results.txt'
    if os.path.isfile(gsea_f):
        df_gsea = pd.read_csv(gsea_f, sep='\t', index_col=0)
        df_gsea.set_index('V1', inplace=True, drop=True)
        df_gsea.rename(columns={'p': 'p_gsea'}, inplace=True)
        df_all = df_all.merge(df_gsea, how='left', left_index=True, right_index=True)
        fpr, tpr, _ = metrics.roc_curve(np.array(df_all['selected']),
                                        -np.array(df_all['p_gsea']))
        auc = metrics.auc(fpr, tpr)
        records.append([N, m, p, auc, 'gsea', k])

    # normal Lasso
    lasso_f = 'normal_lasso.txt'
    if os.path.isfile(lasso_f):
        df_all['lasso_selected'] = 0
        lasso_sel = [l.strip() for l in open(lasso_f).readlines()]
        df_all.loc[lasso_sel, 'lasso_selected'] = 1
        fpr, tpr, _ = metrics.roc_curve(np.array(df_all['selected']),
                                       -np.array(df_all['lasso_selected']))
        auc = metrics.auc(fpr, tpr)
        records.append([N, m, p, auc, 'lasso', k])


os.chdir(odir)
df_out = pd.DataFrame.from_records(records, columns=['N', 'm', 'p', 'AUC', 'Method', 'trial'])
df_out = df_out.round(3)
df_out.to_csv(args.out, sep='\t', index=False)