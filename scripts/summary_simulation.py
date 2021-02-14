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

    # Hisig with id
    hisig_f = 'sim_ms_impact_id.txt'
    if os.path.isfile(hisig_f):
        df_hisig = pd.read_csv(hisig_f, sep='\t')
        df_hisig = df_hisig[['System.names','Selection.pressure', 'p']]
        df_hisig.set_index('System.names', inplace=True, drop=True)
        df_hisig.rename(columns={'p':'hisig.p.id','Selection.pressure':'hisig.score.id'}, inplace=True)
        df_all = df_all.merge(df_hisig, how='left', left_index=True, right_index=True)
        df_all['hisig.p.id'].fillna(1.0, inplace=True)
        df_all['hisig.score.id'].fillna(0, inplace=True)
        fpr, tpr, _ = metrics.roc_curve(np.array(df_all['selected']),
                                          -np.array(df_all['hisig.p.id']))
        auc = metrics.auc(fpr, tpr)
        records.append([N, m, p, auc, 'hisig.p.id', k])
        fpr, tpr, _ = metrics.roc_curve(np.array(df_all['selected']),
                                          np.array(df_all['hisig.score.id']))
        auc = metrics.auc(fpr, tpr)
        records.append([N, m, p, auc, 'hisig.score.id', k])

    # HiSig no id
    hisig_f = 'sim_ms_impact_noid.txt'
    if os.path.isfile(hisig_f):
        df_hisig = pd.read_csv(hisig_f, sep='\t')
        df_hisig = df_hisig[['System.names', 'p']]
        df_hisig.set_index('System.names', inplace=True, drop=True)
        df_hisig.rename(columns={'p':'hisig.p.noid','Selection.pressure':'hisig.score.noid'}, inplace=True)
        df_all = df_all.merge(df_hisig, how='left', left_index=True, right_index=True)
        df_all['hisig.p.noid'].fillna(1.0, inplace=True)
        df_all['hisig.score.noid'].fillna(0, inplace=True)
        fpr, tpr, _ = metrics.roc_curve(np.array(df_all['selected']),
                                          -np.array(df_all['hisig.p.noid']))
        auc = metrics.auc(fpr, tpr)
        records.append([N, m, p, auc, 'hisig.p.noid', k])
        fpr, tpr, _ = metrics.roc_curve(np.array(df_all['selected']),
                                          np.array(df_all['hisig.score.noid']))
        auc = metrics.auc(fpr, tpr)
        records.append([N, m, p, auc, 'hisig.score.noid', k])
    
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

    # normal Lasso id
    lasso_f = 'normal_lasso_w_identity.txt'
    if os.path.isfile(lasso_f):
        df_lasso = pd.read_csv(lasso_f, sep='\t')
        df_lasso.set_index('V1', inplace=True, drop=True)
        df_lasso.rename(columns={'coef':'lasso_id'}, inplace=True)
        df_lasso['lasso_id'] = np.abs(df_lasso['lasso_id'])
        df_all = df_all.merge(df_lasso, how='left', left_index=True, right_index=True)
        fpr, tpr, _ = metrics.roc_curve(np.array(df_all['selected']),
                                       np.array(df_all['lasso_id']))
        auc = metrics.auc(fpr, tpr)
        records.append([N, m, p, auc, 'lasso_id', k])

    # normal Lasso no id
    lasso_f = 'normal_lasso_no_identity.txt'
    if os.path.isfile(lasso_f):
        df_lasso = pd.read_csv(lasso_f, sep='\t')
        df_lasso.set_index('V1', inplace=True, drop=True)
        df_lasso.rename(columns={'coef':'lasso_noid'}, inplace=True)
        df_lasso['lasso_noid'] = np.abs(df_lasso['lasso_noid'])
        df_all = df_all.merge(df_lasso, how='left', left_index=True, right_index=True)
        fpr, tpr, _ = metrics.roc_curve(np.array(df_all['selected']),
                                       np.array(df_all['lasso_noid']))
        auc = metrics.auc(fpr, tpr)
        records.append([N, m, p, auc, 'lasso_noid', k])
    
os.chdir(odir)
df_out = pd.DataFrame.from_records(records, columns=['N', 'm', 'p', 'AUC', 'Method', 'trial'])
df_out = df_out.round(3)
df_out.to_csv(args.out, sep='\t', index=False)
