import os,sys
import pandas as pd
import numpy as np


### input ###
if len(sys.argv) != 5:
    print('usage:1-directory for hisig results; 2-prefix of the ontology file, including the path; 3-output directory; 4-q-values to use')
    exit(0)

dir_hisig_results = sys.argv[1]
ont_prefix = sys.argv[2]
mut_freq_per_term = ont_prefix + '.freq'
size_per_term = ont_prefix + '.genes'
outdir = sys.argv[3]
if not os.path.isdir(outdir):
    os.makedirs(outdir)
###

### global variables ###

qcut = float(sys.argv[4])
hisig_summary_tag = '_ms_impact_summary.tsv'
generank_rescale_tag = '_rescaled.tsv'
n_patients = '/cellar/users/f6zheng/Data/human_cancer_hierarchy_submit/cancerdata/analyzed/mutsigcv_pancanatlas/n_patients.txt'
mut_count_observed = '/cellar/users/f6zheng/Data/human_cancer_hierarchy_submit/cancerdata/analyzed/mutsigcv_pancanatlas/mutsigcv1.4_pancanceratlas_observe.csv'
# f_cohort_labels = ''
# cohorts = [l.strip() for l in open(f_cohort_labels).readlines()]
cohorts = ['blca', 'brca', 'coad', 'gbm', 'hnsc', 'kirc', 'lihc', 'luad', 'lusc', 'ov', 'skcm', 'stad', 'ucec']
###

df_mut_freq_per_term = pd.read_csv(mut_freq_per_term, sep='\t', index_col=0)
dict_mut_freq_per_term = {}
for c in cohorts:
    dict_mut_freq_per_term[c] = df_mut_freq_per_term[c + '_mut'].to_dict()
dict_size_per_term = pd.read_csv(size_per_term, sep='\t', index_col=0, header=None)[1].to_dict()


df_mut_count_observed = pd.read_csv(mut_count_observed, sep='\t', index_col=0)
dict_mut_count_observed = {}
for c in cohorts:
    dict_mut_count_observed[c] = df_mut_count_observed[c+'_observed'].to_dict()

dict_npat = pd.read_csv(n_patients, header=None, sep='\t', index_col=1)[0].to_dict()


### parse hisig outputs, combine results from all cohorts
hisig_summary_dfs = {}
generank_rescale_dfs = {}
combined_dfs = {}
nsig = {}
sig_genes = {}

for c in cohorts:
    hisig_summary_dfs[c] = pd.read_csv(dir_hisig_results+'/'+c+hisig_summary_tag, sep='\t', index_col=0)
    hisig_summary_dfs[c]['System_name'] = hisig_summary_dfs[c]['System_name'].astype(str)
    generank_rescale_dfs[c] = pd.read_csv(dir_hisig_results+'/'+c+generank_rescale_tag, sep='\t', index_col=0)
    generank_rescale_dfs[c].index = generank_rescale_dfs[c].index.astype(str)
    hisig_summary_dfs[c] = hisig_summary_dfs[c].loc[hisig_summary_dfs[c]['q'] < qcut, :]
    # combine two data frames for each cohort
    combined_dfs[c] = hisig_summary_dfs[c].merge(generank_rescale_dfs[c], how ='left',
                                                left_on = 'System_name', right_index=True)
    sig_genes[c] = []
    combined_dfs[c]['inferred_system_mutation_frac'] = combined_dfs[c]['inferred_system_mutation_count']/dict_npat[c]
    combined_dfs[c] = combined_dfs[c].round(3)
    nsig[c] = combined_dfs[c].shape[0]

    cols_keep = ['q', 'Genes', 'inferred_system_mutation_count', 'inferred_system_mutation_frac',
                 'adjusted_ranked_gene_names']
    combined_dfs[c] = combined_dfs[c][['System_name'] + cols_keep]
    combined_dfs[c].rename(columns={cl: c.upper() + '-' + cl for cl in cols_keep}, inplace=True)

    # re-ordered the mutation count according to "adjusted ranked gene names"
    for i, row in combined_dfs[c].iterrows():
        if isinstance(row[c.upper() + '-adjusted_ranked_gene_names'], str):
            genes_reordered = row[c.upper() + '-adjusted_ranked_gene_names'].split('|')
            counts_reordered = '|'.join(['{:d}'.format(int(dict_mut_count_observed[c][g])) for g in genes_reordered])
            combined_dfs[c].loc[i, c.upper() + '-adjusted_ranked_mutation_counts'] = counts_reordered
            combined_dfs[c].loc[i, c.upper() + '-adjusted_top2_genes'] = genes_reordered[0] + ' ' + genes_reordered[1]
            sig_genes[c].extend(genes_reordered[:2])
        combined_dfs[c].loc[i, c.upper() + '-raw_mut_freq'] = dict_mut_freq_per_term[c][row['System_name']]

df_combined = combined_dfs[cohorts[0]].copy()
for i in range(1, len(cohorts)):
    df = combined_dfs[cohorts[i]].copy()
    df_combined = df_combined.merge(df, how='outer', left_on='System_name', right_on='System_name')
###

### summary statistics ###

## number of sig results of each study
df_n_sig = pd.Series(nsig)
df_n_sig.to_csv(outdir + '/num_sig_results_each_study_q{}.tsv'.format(qcut), sep='\t')

## number of cancer for each system
cols_q = [c for c in df_combined.columns if c.endswith('-q')]
df_combined['N cancer'] = np.sum(df_combined[cols_q].notnull(), axis=1)

## activate cancer types
for i, row in df_combined.iterrows():
    active = df_combined.loc[i, cols_q].notnull()
    active = ','.join([c.split('-')[0] for c in active.loc[active==1].index.tolist()])
    df_combined.loc[i, 'Activated in cancer'] = active

## module size
df_combined['Size'] = np.array([dict_size_per_term[t] for t in df_combined['System_name'].tolist()])

###

## gene nominations ###
for c in cohorts:
    sig_genes[c] = sorted(list(set(sig_genes[c])))
    with open(c + '.sig_genes_q{}.txt'.format(qcut), 'w') as fh:
        for g in sig_genes[c]:
            fh.write(g + '\n')
###


### output ###
df_combined.set_index('System_name', inplace=True, drop=True)
highlight_col = ['N cancer', 'Size', 'Activated in cancer']
df_combined = df_combined[highlight_col + [c for c in df_combined.columns if not c in highlight_col]]

df_combined = df_combined.round(3)
df_combined.sort_values(['N cancer', 'Size'], ascending=[False, True], inplace=True)

df_combined.to_csv(outdir+'/mutation_allcohort_q{}.tsv'.format(qcut), sep='\t')


