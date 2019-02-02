import argparse
import pandas as pd
import numpy as np
from ddot import *
from utils import *

par = argparse.ArgumentParser()
par.add_argument('--ont', required=True, help='the ontology file')
par.add_argument('--rout', required=True, nargs = '+', help='the output of R; if there are multiple files concatenate them')
par.add_argument('--nmut', nargs =2, help='gene-wise mutation profile; first file - the actual input of the Lasso model (transformation of mutation); second file - the raw count of mutations')
par.add_argument('--index1',action='store_true', help='if true, then file starts from V1')
par.add_argument('--out', required=True, help='output of this script')
par.add_argument('--sort_by_pval', action='store_true', help='if true, sort by pvalue')
par.add_argument('--min_term_size', type=int, default=2, help='term size lower than this will not be considered; by default is 2, consider all terms')
par.add_argument('--node_attr', help ='if not None, merge with a data frame with extra information')
par.add_argument('--hiview_info', help ='if not None, merge with a data frame with extra information')
args = par.parse_args()

ont = Ontology.from_table(args.ont, clixo_format=True, is_mapping=lambda x:x[2]=='gene')
ont.propagate('forward', inplace=True)
ngenes = len(ont.genes)

nmut_lasso = np.loadtxt(args.nmut[0])
nmut_original = np.loadtxt(args.nmut[1])
nmut_original = nmut_original.astype(int)

assert ngenes == len(nmut_lasso), 'error of nmut: incorrect number of genes'
assert ngenes == len(nmut_original), 'error of nmut: incorrect number of genes'

index_start = 1
if args.index1:
    index_start=0

rout_batches = []
for i in range(len(args.rout)):
    rout_batch = pd.read_table(args.rout[i], sep="\t", header=None, index_col=0)
    rout_batches.append(rout_batch)
real_result = np.array(rout_batches[0][1][index_start:])

# calculate pvalues
rand_weights = []
for i in range(len(args.rout)):
    rand_weights.append(np.array(rout_batches[i].iloc[index_start:, 1:]))
all_rand_weights = np.concatenate(rand_weights, axis=1)
n_shuffle = all_rand_weights.shape[1]
real_result_rep = np.repeat(real_result[:, np.newaxis], n_shuffle, axis=1)
n_reject = 1.0* np.sum(all_rand_weights >= real_result_rep, axis=1)

n_reject[n_reject == 0] = 1.0
pvals = n_reject / n_shuffle


if args.min_term_size == 2:
    multitest = multipletests(pvals[ngenes:], alpha=0.3, method='fdr_bh', is_sorted=False)
    qvals = multitest[1]
    df = printModuleProfile(real_result, ont, nmut_lasso, nmut_raw=nmut_original, pvals = pvals, qvals = qvals, no_gene=True)
else:
    ont_term_used = np.array([i for i in range(len(ont.terms)) if ont.term_sizes[i] >= args.min_term_size])
    multitest = multipletests(pvals[ngenes:][ont_term_used], alpha=0.3, method='fdr_bh', is_sorted=False)
    qvals = np.zeros(len(ont.terms), )
    qvals[ont_term_used] = multitest[1]
    df = printModuleProfile(real_result, ont, nmut_lasso, nmut_raw=nmut_original, pvals = pvals, qvals = qvals, no_gene=True, min_term_size = args.min_term_size)
# if args.sort_by_pval:
#     outstr2 = printModuleProfile(real_result, ont, nmut, pvals, qvals, sort_by_pval=True)
if args.sort_by_pval:
    df = df.loc[df['q'] != '', :]
    df.sort_values(by=['q', 'p', 'Selection pressure'], ascending=[True, True, False], inplace=True)
    df.reset_index(inplace=True, drop=True)

if args.node_attr !=None:
    df2 = pd.read_table(args.node_attr, sep='\t')
    df = df.merge(df2, how='left', left_on='System clixo/louvain name', right_on='System ID')
if args.hiview_info != None:
    df2 = pd.read_table(args.hiview_info, sep='\t')
    df = df.merge(df2, how='left', left_on='System clixo/louvain name', right_on='Old_name')

df =df.round(6)
df.to_csv(args.out, sep='\t')