import argparse
import pandas as pd
import numpy as np
from scipy.sparse import *
from statsmodels.stats.multitest import *
from HiSig.utils import *

def parse_r_output(rout, ont_conn, terms, genes, signal,
                   signal2=None, outf=None, min_term_size=2, node_attr=None):
    '''
    
    :param ont: 
    :param rout: 
    :param signal: 
    :param signal2: 
    :param outf: 
    :param min_term_size: 
    :param node_attr: 
    :return: 
    '''

    coo_gene_term = np.loadtxt(ont_conn).astype(int)
    mat_gene2term = coo_matrix((np.ones(coo_gene_term.shape[0],), (coo_gene_term[:, 0], coo_gene_term[:, 1])),
                               )
    mat_gene2term = mat_gene2term.tocsr()

    gene_names = [l.strip() for l in open(genes).readlines()]
    term_names = [l.strip() for l in open(terms).readlines()]
    term_sizes = np.sum(mat_gene2term, axis=0)

    ngenes = len(gene_names)

    assert ngenes == len(signal), 'error in [signal]: incorrect number of genes'
    if isinstance(signal2, np.ndarray):
        assert ngenes == len(signal2), 'error in [signal2]: incorrect number of genes'
    assert isinstance(rout, list), 'error: [rout] needs to be a list'

    rout_batches = []
    for i in range(len(rout)):
        rout_batch = pd.read_csv(rout[i], sep="\t", header=None, index_col=0)
        rout_batches.append(rout_batch)
    real_result = np.array(rout_batches[0][1])

    # calculate pvalues
    rand_weights = []
    for i in range(len(rout)):
        rand_weights.append(np.array(rout_batches[i].iloc[:, 1:]))
    all_rand_weights = np.concatenate(rand_weights, axis=1)
    n_shuffle = all_rand_weights.shape[1]
    real_result_rep = np.repeat(real_result[:, np.newaxis], n_shuffle, axis=1)
    n_reject = 1.0 * np.sum(all_rand_weights >= real_result_rep, axis=1)

    n_reject[n_reject == 0] = 1.0
    pvals = n_reject / n_shuffle

    if min_term_size == 2:
        multitest = multipletests(pvals[ngenes:], alpha=0.3, method='fdr_bh', is_sorted=False)
        qvals = multitest[1]
        df = printModuleProfile(real_result, mat_gene2term, term_names, gene_names, signal, signal2=signal2, pvals=pvals, qvals=qvals,
                                no_gene=True)
    else:
        term_used = np.array([i for i in range(len(term_names)) if term_sizes[i] >= min_term_size])
        multitest = multipletests(pvals[ngenes:][term_used], alpha=0.3, method='fdr_bh', is_sorted=False)
        qvals = np.zeros(len(terms), )
        qvals[term_used] = multitest[1]
        term_names = term_used[term_used]
        mat_gene2term = mat_gene2term[:, term_used]
        df = printModuleProfile(real_result, mat_gene2term, term_names, gene_names, signal, signal2=signal2, pvals=pvals, qvals=qvals,
                                no_gene=True, min_term_size=min_term_size)

    df.sort_values(by=['q', 'p', 'Selection_pressure'], ascending=[True, True, False], inplace=True)
    df.reset_index(inplace=True, drop=True)

    if node_attr != None:
        df2 = pd.read_table(node_attr, sep='\t')
        df = df.merge(df2, how='left', left_on='System_name', right_index=True)

    df = df.round(6)

    if outf !=None:
        df.to_csv(outf, sep='\t')
    return df

if __name__ == "__main__":

    par = argparse.ArgumentParser()
    par.add_argument('--ont_conn', required=True, help='the ontology file')
    par.add_argument('--rout', required=True, nargs = '+', help='the output of R; if there are multiple files concatenate them')
    par.add_argument('--terms', required=True)
    par.add_argument('--genes', required=True)
    par.add_argument('--signal', required=True, help='per gene mutation signal - the actual input of the Lasso regression (transformation of mutation);')
    par.add_argument('--signal2', help='another per gene signal to help interpretation; in this case can be the raw rount of mutations')
    par.add_argument('--out', required=True, help='output of this script')
    par.add_argument('--min_term_size', type=int, default=2, help='term size lower than this will not be considered; by default is 2, i.e. consider all terms')
    par.add_argument('--node_attr', help ='if not None, merge with a data frame with extra information')
    args = par.parse_args()

    signal = np.loadtxt(args.signal)
    signal2 = None
    if args.signal2 != None:
        signal2 = np.loadtxt(args.signal2).astype(int)

    parse_r_output(args.rout, args.ont_conn, args.terms, args.genes, signal,
                   signal2=signal2, outf=args.out, min_term_size=args.min_term_size, node_attr=args.node_attr)