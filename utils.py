import numpy as np
import pandas as pd
from scipy.sparse import *


def printModuleProfile(coef, ont, nmut, nmut_raw = None, pvals = (), qvals = (), top=None, print_limit=50, sort_by_pval = False, min_term_size=2, no_gene=False):
    '''

    :param coef: selection score
    :param ont: ontology object
    :param nmut: input of the lasso model
    :param nmut_raw: number of mutations for each gene
    :param pvals: p values of each system
    :param qvals: q values
    :param top: if not None, only print top N results
    :param print_limit: for terms bigger than this limit do not print the exact gene names
    :param sort_by_pval: if True, sort results by p-value instead of selection score
    :param no_gene: if True, the results only contain terms
    :return:
    '''
    ngenes = len(ont.genes)
    if sort_by_pval:
        idx = np.argsort(pvals)
    else:
        idx = np.argsort(coef)[::-1]
    if top!=None:
        idx = idx[:top]
    outstr = []
    nmut_sort = rankdata(-nmut, method='max')

    # # multiple-test correction
    for i in idx:
        if coef[i] == 0:
            break
        if i >= ngenes: # if the system is a term
            t = ont.terms[i-ngenes]
            if ont.term_sizes[i-ngenes] <  min_term_size:
                continue

            g_in_t_ind = [gi for gi in ont.term_2_gene[t]]
            g_in_t_ind_sort = sorted(g_in_t_ind, key = lambda x:nmut_sort[x])

            g_in_t = [ont.genes[gi] for gi in g_in_t_ind_sort]

            g_nmut = [nmut[gi] for gi in g_in_t_ind_sort]
            g_rank = [nmut_sort[gi] for gi in g_in_t_ind_sort]
            g_nmut_raw = [nmut_raw[gi] for gi in g_in_t_ind_sort]

            if len(g_in_t) > print_limit:
                outlist = [i+1, t, '{} genes including: {}'.format(len(g_in_t), ','.join(g_in_t[:print_limit])), '{:.6f}'.format(coef[i])]
                if (len(pvals) > 0) and (len(qvals) > 0):
                    outlist.extend(['{:.4f}'.format(pvals[i]), '{:.4f}'.format(qvals[i-ngenes])])
                else:
                    outlist.extend(['', ''])
                outlist.extend(['|'.join(map(str, g_nmut[:print_limit])),
                                '|'.join(map(str, g_rank[:print_limit])),
                                '|'.join(map(str, g_nmut_raw[:print_limit]))])
            else:
                outlist = [i+1, t, ','.join(g_in_t), '{:.6f}'.format(coef[i])]
                if (len(pvals) > 0) and (len(qvals) > 0):
                    outlist.extend(['{:.4f}'.format(pvals[i]), '{:.4f}'.format(qvals[i-ngenes])])
                else:
                    outlist.extend(['', ''])
                # outstr.append([i, coef[i], ','.join(g_in_t), ';'.join(map(str, nmut_in_t))])
                outlist.extend(['|'.join(map(str, g_nmut)),
                                '|'.join(map(str, g_rank)),
                                '|'.join(map(str, g_nmut_raw))])
            outstr.append(outlist)
        elif no_gene == False:
            outlist = [i+1, ont.genes[i], '{:.6f}'.format(coef[i]), '', '', nmut[i], nmut_sort[i], nmut_raw[i]]
            outstr.append(outlist) # for genes, p values are quite meaningless

    colnames = ['System index', 'System clixo/louvain name', 'Genes', 'Selection pressure', 'p', 'q',
                    'Mutation model input', 'Rank of model input', 'Mutation count']
    df_out = pd.DataFrame.from_records(outstr, columns=colnames)
    return df_out


def redistribute_gene_score(coef, mat, signal, out):
    '''
    calculate the redistributed gene score after the signal on genes was redistributed by regression
    :param coef: path to input file, a data frame, with dims [n_features, n_lambdas]
    :param mat: a binary sparse matrix in text format, with dims [n_genes, n_features]
    :param signal: a [n_genes, ] vector of the original signal on gene
    :return: 
    '''
    mat_feature_lambda = np.array(pd.read_table(coef, index_col=0))
    coo_gene_feature = np.loadtxt(mat).astype(int)
    vsignal = np.loadtxt(signal)
    with open(out, 'w') as ofh:
        for i in range(mat_feature_lambda.shape[1]):
            mat_gene_feature = coo_matrix((mat_feature_lambda[:, i][coo_gene_feature[:, 1]],
                                 (coo_gene_feature[:,0], coo_gene_feature[:,1])))
            # rescale this
            for x, y, d in zip(mat_gene_feature.row, mat_gene_feature.col, mat_gene_feature.data):
                if (d == 0) or (vsignal[x]==0):
                    continue
                ofh.write('{}\t{}\t{}\t{:.6f}\n'.format(i, x, y, d * vsignal[x] / float(np.sum(mat_gene_feature, axis=1)[x])))