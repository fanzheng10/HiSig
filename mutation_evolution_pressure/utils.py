import os
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.stats import rankdata
from ddot import *
from sklearn.linear_model import lasso_path
from statsmodels.stats.multitest import *

def calculateMutationNumber(mutf, genes, outfile=None, gene_length=None, log_pc=None):
    '''
    calculate the number of alteration events of genes in cancer cohort
    :param mutf: mutation frequency; in in-house "ONP" format
    :param genes: a list of gene symbols
    :param outfile: if not None, write a txt file
    :param gene_length: if not None, being a dictionary of gene length
    :param log_pc: if not None, take a log transform and add this as pseudocount
    :return: a vector for the numeber of events of each gene
    '''
    df_mut = pd.read_table(mutf, sep='\t', index_col=0, low_memory=False)
    genes_mutated = df_mut.index.tolist()
    nmut = np.zeros(len(genes))
    for i in range(len(genes)):
        if genes[i] in genes_mutated:
            nmut[i] = np.sum(df_mut.loc[genes[i], :].notnull())

    if not gene_length == None:
        assert isinstance(gene_length, dict)
        len_median = np.median(np.array(gene_length.values()))
        for i in range(len(genes)):
            if genes[i] in gene_length:
                nmut[i] /= (gene_length[genes[i]] / len_median)

    if not log_pc == None:
        nmut = np.log(nmut + log_pc)

    if outfile != None:
        np.savetxt(outfile, nmut.T, fmt="%.3f")
    return nmut


def parseCliXOparamter(f):
    '''
    parse clixo parameters for sweeping
    :param f: a text flle; each row starts with parameter name (-a, -b, -m, -z), followed by values to test
    :return: a dictionary with key being parameter name and values being tested parameters
    '''
    pdict = {x:[] for x in ['a', 'b', 'm', 'z']}
    with open(f) as fh:
        for l in fh:
            info = l.strip().split()
            if info[0] in pdict:
                pdict[info[0]] = info[1:]
    return pdict


def shuffleOntology(ont_file, local_path, id_range = (1,2)):
    '''
    shuffle ontology, shuffling gene labels
    :param ont_file:
    :param n:
    :return:
    '''
    ali_path = '/cellar/users/f6zheng/tools/alignOntology'
    # os.system('grep gene ' + ont_file + ' | cut -f 2 | sort -u > all_genes')
    if not os.path.isdir(local_path):
        os.makedirs(local_path)

    rand_onts = []
    for i in range(id_range[0], id_range[1]):
        rand_dict = local_path+ '/rand_dict_{}'.format(i)
        os.system('shuf all_genes | paste all_genes - > {}'.format(rand_dict))
        rand_ont_f = local_path+ '/rand_ont_{}'.format(i)
        os.system(ali_path + '/replace_words_using_dicts.pl {} {} > {}'.format(ont_file, rand_dict, rand_ont_f))
        rand_onts.append(rand_ont_f)
    return rand_onts


def maxfrac(coef_all, n_alpha):
    '''
    calculate the maximum contribution of each term's coefficient
    :param coef_all: a numpy matrix, first dimension is the number of terms
    :param n_alpha: number of different lambda values tested in Lasso
    :return: the maximum contribution of each term; this is like a score of each term
    '''
    coef_all_frac = np.zeros_like(coef_all)

    for i in range(n_alpha):
        if np.sum(coef_all[:, i]) == 0:
            continue
        coef_all_frac[:, i] = coef_all[:, i] / np.sum(coef_all[:, i])

    coef_max_frac = np.max(coef_all_frac, axis=1)
    return coef_max_frac

def create_data_blocks(dims_big, dims_small, data, data_num, data_dims):
    mat = -1 * np.ones((dims_big[0] * dims_small[0], dims_big[1] * dims_small[1]))
    columns = [d[0] for d in data_dims]
    assert data_num in data.columns
    for c in columns:
        assert c in data.columns
    for i, row in data.iterrows():
        num = row[data_num]
        coord = [data_dims[k].index(row[columns[k]]) - 1 for k in range(len(columns))]
        x = coord[0] * dims_small[0] + coord[2]
        y = coord[1] * dims_small[1] + coord[3]
        mat[x, y] = num
    return mat


def run_lasso(args):
    '''
    args condense multiple parameters for the need to pass to mp.map
    :param args: args[0] is the connectivity file, args[1] is mutation count
    :return: the bottleneck score of terms
    '''

    # n_mutation_of_genes = calculateMutationNumber(args[1], ont.genes)
    clf_lasso_path = lasso_path(args[0], args[1], positive=True)
    coef_max_frac = maxfrac(clf_lasso_path[1], len(clf_lasso_path[0]))
    return coef_max_frac


def printModuleProfile(coef, ont, nmut, nmut_raw, pvals = (), qvals = (), top=None, print_limit=50, sort_by_pval = False, min_term_size=2, no_gene=False):
    '''

    :param coef: selection score
    :param ont: ontology object
    :param nmut: input of the lasso model
    :param nmut_raw: number of mutations for each gene
    :param pvals: p values of each system
    :param qvals: q values
    :param top: if not None, only print top N results
    :param print_limit: for terms bigger than this limit do not print the exact gene names
    :param sort_by_pval: if True, sort results by p-value instead of bottleneck score
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

def ont_conn_to_text(ontf, outfile):
    '''

    :param ontf: file of ontology
    :param outfile: output text file
    :return:
    '''
    ont = Ontology.from_table(ontf, clixo_format=True, is_mapping=lambda x:x[2]=='gene')
    ont.propagate('forward', inplace=True)
    ngenes = len(ont.genes)
    ont_conn = ont.connected()[:len(ont.genes), :].astype(int)
    np.savetxt(outfile, ont_conn, fmt="%d", delimiter="\t")

def removeSmallSystems(inpf, outf, cut=5):
    conn = np.loadtxt(inpf)
    # remove columns with small number of hits
    colSum = np.sum(conn, axis=1)
    conn_out = conn[:, colSum > cut]
    np.savetxt(outf, conn_out, fmt="%d")
    print('Finish ' + inpf)