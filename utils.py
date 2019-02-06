import numpy as np
import pandas as pd
import itertools
from scipy.sparse import *
from scipy.stats import rankdata


def printModuleProfile(coef, ont, signal, pvals, qvals,
                       signal2 = None, top=None, print_limit=25, min_term_size=2, no_gene=False, signal2_name='Mutation_count'):
    '''

    :param coef: selection score
    :param ont: ontology object
    :param signal: input of the lasso model
    :param signal2: number of mutations for each gene
    :param pvals: p values of each system
    :param qvals: q values
    :param top: if not None, only print top N results
    :param print_limit: for terms bigger than this limit do not print the exact gene names
    :param no_gene: if True, the results only contain terms
    :return:
    '''
    ngenes = len(ont.genes)
    idx = np.argsort(pvals)

    if top!=None:
        idx = idx[:top]
    outstr = []
    signal_sort = rankdata(-signal, method='max') # tied ranks with receive the maximum index

    # # multiple-test correction
    for i in idx:
        if coef[i] == 0:
            break
        if i >= ngenes: # if the system is a term
            t = ont.terms[i-ngenes]
            if ont.term_sizes[i-ngenes] <  min_term_size:
                continue

            g_in_t_ind = [gi for gi in ont.term_2_gene[t]]
            g_in_t_ind_sort = sorted(g_in_t_ind, key = lambda x:signal_sort[x])

            g_in_t = [ont.genes[gi] for gi in g_in_t_ind_sort]

            g_signal = [signal[gi] for gi in g_in_t_ind_sort]
            g_rank = [signal_sort[gi] for gi in g_in_t_ind_sort]
            g_signal_raw = [signal2[gi] for gi in g_in_t_ind_sort] # TODO: right now if not having this vector, the program will fail. Will make this optional

            if len(g_in_t) > print_limit:
                outlist = [i+1, t, '{} genes including: {}'.format(len(g_in_t), ','.join(g_in_t[:print_limit])), '{:.6f}'.format(coef[i])]
                if (len(pvals) > 0) and (len(qvals) > 0):
                    outlist.extend(['{:.4f}'.format(pvals[i]), '{:.4f}'.format(qvals[i-ngenes])])
                else:
                    outlist.extend(['', ''])
                outlist.extend(['|'.join(map(str, g_signal[:print_limit])),
                                '|'.join(map(str, g_rank[:print_limit])),
                                '|'.join(map(str, g_signal_raw[:print_limit]))])
            else:
                outlist = [i+1, t, ','.join(g_in_t), '{:.6f}'.format(coef[i])]
                if (len(pvals) > 0) and (len(qvals) > 0):
                    outlist.extend(['{:.4f}'.format(pvals[i]), '{:.4f}'.format(qvals[i-ngenes])])
                else:
                    outlist.extend(['', ''])
                # outstr.append([i, coef[i], ','.join(g_in_t), ';'.join(map(str, signal_in_t))])
                outlist.extend(['|'.join(map(str, g_signal)),
                                '|'.join(map(str, g_rank)),
                                '|'.join(map(str, g_signal_raw))])
            outstr.append(outlist)
        elif no_gene == False:
            outlist = [i+1, ont.genes[i], '{:.6f}'.format(coef[i]), '', '', signal[i], signal_sort[i], signal2[i]]
            outstr.append(outlist) # for genes, p values are quite meaningless # TODO: recosnider this decision

    colnames = ['System_index', 'System_name', 'Genes', 'Selection_pressure', 'p', 'q',
                    'Model_input', 'Rank_of_model_input', signal2_name]
    df_out = pd.DataFrame.from_records(outstr, columns=colnames)
    return df_out


def redistribute_gene_score(coef, mat, signal, exponential=False):
    '''
    calculate the redistributed gene score after the signal on genes was redistributed by regression
    :param coef: path to input file, a data frame, with dims [n_features, n_lambdas]
    :param mat: a binary sparse matrix in text format, with dims [n_genes, n_features]
    :param signal: a [n_genes, ] vector of the original signal on gene
    :return: outputs: a dictionary, lambda being the key, and each value is a [n_genes, n_systems] (in sparse form, so 3 lists)
    '''
    mat_feature_lambda = np.array(pd.read_table(coef, index_col=0))
    coo_gene_feature = np.loadtxt(mat).astype(int)
    vsignal = np.loadtxt(signal)
    # outputs = ([], [], [], [])
    outputs = {}

    for i in range(mat_feature_lambda.shape[1]):
        mat_gene_feature = coo_matrix((mat_feature_lambda[:, i][coo_gene_feature[:, 1]],
                             (coo_gene_feature[:,0], coo_gene_feature[:,1])))
        # rescale this
        scale = np.sum(mat_gene_feature, axis=1)
        idx = (mat_gene_feature.data > 0) & (vsignal[mat_gene_feature.row] >0)
        if np.sum(idx) == 0:
            continue
        # out = mat_gene_feature.data[idx] * vsignal[mat_gene_feature.row[idx]] / scale[0]
        out = mat_gene_feature.data[idx] / np.asarray(scale[mat_gene_feature.row[idx]].T)[0] # note that I decided to output the rescaled beta, not yet multiply to signal
        if exponential:
            mat_gene_feature_2 = coo_matrix((np.exp(out), (mat_gene_feature.row[idx], mat_gene_feature.col[idx])))
            scale_2 = np.sum(mat_gene_feature_2, axis=1)
            idx2 = (mat_gene_feature_2.data > 0) & (vsignal[mat_gene_feature_2.row] >0)
            if np.sum(idx2) == 0:
                continue
            out2 = mat_gene_feature_2.data[idx2] / np.asarray(scale_2[mat_gene_feature_2.row[idx2]].T)[0]
            outputs[i] = (mat_gene_feature_2.row[idx2], mat_gene_feature_2.col[idx2], out2)
        else:
            outputs[i] = (mat_gene_feature.row[idx], mat_gene_feature.col[idx], out)

    return outputs

def feature_best_lambda(coef):
    '''
    for each feature, get the lambda when it achieves the highest importance
    :param coef: a [n_feature, n_lambda] matrix
    :return: two vector; feature index and the best lambda of these feature. Features that are never positive are omitted
    '''
    mat_feature_lambda = np.array(pd.read_table(coef, index_col=0))
    idx = np.sum(mat_feature_lambda, axis=0) > 0
    mat_feature_lambda_scaled = mat_feature_lambda[:, idx] / (np.sum(mat_feature_lambda, axis=0)[idx])
    idy = np.where(np.sum(mat_feature_lambda_scaled, axis=1) != 0)[0]
    mat_feature_argmax = np.argmax(mat_feature_lambda_scaled, axis=1)
    return idy, mat_feature_argmax[idy] + (mat_feature_lambda.shape[1] -np.sum(idx))

def term_specific_rank(ont, coef_adjust, term_names, term_best_lambda, outf=None, print_limit=25):
    '''
    
    :param ont: an Ontology object  
    :param coef_adjust: adjusted regression coefficient, a dictionary that is the output from the {redistribute_gene_score} function
    :param term_names: terms names of interest
    :param term_best_lambda: a dictionary for terms' best lambda
    :param outf: output dataframe path
    :return: df_out: output dataframe
    '''
    # df_coef_adjust = pd.DataFrame.from_dict({'gene':coef_adjust[0], 'feature':coef_adjust[1], 'weight':coef_adjust[2], 'lambda':coef_adjust[3]})
    out_records = []
    for t in term_names:
        tid = ont.terms.index(t)
        lam = term_best_lambda[tid+len(ont.genes)]
        coef_adjust_i = coef_adjust[lam]
        df_coef_adjust = pd.DataFrame.from_dict(
            {'gene': coef_adjust_i[0], 'feature': coef_adjust_i[1], 'weight': coef_adjust_i[2]})

        df_coef_sub = df_coef_adjust.loc[(df_coef_adjust['feature'] == tid+ len(ont.genes)), :]
        df_coef_sub.sort_values(by='weight', ascending=False, inplace=True)
        if df_coef_sub.shape[0] > print_limit:
            df_coef_sub = df_coef_sub.iloc[:print_limit, :]
        df_coef_sub['gene_name'] = np.array([ont.genes[x] for x in df_coef_sub['gene'].tolist()])
        df_coef_sub = df_coef_sub.round(3)
        df_coef_sub['weight'] = df_coef_sub['weight'].astype(str)
        out_records.append((t, '|'.join(df_coef_sub['gene_name'].tolist()), '|'.join(df_coef_sub['weight'].tolist())))
    df_out = pd.DataFrame.from_records(out_records)
    if outf!=None:
        df_out.to_csv(outf, sep='\t', index=False, header=False)
    return df_out


def estimate_nsamples_per_term(ont, coef_adjust_expo, term_names, term_best_lambda, tumor_profile, signal=None, outf=None):
    '''
    to estimate number of samples mutated by each system
    :param ont: an Ontology object
    :param coef_adjust: "exponentially adjusted" regression coefficient, a dictionary that is the output from the {redistribute_gene_score} function
    :param term_names: terms names of interest
    :param test_best_lambda: a dictionary for terms' best lambda
    :param tumor_profile: a [n_gene, n_patients] data frame (in-house ONP format)
    :param signal #TODO: think about this
    :return: df_out: output dataframe
    '''
    dict_system_count_inferred = {}
    for tumid in tumor_profile.columns.tolist():
        genes_mutated_in_sample = []
        for g in tumor_profile[tumid].notnull().index.tolist():
            if g in ont.genes:
                genes_mutated_in_sample.append(ont.genes.index(g))
        for term_name in term_names:
            term_id = ont.terms.index(term_name)
            lam = term_best_lambda[term_id + len(ont.genes)]
            coef_adjust_sel = coef_adjust[lam] # [n_genes, n_system] rescaled beta matrix
            # genes_in_term = [ont.genes[gid] for gid  in ont.term_2_gene[term_name]]
            gene_ids_in_term = ont.term_2_gene[term_name]
            genes_hit = set(genes_mutated_in_sample).intersection(gene_ids_in_term)
            multiplicon = 1
            for g_hit in genes_hit:
                beta_scaled = coef_adjust_sel[g_hit, len(ont.genes) + term_id]
                if coef_adjust_sel[g_hit, len(ont.genes) + term_id] == 1: # this gene is only mutated in this system
                    # TODO: wrong, still in sparse form
                    multiplicon = 0
                    break
                else:
                    multiplicon *= 1-coef_adjust_sel

            genes_activated = list(set(tumor_profile[tumid
                                       ].notnull().index.tolist()).intersection(genes_in_term))



    for t in term_names:
        tid = ont.terms.index(t)
        lam = term_best_lambda[tid+len(ont.genes)]
        coef_adjust_i = coef_adjust[lam] # [n_genes, n_system] rescaled beta matrix







