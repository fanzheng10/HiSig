import numpy as np
import argparse
from ddot import *
from utils import *
from sklearn.metrics import average_precision_score, roc_auc_score

def rerank_genes(coef, conn, genes):
    '''

    :param coef: a file of [n_genes, n_lambda] matrix
    :param conn: a file of [n_genes, n_features] sparse matrix
    :param genes: a list of gene names
    :return: rank_gl_alllambda: list of lists, [n_lambda, n_genes]
    '''
    mat_feature_lambda = np.array(pd.read_table(coef, index_col=0))
    coo_gene_feature = np.loadtxt(conn).astype(int)

    predictions = []
    for i in range(mat_feature_lambda.shape[1]):
        mat_gene_feature = coo_matrix(
            (mat_feature_lambda[:, i][coo_gene_feature[:, 1]], (coo_gene_feature[:, 0], coo_gene_feature[:, 1])))
        # prediction is simply row sum
        if np.sum(mat_gene_feature.data) == 0:
            continue
        prediction = np.sum(mat_gene_feature, axis=1)
        predictions.append(prediction)

    return predictions
    # rank_gl_alllambda = []
    # for pred in predictions[::-1]:
    #     rank_gl = np.array([genes[gi] for gi in np.argsort(np.asarray(pred).ravel())[::-1]])
    #     rank_gl_alllambda.append(rank_gl)
    #
    # return rank_gl_alllambda


if __name__ == "__main__":
    par = argparse.ArgumentParser()
    par.add_argument('--ont', required=True, help='the ontology file')
    par.add_argument('--coef', required=True, help='sparse matrix, for gene-system membership')
    par.add_argument('--conn', required=True, help='the output of R; if there are multiple files concatenate them')
    par.add_argument('--pos_genes', required=True, help='a text file containing gene names of the positive set')
    par.add_argument('--metric', choices=['ROC', 'PRC'],  help='the metric to select the best lambda - choose between ROC and PRC')
    par.add_argument('--out', required=True, help='output file')
    args = par.parse_args()

    ont = Ontology.from_table(args.ont, clixo_format=True, is_mapping = lambda x:x[2] =='gene')
    ont.propagate('forward', inplace=True)
    predictions_all_lambda = rerank_genes(args.coef, args.conn, ont.genes)

    pos_genes = []
    with open(args.pos_genes) as fh:
        for l in fh:
            pos_genes.append(l.strip())

    y_true = np.zeros(len(ont.genes), )
    for g in pos_genes:
        if g in ont.genes:
            y_true[ont.genes.index(g)] = 1

    auc_array = []
    if args.metric == 'PRC':
        for i in range(len(predictions_all_lambda)):
            auc_array.append(average_precision_score(y_true, predictions_all_lambda[i]))
    if args.metric == 'ROC':
        for i in range(len(predictions_all_lambda)):
            auc_array.append(roc_auc_score(y_true, predictions_all_lambda[i]))

    best_lam = np.argmax(auc_array)
    best_auc = np.max(auc_array)
    predictions_best = predictions_all_lambda[best_lam]

    gene_ind_order = np.argsort(np.asarray(predictions_best.ravel())[0])[::-1]
    with open(args.out, 'w') as fh:
        fh.write('# Best AU{}: {:4f}; found at lambda index {:d}\n'.format(args.metric, best_auc, best_lam))
        for i in range(len(gene_ind_order)):
            if predictions_best[gene_ind_order[i]] == 0:
                break
            fh.write('{}\t{:.3f}\t{}\n'.format(ont.genes[gene_ind_order[i]], float(predictions_best[gene_ind_order[i]]), int(y_true[gene_ind_order[i]])))