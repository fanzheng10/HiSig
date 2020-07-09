import pandas as pd
import numpy as np
import argparse

def prepare_input(ont, sig, outf_conn, outf_sig):
    '''
    
    :param ont: 
    :param sig: 
    :param outf_conn: 
    :param outf_sig: 
    :param file_exist: 
    :return: 
    '''
    df = pd.read_csv(ont, '\t', header=None)
    df_terms = df.loc[df[2] != 'gene', :]
    df_genes = df.loc[df[2] == 'gene', :]

    genes = sorted(set(df_genes[1].tolist()))
    terms = sorted(set(df[0].tolist() + df_terms[1].tolist()))
    genes_idx = {genes[i]:i for i in range(len(genes))}
    terms_idx = {terms[i]:i for i in range(len(terms))}

    mat_terms = np.eye(len(terms))
    mat_gene2term = np.zeros((len(genes), len(terms)))
    for i, row in df_terms.iterrows():
        mat_terms[terms_idx[row[1]], terms_idx[row[0]]] = 1
    for i, row in df_genes.iterrows():
        mat_gene2term[genes_idx[row[1]], terms_idx[row[0]]] = 1

    # propagate
    while True:
        mat_gene2term_new = np.dot(mat_gene2term, mat_terms)
        if np.sum(mat_gene2term_new > 0) == np.sum(mat_gene2term > 0):
            mat_gene2term = mat_gene2term_new
            break
        else:
            mat_gene2term = mat_gene2term_new

    with open(outf_conn, 'w') as ofh:
        for i in range(len(genes)):
            ofh.write('{}\t{}\n'.format(i, i))
        row, col = np.where(mat_gene2term > 0)
        for i in range(len(row)):
            ofh.write('{}\t{}\n'.format(row[i], len(genes) + col[i]))

    with open('genes.txt', 'w') as fh:
        for i in range(len(genes)):
            fh.write(genes[i] + '\n')
    with open('terms.txt', 'w') as fh:
        for i in range(len(terms)):
            fh.write(terms[i] + '\n')

    signal = {}
    with open(sig) as fh, open(outf_sig, 'w') as ofh:
        for l in fh:
            g, s = l.strip().split()
            s = float(s)
            signal[g] = s
        for g in genes:
            if g in signal:
                ofh.write('{}\n'.format(signal[g]))
            else:
                ofh.write('0.0\n')


if __name__ == "__main__":
    par = argparse.ArgumentParser()
    par.add_argument('--ont', required=True, help = 'an ontology file')
    par.add_argument('--sig', required=True, help = 'a text file for the signal on gene (leaf nodes)')
    par.add_argument('--out', required=True, help = 'output prefix')
    par.add_argument('--exist', action='store_true', help='if true, skip generating sparse matrix (redundant)')
    args = par.parse_args()


    prepare_input(args.ont, args.sig, args.out+ '_conn.txt', args.out + '_signal.txt')