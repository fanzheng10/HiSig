import glob, argparse
import pandas as pd
from ddot import *

# def ont_conn_to_text(ontf, outf):
#     '''
#
#     :param ontf: file of ontology
#     :param outfile: output text file
#     :return:
#     '''
#     #TODO: through ddot takes too much memory.
#     ont = Ontology.from_table(ontf, clixo_format=True, is_mapping=lambda x:x[2]=='gene')
#     ont.propagate('forward', inplace=True)
#     ngenes = len(ont.genes)
#     ont_conn = ont.connected()[:len(ont.genes), :].astype(int)
#     np.savetxt(outf, ont_conn, fmt="%d", delimiter="\t")
#     # print('Finish '+args[1])

if __name__ == "__main__":
    par = argparse.ArgumentParser()
    par.add_argument('--ont', required=True, help = 'ontology file')
    par.add_argument('--out', required=True, help = 'output file')
    par.add_argument('--sparse', action='store_true', help = 'if true, write sparse form')
    args = par.parse_args()

    clixo_path = '/cellar/users/f6zheng/Code/CliXO/'
    ont_genes = args.ont + '.genes'
    if not os.path.isfile(ont_genes):
        os.system(clixo_path + '/ontologyTermStats {} genes > {}'.format(args.ont, ont_genes))
    df_genes = pd.read_table(ont_genes, sep="\t", header=None)


    # get all the genes and terms
    all_genes, all_terms, dict_gene_2_terms = [], [], {}
    for i, row in df_genes.iterrows():
        t = row[0]
        all_terms.append(t)
        t_genes = row[2].split(',')[:-1]
        for g in t_genes:
            all_genes.append(g)
            if not g in dict_gene_2_terms:
                dict_gene_2_terms[g] = []
            dict_gene_2_terms[g].append(t)
    all_genes = list(set(all_genes))
    all_genes.sort()
    all_terms.sort()

    term_indices = {all_terms[i]:i for i in range(len(all_terms))}

    with open(args.out, 'w') as fh:
        if not args.sparse:
            for i in range(len(all_genes)):
                g = all_genes[i]
                outarray = [0 for j in range(len(all_genes))]
                outarray[i] = 1
                outarray2 = [0 for k in range(len(all_terms))]
                for t in dict_gene_2_terms[g]:
                    outarray2[term_indices[t]] = 1
                outarray.extend(outarray2)
                outstr = "\t".join(map(str, outarray))
                fh.write(outstr + '\n')
        else:
            for i in range(len(all_genes)):
                fh.write('{}\t{}\n'.format(i, i))
                for t in dict_gene_2_terms[all_genes[i]]:
                    j = len(all_genes) + term_indices[t]
                    fh.write('{}\t{}\n'.format(i, j))


    with open(args.ont +'.terms', 'w') as fh:
        for t in all_terms:
            fh.write(t +'\n')