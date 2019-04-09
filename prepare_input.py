import argparse
from ddot import *

def prepare_input(ont, sig, outf_conn, outf_sig, scut=2):
    '''
    
    :param ont: 
    :param signal: 
    :param outf: 
    :param scut: 
    :return: 
    '''
    ont = Ontology.from_table(ont, clixo_format=True, is_mapping = lambda x:x[2]=='gene')
    ont.propagate('forward', inplace=True)
    if len(ont.get_roots())> 1:
        ont.add_root('ROOT', inplace=True)

    with open(outf_conn, 'w') as ofh:
        for i in range(len(ont.genes)):
            g = ont.genes[i]
            terms = ont.gene_2_term[g]
            ofh.write('{}\t{}\n'.format(i, i))
            for t in terms:
                if ont.term_sizes[t] < scut:
                    continue
                ofh.write('{}\t{}\n'.format(i, t+len(ont.genes)))

    signal = {}
    with open(sig) as fh, open(outf_sig, 'w') as ofh:
        for l in fh:
            g, s = l.strip().split()
            s = float(s)
            signal[g] = s
        for g in ont.genes:
            if g in signal:
                ofh.write('{}\n'.format(signal[g]))
            else:
                ofh.write('0.0\n')

    with open('terms.txt', 'w') as fh:
        for i in range(len(ont.terms)):
            if ont.term_sizes[i] < scut:
                continue
            fh.write(ont.terms[i] +'\n')


if __name__ == "__main__":
    par = argparse.ArgumentParser()
    par.add_argument('--ont', required=True, help = 'an ontology file')
    par.add_argument('--sig', required=True, help = 'a text file for the signal on gene (leaf nodes)')
    par.add_argument('--out', nargs=2, required=True, help = 'output file')
    par.add_argument('--scut', default=2, type=int, help='cutoff; systems smaller than this will be removed')
    args = par.parse_args()


    prepare_input(args.ont, args.sig, args.out[0], args.out[1], args.scut)