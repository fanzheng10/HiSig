import argparse
from ddot import *


if __name__ == "__main__":
    par = argparse.ArgumentParser()
    par.add_argument('--ont', required=True, help = 'an ontology file')
    par.add_argument('--sig', required=True, help = 'a text file for the signal on gene (leaf nodes)')
    par.add_argument('--out', required=True, help = 'output file')
    par.add_argument('--scut', required=True, type=int, help='cutoff; systems smaller than this will be removed')
    args = par.parse_args()

    ont = Ontology.from_table(args.ont, clixo_format=True, is_mapping = lambda x:x[2]=='gene')
    ont.propagate('forward', inplace=True)

    with open(args.out, 'w') as ofh:
        for i in range(len(ont.genes)):
            g = ont.genes[i]
            terms = ont.gene_2_term[g]
            for t in terms:
                ofh.write('{}\t{}\n'.format(i, t+len(ont.genes)))

    signal = {}
    with open(args.sig) as fh, open('gene_signals.txt', 'w') as ofh:
        for l in fh:
            g, s = l.strip().split()
            s = float(s)
            signal[g] = s
        for g in ont.genes:
            if g in signal:
                ofh.write('{}\t{}\n'.format(g, signal[g]))
            else:
                ofh.write('{}\t0.0\n'.format(g))

    with open('terms.txt', 'w') as fh:
        for t in ont.terms:
            fh.write(t +'\n')