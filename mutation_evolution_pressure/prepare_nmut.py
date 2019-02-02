import numpy as np
import pandas as pd
import argparse
import multiprocessing as mp
from ddot import *

def calculateMutationNumber2(args):
    df_mut = pd.read_table(args[0], sep='\t', index_col=0, low_memory=False)
    genes_mutated = df_mut.index.tolist()
    nmut = np.zeros(len(genes))
    for i in range(len(genes)):
        if genes[i] in genes_mutated:
            nmut[i] = np.sum(df_mut.loc[genes[i], :].notnull())
    if args[1] != None:
        np.savetxt(args[1], nmut.T, fmt="%.3f")
    print "Finish " + args[1]
    # return nmut

if __name__ == "__main__":
    par = argparse.ArgumentParser()
    par.add_argument('--alt', nargs= '+', required=True, help='alteration file (ONP)')
    par.add_argument('--ont', required=True, help='the ontology file')
    par.add_argument('--mcore', type=int, default=1, help='if not 1, use multiple cores')
    # par.add_argument('--out', required=True, help='output files')
    args = par.parse_args()

    ont = Ontology.from_table(args.ont, clixo_format=True, is_mapping=lambda x:x[2]=='gene')
    genes = ont.genes

    # read the mutation fle
    if args.mcore==1:
        inp = (args.alt[0], os.path.basename(args.alt[0]) + '.y.txt')
        calculateMutationNumber2(inp)
    else:
        pool = mp.Pool(processes=min(args.mcore, mp.cpu_count()-2))

        inps = []
        for i in range(len(args.alt)):
            inps.append((args.alt[i], os.path.basename(args.alt[i]) + '.y.txt'))
        pool.map(calculateMutationNumber2, inps)