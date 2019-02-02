import numpy as np
import pandas as pd
import argparse
import multiprocessing as mp
from utils import *

def removeSmallSystems(argv):
    conn = np.loadtxt(argv[0], dtype=int, delimiter="\t")
    # remove columns with small number of hits
    colSum = np.sum(conn, axis=0)
    conn_out = conn[:, colSum > args.cut]
    np.savetxt(argv[1], conn_out, fmt="%d", delimiter='\t')
    print('Finish ' + argv[0])

def removeSmallSystems_from_sparse(argv):
    conn = pd.read_table(argv[0], header=None, sep="\t")
    counts = conn[1].value_counts()
    keep_index = counts[counts > args.cut].index.tolist()
    keep_index.sort()
    conn_out = conn.loc[conn[1].isin(keep_index), :]
    dict_reindex = {keep_index[i]:i for i in range(len(keep_index))}
    records = []
    for i, row in conn_out.iterrows():
        records.append((row[0], dict_reindex[row[1]]))
    conn_out2 = pd.DataFrame.from_records(records)
    conn_out2.to_csv(argv[1], sep="\t", header=False, index=False)
    print('Finish ' + argv[0])

par = argparse.ArgumentParser()
par.add_argument('--cut', required=True, type=int, help= 'cutoff; systems smaller than this will be removed')
par.add_argument('--conn_fs', nargs='+', help='connectivity files')
par.add_argument('--sparse', action='store_true', help='if true, apply on sparse file')
par.add_argument('--keep_gene', action='s')
args = par.parse_args()


pool = mp.Pool(processes=10)

if args.sparse == False:
    inps = []
    for i in range(len(args.conn_fs)):
        inps.append((args.conn_fs[i], args.conn_fs[i].replace('conn', 'conn_gt{}'.format(args.cut))))
    pool.map(removeSmallSystems, inps)
else:
    inps = []
    for i in range(len(args.conn_fs)):
        inps.append((args.conn_fs[i], args.conn_fs[i].replace('conn', 'conn_gt{}'.format(args.cut))))
    pool.map(removeSmallSystems_from_sparse, inps)