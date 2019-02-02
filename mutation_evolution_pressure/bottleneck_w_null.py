import argparse, shutil
import multiprocessing as mp
from ddot import *
from utils import *

if __name__ == "__main__":
    par = argparse.ArgumentParser()
    par.add_argument('--ont', required=True)
    par.add_argument('--alt', required=True, help='ONP file')
    par.add_argument('--n', type=int, default=500, help='minimum number of ontology shuffling')
    par.add_argument('--fdr', type=float, default=0.3, help='fdr cutoff')
    par.add_argument('--log', action="store_true", help='if true, log transform mutation rate')
    par.add_argument('--save_more', action="store_true", help='if true, save the coefficients of random samples as a numpy matrix')
    args = par.parse_args()

    # local_path = '/nrnb/users/f6zheng/mutation_pressure/' + os.path.basename(args.ont) + '_' + os.path.basename(args.alt)
    # if not os.path.isdir(local_path):
    #     os.makedirs(local_path)

    os.system('grep gene ' + args.ont + ' | cut -f 2 | sort -u > all_genes')
    ont = Ontology.from_table(args.ont, clixo_format=True, is_mapping=lambda x:x[2]=='gene')
    ont.propagate('forward', inplace=True)
    ngenes = len(ont.genes)
    ont_conn = ont.connected()[:len(ont.genes), :].astype(int)

    #
    nmut_real = calculateMutationNumber(args.alt, ont.genes)
    if args.log:
        nmut_real = np.log(nmut_real + 1)

    # run lasso (randomized ones are run in parallel)
    real_result = run_lasso((ont_conn, nmut_real))

    batch_size = 100
    n_shuffle = max(args.n, int(len(ont.terms) / args.fdr))
    # n_shuffle = 10
    batch_ranges = np.arange(1, n_shuffle, batch_size)
    batch_ranges = np.append(batch_ranges, n_shuffle+1)
    print "Need {} permutations".format(n_shuffle)
    results_all, pvals, qvals = None, None, None
    nreject = np.zeros(len(ont.genes) +len(ont.terms),)
    nmut_rand = [np.random.permutation(nmut_real) for k in range(n_shuffle)]

    for i in range(len(batch_ranges)-1):
        # rand_onts = shuffleOntology(args.ont, local_path, (batch_ranges[i], batch_ranges[i+1]))
        pool = mp.Pool(processes=6)
        # map_args = [(r, nmut_real) for r in rand_onts]
        map_args = [(ont_conn, nmut_rand[k-1]) for k in range(batch_ranges[i], batch_ranges[i+1])]
        results = pool.map(run_lasso, map_args)
        results_stack = np.vstack(results)
        # if i > 0:
        #     results_all = np.append(results_all, results_stack, axis=0)
        # else:
        #     results_all = results_stack

        # calculate p-value on the fly
        real_result_rep = np.repeat(real_result[:, np.newaxis], batch_ranges[i+1]-batch_ranges[i], axis=1)
        nreject += np.sum(results_stack >= real_result_rep.T, axis=0)
        # calculate q-value
        # pvals = nreject / (batch_ranges[i+1]-1)


        # multitest = multipletests(pvals[ngenes:], alpha=0.3, method='fdr_bh', is_sorted=False)
        # n_sig = np.sum(multitest[0])
        # qvals = multitest[1]

        print "Finish {} permutations".format(batch_ranges[i+1]-1)
        # if n_sig == 0: # nothing is significant even with lower number of permutation
        #     break
    # don't do early termination
    pvals = nreject / n_shuffle
    multitest = multipletests(pvals[ngenes:], alpha=0.3, method='fdr_bh', is_sorted=False)
    qvals = multitest[1]

    tag = args.ont.replace('.ont', '') + '_' + os.path.basename(args.alt).replace('.onp', '')
    if args.log:
        tag = tag + '.log'
    if args.save_more:
        np.save(tag + '.pval', pvals)
        np.save(tag + '.coef', real_result)
        # np.save(tag + '.coef_rand', results_all)

    # print results
    outstr = printModuleProfile(real_result, ont, nmut_real, pvals, qvals)
    with open(tag + '.bn.txt', 'w') as fh:
        for s in outstr:
            fh.write('\t'.join(map(str,s)) + '\n')

    # # clean_up
    # shutil.rmtree(local_path)
