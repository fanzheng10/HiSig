import argparse
import numpy as np
from ddot import *
import Cluster
import Database

def parseCliXOparamter(f):
    pdict = {x:[] for x in ['a', 'b', 'm', 'z']}
    with open(f) as fh:
        for l in fh:
            info = l.strip().split()
            if info[0] in pdict:
                pdict[info[0]] = info[1:]
    return pdict

par = argparse.ArgumentParser()
par.add_argument('--o', required=True, help='ontology 3-col file')
par.add_argument('--s', required=True, help='similarity score file (numpy array)')
par.add_argument('--r', type=int, nargs=2, default=(100,700), help='the term size range to run CliXO')
par.add_argument('--p', required=True, help='parameter file; clixo parameters to test')
par.add_argument('--clixo_path', default='/cellar/users/f6zheng/Code/CliXO/clixo')
par.add_argument('--rename', action='store_true', help='whether need to reindex the term names')
args = par.parse_args()

ont = Ontology.from_table(args.o, clixo_format=True, is_mapping=lambda x:x[2]=='gene')
ont.propagate('forward', inplace=True)

# reindex terms by their size
if args.rename:
    term_names = sorted(ont.terms, key=lambda x: ont.term_sizes[ont.terms_index[x]])
    term_names = map(str, term_names)
    term_rename = {t: 'S' + str(len(term_names) - term_names.index(t)) for t in term_names}
    ont.rename(terms=term_rename, inplace=True)
    ont.to_table(args.o.replace('.ont', '.clean.ont'), clixo_format=False, header=False)
# ont_conn = ont.connected()

terms = [t for t in ont.terms if (ont.term_sizes[ont.terms_index[t]] >= args.r[0]) and (ont.term_sizes[ont.terms_index[t]] <= args.r[1])]
to_delete = []
for i in range(len(terms)):
    for ct in ont.parent_2_child[terms[i]]:
        to_delete.append(ct)
#     for j in range(i+1, len(terms)):
#         if ont_conn[len(ont.genes)+terms[i], len(ont.genes)+terms[j]]:
#             to_delete.append(terms[i])

terms = [t for t in terms if not t in to_delete]

# create score_file
geneind = Database.geneIndex()
rf_score = np.load(args.s)

clixo_pars = parseCliXOparamter(args.p)

input_scores = ['score_{}.txt'.format(t) for t in terms]

for i in range(len(input_scores)):
    n = 0
    dict_arr = {x:[] for x in ['a', 'b', 'm', 'z']}
    dirname = input_scores[i].replace('.txt', '')
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    Database.subset_pairwise_score(rf_score, [ont.genes[g] for g in ont.term_2_gene[terms[i]]], dirname+'/'+input_scores[i])

    for a, b, m, z in itertools.product(clixo_pars['a'], clixo_pars['b'], clixo_pars['m'], clixo_pars['z']):
        dict_arr['a'].append(a)
        dict_arr['b'].append(b)
        dict_arr['m'].append(m)
        dict_arr['z'].append(z)
        n += 1

    cmds = ["export as=({})".format(' '.join(dict_arr['a'])),
            "export bs=({})".format(' '.join(dict_arr['b'])),
            "export ms=({})".format(' '.join(dict_arr['m'])),
            "export zs=({})".format(' '.join(dict_arr['z']))]
    cmds.extend(["export a=${as[$((SGE_TASK_ID-1))]}",
                 "export b=${bs[$((SGE_TASK_ID-1))]}",
                 "export m=${ms[$((SGE_TASK_ID-1))]}",
                 "export z=${zs[$((SGE_TASK_ID-1))]}"])

    out = dirname+ '/'+input_scores[i].replace('.txt', '.$SGE_TASK_ID.out')
    cmds.append(' '.join([args.clixo_path, '-i', dirname+'/'+input_scores[i],
               '-a $a -b $b -m $m -z $z -s 0.3 > ' + out]))
    Cluster.qsub(cmds, fileName=input_scores[i].replace('.txt', '.clixo.sh'), mem=8, arrjob = (1, n), opts=["#$ -tc 50"], dry=True)


