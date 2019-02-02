import os,argparse, itertools
import Cluster

par = argparse.ArgumentParser()
par.add_argument('--conn', nargs='+', help = "all the structure matrices of ontology")
par.add_argument('--nmut', required=True, help = "a list of all mutation files")
par.add_argument('--sh', required=True, help = "the output sh file")
par.add_argument('--rand', help = "if none, number of permutation")
args = par.parse_args()

yfile = args.nmut
src = '/cellar/users/f6zheng/Code/HumanOnt/CaCHe_pipelines/mutation_evolution_pressure/R/biglasso_once_1core.R'

if args.rand != None:
    n_rand = int(args.rand)
    tmpdir = '/nrnb/users/f6zheng/biglasso/'
    cmds = Cluster.arraryJobConstructor(('shuf', 'conn'), (map(str, list(range(n_rand))),args.conn))
    outfile = yfile + '.$conn.bn.$shuf'
    njob = len(args.conn) * n_rand
    cmds.append('R -f {} --args $conn {} {} 1 1'.format(src, yfile, outfile))

else:
    cmds = Cluster.arraryJobConstructor(('conn',), (args.conn,))
    outfile = yfile + '.$conn.bn'
    njob = len(args.conn)
    cmds.append('R -f {} --args $conn {} {} 1 0'.format(src, yfile, outfile))

# different jobs copying files simultaneously cause severe I/O issue

# cmds.append('odir=$PWD')

# create a local directory
# somehow very slow if it is a nrnb location. Just use Data location is fine
# ldir_path = '/nrnb/users/f6zheng/biglasso_once/' + str(os.getpid()) + '.$SGE_TASK_ID'

# ldir_path = '/cellar/users/f6zheng/Data/tmp/biglasso_once/' + str(os.getpid()) + '.$SGE_TASK_ID'
# cmds.append("mkdir {}".format(ldir_path))
# # copy the input files to this path
# cmds.append('cp $conn $nmut {}'.format(ldir_path))
# # change directory
# cmds.append('cd '+ldir_path)
# run R



# copy results back
# cmds.append('cp -r $nmut.$conn.bn* $odir')
# # remove the local directory
# cmds.append('chdir $odir')
# cmds.append('rm -r '+ldir_path)

Cluster.qsub(cmds, mem=6, fileName=args.sh, arrjob=(1, njob), opts=["#$ -tc 500"], dry=True)

