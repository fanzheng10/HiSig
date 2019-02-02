import itertools
import Cluster
from ddot import *
import matplotlib.pyplot as plt
import seaborn as sns

# prepare before glmnet

def writegenes(n, sdir):
    for i in range(n):
        ont = Ontology.from_table(sdir+'/nest{}_till_0.15.ont'.format(i+1),
                                 clixo_format=True, is_mapping=lambda x:x[2]==True)
        with open('fuse{}/nest{}_gene.txt'.format(i+1, i+1), 'w') as fh:
            for g in ont.genes:
                fh.write(g+'\n')
    
def write_lm_analysis_input(n, df_diff, ctypes):
    diff_all = []
    for j in range(len(ctypes)):
        diff = np.array(df_diff['{}_diff_log_ps5'.format(ctypes[j])])
        diff_all.append(diff)
    
    for i in range(n):
        gf = 'fuse{}/nest{}_gene.txt'.format(i+1, i+1)
        if not os.path.isfile(gf):
            continue
        genes = [l.strip() for l in open(gf).readlines()]
        ind_dic = {k:v for k,v in zip(df_diff.index.tolist(), np.arange(df_diff.shape[0]))}
        ind = np.array([ind_dic[g] if g in ind_dic else -1 for g in genes]) 
        for j in range(len(ctypes)):
            y = np.zeros(len(genes),)
            y[ind>=0] = diff_all[j][ind[ind>=0]]
            np.savetxt('fuse{}/build_20181101_{}_diff_nneg_log_ps5.txt'.format(i+1, ctypes[j]), y, 
                       fmt="%.3f")

            
def write_mut_observed(n, df_obs, ctypes):
    for i in range(n):
        gf = 'fuse{}/nest{}_gene.txt'.format(i+1, i+1)
        if not os.path.isfile(gf):
            continue
        genes = [l.strip() for l in open(gf).readlines()]
        ind_dic = {k:v for k,v in zip(df_obs.index.tolist(), np.arange(df_obs.shape[0]))}
        ind = np.array([ind_dic[g] if g in ind_dic else -1 for g in genes]) 
        for j in range(len(ctypes)):
            y = np.zeros(len(genes),)
            y[ind>=0] = np.array(df_obs[ctypes[j] + '_observed'])[ind[ind>=0]]
            np.savetxt('fuse{}/build_20181101_{}_observed.txt'.format(i+1, ctypes[j]), y, fmt='%d')


# prepare array jobs
            
            
def batch_prepareconn(n, sdir):
    ont = sdir + '/nest$SGE_TASK_ID\_till_0.15.ont'
    cmd = ['python', '/cellar/users/f6zheng/Code/HumanOnt/CaCHe_pipelines/mutation_evolution_pressure/prepare_conn.py',
          '--ont', sdir + '/nest$SGE_TASK_ID\_till_0.15.ont',
          '--out', 'fuse$SGE_TASK_ID/nest$SGE_TASK_ID\_conn.txt',
          '--sparse']
    cmd = ' '.join(cmd)
    Cluster.qsub([cmd], mem=32, fileName = 'prepare_conn.sh', arrjob=(1, n), short=True, dry=True)

    
def batch_glmnet(n, ctypes, labels, outname='glmnet_all.sh', tc=100):
    if isinstance(n, int):
        ids = [i+1 for i in range(n)]
    else:
        ids = n
    
    assert len(labels) == 2
    
    dirs = []
    args1 = []
    args2 = []
    args3 = []

    for c,i in itertools.product(ctypes, ids):
        dirs.append('fuse{}'.format(i))
        args1.append('nest{}_conn.txt'.format(i))
        args2.append('build_20181101_{}_{}.txt'.format(c, labels[0]))
        args3.append('nest{}.{}_{}.glmnet-ms-impact'.format(i, c, labels[1]))


    cmds = []
    cmds.append('export folders=({})'.format(' '.join(dirs)))
    cmds.append('export args1=({})'.format(' '.join(args1)))
    cmds.append('export args2=({})'.format(' '.join(args2)))
    cmds.append('export args3=({})'.format(' '.join(args3)))

    cmds.append('export folder=${folders[$((SGE_TASK_ID - 1))]}')
    cmds.append('export arg1=${args1[$((SGE_TASK_ID - 1))]}')
    cmds.append('export arg2=${args2[$((SGE_TASK_ID - 1))]}')
    cmds.append('export arg3=${args3[$((SGE_TASK_ID - 1))]}')

    cmds.append('cd $folder')
    Rcmd = ['R', '-f /cellar/users/f6zheng/Code/HumanOnt/CaCHe_pipelines/mutation_evolution_pressure/R/glmnet.R',
          '--args $arg1 $arg2 $arg3 2']
    Rcmd = ' '.join(Rcmd)
    cmds.append(Rcmd)
    Cluster.qsub(cmds, fileName =outname, 
                 cpu=8, mem=8, arrjob=(1,len(dirs)), opts=['#$ -tc {}'.format(tc)], dry=True)

def batch_glmnet_cv(n, ctypes, labels, outname='glmnet_all.sh', tc=100, cv=5):
    if isinstance(n, int):
        ids = [i+1 for i in range(n)]
    else:
        ids = n
    assert len(labels) == 2
    
    dirs = []
    args1 = []
    args2 = []
    args3 = []
    args4 = []

    for c,i in itertools.product(ctypes, ids):
        dirs.append('fuse{}'.format(i))
        args1.append('nest{}_conn.txt'.format(i))
        args2.append('nest{}_conn_gt4.txt'.format(i))
        args3.append('build_20181101_{}_{}.txt'.format(c, labels[0]))
        args4.append('nest{}.{}_{}.glmnet-cv'.format(i, c, labels[1]))


    cmds = []
    cmds.append('export folders=({})'.format(' '.join(dirs)))
    cmds.append('export args1=({})'.format(' '.join(args1)))
    cmds.append('export args2=({})'.format(' '.join(args2)))
    cmds.append('export args3=({})'.format(' '.join(args3)))
    cmds.append('export args4=({})'.format(' '.join(args4)))

    cmds.append('export folder=${folders[$((SGE_TASK_ID - 1))]}')
    cmds.append('export arg1=${args1[$((SGE_TASK_ID - 1))]}')
    cmds.append('export arg2=${args2[$((SGE_TASK_ID - 1))]}')
    cmds.append('export arg3=${args3[$((SGE_TASK_ID - 1))]}')
    cmds.append('export arg4=${args4[$((SGE_TASK_ID - 1))]}')

    cmds.append('cd $folder')
#     cmds.append('python ~/Code/HumanOnt/CaCHe_pipelines/mutation_evolution_pressure/prepare_truncate_conn.py --cut 4 --conn_fs $arg1 --sparse')
    Rcmd = ['R', '-f /cellar/users/f6zheng/Code/HumanOnt/CaCHe_pipelines/mutation_evolution_pressure/R/glmnet_cv.R',
          '--args $arg1 $arg2 $arg3 $arg4 {}'.format(cv)]
    Rcmd = ' '.join(Rcmd)
    cmds.append(Rcmd)
    Cluster.qsub(cmds, fileName =outname, 
                 cpu=8, mem=8, arrjob=(1,len(dirs)), opts=['#$ -tc {}'.format(tc)], dry=True)

    
def batch_summary(n, ctypes, sdir, labels, outname = 'summary_all.sh', tc=500):
    if isinstance(n, int):
        ids = [i+1 for i in range(n)]
    else:
        ids = n
    
    assert len(labels) == 5
    dirs = []
    args = [[], [], [], [], []]

    for c,i in itertools.product(ctypes, ids):
        dirs.append('fuse{}'.format(i))
        args[0].append(sdir + '/nest{}_till_0.15.ont'.format(i))
        args[1].append('nest{}.{}_{}.{}.txt'.format(i,c, labels[0], labels[1]))
        args[2].append('build_20181101_{}_{}.txt'.format(c, labels[2]))
        args[3].append('build_20181101_{}_{}.txt'.format(c, labels[3]))
        args[4].append('nest{}.{}_{}.{}.txt'.format(i,c, labels[0], labels[4]))


    cmds = []
    cmds.append('export folders=({})'.format(' '.join(dirs)))
    cmds.append('export args0=({})'.format(' '.join(args[0])))
    cmds.append('export args1=({})'.format(' '.join(args[1])))
    cmds.append('export args2=({})'.format(' '.join(args[2])))
    cmds.append('export args3=({})'.format(' '.join(args[3])))
    cmds.append('export args4=({})'.format(' '.join(args[4])))

    cmds.append('export folder=${folders[$((SGE_TASK_ID - 1))]}')
    cmds.append('export arg0=${args0[$((SGE_TASK_ID - 1))]}')
    cmds.append('export arg1=${args1[$((SGE_TASK_ID - 1))]}')
    cmds.append('export arg2=${args2[$((SGE_TASK_ID - 1))]}')
    cmds.append('export arg3=${args3[$((SGE_TASK_ID - 1))]}')
    cmds.append('export arg4=${args4[$((SGE_TASK_ID - 1))]}')

    cmds.append('cd $folder')

    cmd = ['python', '/cellar/users/f6zheng/Code/HumanOnt/CaCHe_pipelines/mutation_evolution_pressure/r_output_parser.py',
          '--ont', '$arg0',
          '--rout', '$arg1',
           '--nmut', '$arg2', '$arg3',
           '--index1', 
          '--out', '$arg4',
          '--sort_by_pval',
          '--min_term_size 4']

    cmd = ' '.join(cmd)
    cmds.append(cmd)
    Cluster.qsub(cmds, fileName =outname, mem=16,
                    arrjob=(1,len(dirs)), opts = ['#$ -tc {}'.format(tc)], short=True, dry=True)

# for the heatmap

def combine_par(df, pars =['a', 'b', 'm', 'z']):
    df[pars[0]] = df[pars[0]].astype(str)
    df['_'.join(pars)] = df[pars[0]].astype(str)
    for p in pars[1:]:
        df[p] = df[p].astype(str)
        df['_'.join(pars)] += '_'
        df['_'.join(pars)] += df[p]
    return df

def create_data_blocks(dims_big, dims_small, data, data_num, data_dims):
    mat = -1 * np.ones((dims_big[0] * dims_small[0], dims_big[1] * dims_small[1]))
    columns = [d[0] for d in data_dims]
    assert data_num in data.columns
    for c in columns:
        assert c in data.columns
    for i, row in data.iterrows():
        num = row[data_num]
        coord = [data_dims[k].index(row[columns[k]]) - 1 for k in range(len(columns))]
        x = coord[0] * dims_small[0] + coord[2]
        y = coord[1] * dims_small[1] + coord[3]
        mat[x, y] = num
    return mat

def get_nsys(n, df_config, sdir, size_cut = (2, 20000)):
    dic_size = {}
    for i in range(n):
        df = pd.read_table(sdir+'/nest{}_till_0.15.ont.genes'.format(i+1), sep='\t', header=None)
        df = df.loc[(df[1]>=size_cut[0]) & (df[1]<=size_cut[1])]
        dic_size[i] = df.shape[0]
    df_config['Size'] = pd.Series(dic_size)
    return df_config
    
    
def count_nsigsys(n, df_config, ctypes, label, sdir, q_cut=[0.3]):
    df_out = df_config.copy()
    for q,c in itertools.product(q_cut, ctypes):
        n_sig = []
        for i in range(1, n+1):
            f = sdir + "/fuse{}/nest{}.{}_mutsig_diff_nneg_log_ps5.{}.txt".format(i, i, c, label)
            if not os.path.isfile(f):
                n_sig.append(-1)
                continue
            df = pd.read_table(f,sep='\t')
            n_sig.append(np.sum(np.array(df["q"]) < q))
        df_out['{}-mutsig-diff-ms-q{}'.format(c, q)] = np.array(n_sig)
    return df_out
        
    
def heatmap_4d(df, ctypes, var, figname=None, cut=0.3, div_by_size = False, figsize=(8, 7)):
    assert len(var) ==4
    
    nvar = [len(v)-1 for v in var]
    mat = np.zeros((nvar[0] * nvar[2], nvar[1]*nvar[3]))

    for ct in ctypes:
        col = '{}-mutsig-diff-ms-q{}'.format(ct.lower(), cut)
        mat_qs = create_data_blocks((nvar[0], nvar[1]), (nvar[2],nvar[3]), df, col, var)
        mat += mat_qs
    mat /= len(ctypes)
    
    mat_size = create_data_blocks((nvar[0],nvar[1]), (nvar[2],nvar[3]), df, 'Size', var)

    if div_by_size:
        mat/= mat_size
        
    label = 'No. systems q < {}'.format(cut)

    fig, ax = plt.subplots(1,1, figsize=figsize)
    #     plt.subplots_adjust(wspace=0.3)
    
    if div_by_size:
        mat *=100
        sns.heatmap(mat, ax =ax, xticklabels =False, yticklabels=False,  
                square=True, cmap='YlOrRd', annot=True, fmt=".3f")
    #            cbar_kws = {'ticks'})
    else:
        sns.heatmap(mat, ax =ax, xticklabels =False, yticklabels=False,  
                square=True, cmap='YlOrRd', annot=True, fmt=".1f")

    dyin, dxin = len(var[2])-1, len(var[3])-1 
    dy, dx = len(var[0])-1, len(var[1])-1 
    cb_label = col

    ax.set_yticks([dyin*(i-0.5) for i in range(1, dy+1)], minor= True)
    ax.set_yticks([dyin*i for i in range(dy)], minor= False)

    ax.set_xticks([dxin*(i-0.5) for i in range(1, dx+1)], minor= True)
    ax.set_xticks([dxin*i for i in range(dx)], minor= False)


    ax.yaxis.grid(True, which='major', linewidth=4, color='white')
    ax.xaxis.grid(True, which='major', linewidth=4, color='white')


    ax.set_xlabel(var[1][0], fontsize=16)
    ax.set_ylabel(var[0][0], fontsize=16)
    ax.set_xticklabels(var[1][1:], minor=True, fontsize=12)
    ax.set_yticklabels(var[0][1:], minor=True, fontsize=12)
    ax.set_title(r'Each {}$\times${} grid:'.format(nvar[2], nvar[3]) + 
                 '\n' +r'vertical: $\alpha$ ({})'.format(','.join(map(str,var[2])))+
                 '\n' +r'horizontal: $\beta$ ({})'.format(','.join(map(str,var[3]))))

    cbar = ax.collections[0]
    cbar.colorbar.set_label(label, fontsize=16)
    if figname!=None:
        plt.savefig(figname, bbox_inches='tight')
    
    return mat, mat_size

    