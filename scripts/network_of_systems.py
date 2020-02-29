import argparse
import pandas as pd
import sys
sys.path.append('/cellar/users/f6zheng/tools/lib/python2.7/site-packages') # TODO: remove this
from ddot import *
import seaborn as sns
import networkx as nx
from myutils import * # TODO: remove this


def selected_system_diagram(ont, df, key, used_cols, genes=None):
    '''
    '''
    # ont = Ontology.from_table(ontf, clixo_format=True, is_mapping=lambda x: x[2] == True)
    # if len(ont.get_roots()) > 1:
    #     ont.add_root('ROOT', inplace=True)
    # ont.propagate('forward', inplace=True)
    chosen_terms = df[key].tolist()

    ont_sel = ont.delete(to_keep=chosen_terms +['ROOT']+ ont.genes)

    G = ont_sel.to_NdexGraph(layout='bubble')
    edges_to_remove = []
    for e in G.edges():
        if len(list(nx.all_simple_paths(G, e[0], e[1]))) > 1:
            edges_to_remove.append((e[0], e[1]))
    G.remove_edges_from(edges_to_remove)

    mapping = {}
    for n in G.nodes(data=True):
        if (n[1]['NodeType'] == 'Gene') and ((n[1]['Label'] in genes)==False):
            G.remove_node(n[0])
        mapping[n[0]] = n[1]['Label']
    H = nx.relabel_nodes(G, mapping)


    # add node attributes
    df2 = df.set_index(key)
    for c in used_cols:
        if c != key:
            nx.set_node_attributes(H, c, df2[c].to_dict())
    return H


par = argparse.ArgumentParser()
par.add_argument('--input', required=True, help='dataframe, same format to the output of parse.py') # TODO: combine multiple files
par.add_argument('--node_attr', help='a dataframe with all the other information of systems')
par.add_argument('--join', default='System_name', help='use this string to join the two dataframe')
par.add_argument('--exclude', default=[], nargs='*', help='if specified, excluding certain systems from netwokrs')
par.add_argument('--ont', required=True, help='the hierarchy 3-col file')
par.add_argument('--genes', default=[], help='a string of gene symbols separated by comma')
par.add_argument('--gene_attr', help='a dataframe with gene attributes')
par.add_argument('--name', help='name for the NDEx network')
par.add_argument('--out', required=True, help='write an output dataframe')
args = par.parse_args()

if __name__ == "__main__":
    # TODO: make NDEx login optional
    ont = Ontology.from_table(args.ont,
                              is_mapping=lambda x:x[2]=='gene', clixo_format=True)
    if not 'ROOT' in ont.terms:
        ont.add_root('ROOT', inplace=True)
    ont.propagate('forward', inplace=True)

    df_use = pd.read_table(args.input, sep='\t')
    assert args.join in df_use.columns.tolist()

    if args.node_attr != None:
        df_node = pd.read_table(args.node_attr, sep='\t')
    # remove the redundant columns
        df_use = df_use[[c for c in df_use.columns.tolist() if (c==args.join) or (not c in df_node.columns.tolist())]]
        assert args.join in df_node.columns.tolist()

    terms = df_use[args.join].tolist()

    # add auxiliary "integrator" nodes
    ont_conn = ont.connected(ont.terms)
    records = []
    ind_sel_terms = np.array([ont.terms.index(t) for t in terms if t in ont.terms])
    for ti in range(len(ont.terms)):
        x = np.sum(ont_conn[ind_sel_terms, ti])
        if (x> 1) and (not ont.terms[ti] in terms):
            records.append([ont.terms[ti], x])

    df_ancestor = pd.DataFrame(records)
    df_ancestor.sort_values(by=1, ascending=False, inplace=True)
    df_ancestor = df_ancestor.loc[(df_ancestor[0] != 'ROOT') & (df_ancestor[0].isin(args.exclude) == False), :]
    df_ancestor = df_ancestor[[0]]
    df_ancestor.rename(columns = {0:'System_name'}, inplace=True)

    # TODO: what is this part
    for c in df_use.columns:
        if not c in df_use.columns:
            if (df_use[c].dtype == 'int64') or (df_use[c].dtype == 'float64'):
                df_use[c] = -1
            else:
                df_use[c] = ''

    df_use['Integrator'] = False
    df_ancestor['Integrator'] = True
    df_use = pd.concat([df_use, df_ancestor])
    df_use.reset_index(inplace=True, drop=True)

    ### aggregate cancer types from descendents
    dict_ctypes = {}
    dict_ctypes_collect = {}

    ind_used_terms = np.array([ont.terms.index(t) for t in df_use['System_name'].tolist()])
    (indx, indy) = np.where(ont_conn)

    for i, row in df_use.iterrows():
        sys_name, ctypes = row['System_name'], row['Activated in cancer']
        if row['Integrator']:
            continue
        dict_ctypes[sys_name] = ctypes.split(' ')

    for i in range(len(indx)):
        if (not indx[i] in ind_sel_terms) or (not indy[i] in ind_used_terms):
            continue
        t_des, t_ans = ont.terms[indx[i]], ont.terms[indy[i]]
        if not t_ans in dict_ctypes_collect:
            dict_ctypes_collect[t_ans] = []
        dict_ctypes_collect[t_ans].extend(dict_ctypes[t_des])

    df_use['Activate in cancer collected'] = np.array([' '.join(sorted(list(set(dict_ctypes_collect[row['System_name']]))))  for i, row in df_use.iterrows()])
    df_use['N cancer collected'] = np.array([len(set(dict_ctypes_collect[row['System_name']])) for i, row in df_use.iterrows()])

    if args.node_attr != None:
        df_use = df_use.merge(df_node, how='left', left_on =args.join,right_on=args.join)

    # add term size
    term_sizes = {i: ont.term_sizes[ont.terms.index(row['System_name'])] for i, row in df_use.iterrows()}
    df_use['Size'] = pd.Series(term_sizes)
    term_genes = {i: ' '.join(sorted([ont.genes[g] for g in ont.term_2_gene[row['System_name']]])) for i, row in df_use.iterrows()}
    df_use['Genes'] = pd.Series(term_genes)

    # TODO: remove unneeded columns

    # add enhanced graph (for RBVI)
    df_combined_eh = df_use.copy()
    cols_q = [c for c in df_use.columns.tolist() if c.endswith('-q')]

    df_combined_eh['logSize'] = np.log2(df_combined_eh['Size'])

    df_combined_eh = df_combined_eh.loc[df_combined_eh[args.join].isin(args.exclude) ==False, :]
    df_combined_eh = df_combined_eh.round(3)

    if args.gene_attr != None:
        df_gene_attr = pd.read_table(args.gene_attr, index_col=0)
        df_gene_attr.reset_index()
        df_gene_attr.rename(columns={'Unnamed: 0':'System_name'}, inplace=True)
        df_gene_attr = df_gene_attr[[c for c in df_gene_attr.columns.tolist() if c in df_combined_eh.columns.tolist()]]
        df_combined_eh = pd.concat([df_combined_eh, df_gene_attr])

    # finally, upload to NDEx
    df_combined_eh = df_combined_eh.loc[df_combined_eh['System_name'].isin(ont.terms + ont.genes), :]
    # create nest IDs for systems
    df_combined_eh.sort_values(by ='Size', ascending=False, inplace=True)
    df_combined_eh.reset_index(inplace=True, drop=True)
    df_combined_eh['Nest ID'] = np.array(['NEST:{}'.format(i+1) for i in range(df_combined_eh.shape[0])])
    df_combined_eh.to_csv(args.out, sep='\t')
    
    # also create a new .ont object
    if args.name !=None:
        G =selected_system_diagram(ont, df_combined_eh, 'System_name', df_combined_eh.columns.tolist(), args.genes)

        Gnx = nx_to_NdexGraph(G)
        Gnx.set_name(args.name)

        cx = Gnx.to_cx()
        counter = 0
        for entry in cx:
            if 'nodeAttributes' in entry:
                for subentry in entry['nodeAttributes']:
                    if 'd' in subentry and subentry['d'] == 'unknown':
                        subentry['d'] = 'double'
                        counter = counter + 1

        G = NdexGraph(cx=cx)
        serv, usr, passwd = ndex_login('http://www.ndexbio.org')
        G.upload_to(serv, usr, passwd)
