import argparse
import numpy as np
import pandas as pd

par = argparse.ArgumentParser()
par.add_argument('--coef', required=True)
par.add_argument('--size', required=True)
par.add_argument('--ngenes', type=int, required=True)
par.add_argument('--bins', nargs = '+', type=int)
par.add_argument('--out', required=True)
args = par.parse_args()

# bins = []
# for i in range(len(args.bins)-1):
#     bins.append((args.bins[i], args.bins[i+1]-1))
# bins.append((args.bins[-1], 20000))

df_coef = pd.read_table(args.coef, sep="\t")
## TODO: temporary; should solve this later
df_coef = df_coef[[c for c in df_coef.columns.tolist() if not c.endswith('.1')]]
df_coef.reset_index(inplace=True)
coef_melt = pd.melt(df_coef, id_vars=['index'])

df_size = pd.read_table(args.size, sep="\t", index_col=0, header=None)
df_size.sort_index(inplace=True)
#
df_size.reset_index(inplace=True)
df_size.rename(columns = {0: 'term_name', 1:"Size"}, inplace=True)
df_size.index = ['V' + str(args.ngenes + i+1) for i in df_size.index]

coef_melt = coef_melt.merge(df_size, how='left', left_on = 'index', right_index=True)
coef_melt.fillna(1, inplace=True)
coef_melt.loc[coef_melt['index'] =='(Intercept)', 'Size'] = 0
coef_melt.loc[coef_melt['value'] <0, 'value'] = 0
coef_melt = coef_melt.loc[coef_melt['value']>0, :]

lambdas = coef_melt['variable'].unique()
lambdas.sort()
lambdas = lambdas.tolist()
lambdas = sorted(lambdas, key = lambda x:float(x))

mat = np.zeros((len(lambdas), len(args.bins)))

for i, row in coef_melt.iterrows():
    if row['value'] == 0:
        continue
    l = row['variable']
    l_ind = lambdas.index(l)
    b_ind = np.digitize(row['Size'], args.bins)-1
    mat[l_ind, b_ind] += row['value']
    # if i % 10000 == 0:
    #     print(i)


with open(args.out, 'w') as fh:
    fh.write('{}\t{}\t{}\n'.format('Lambda', 'Size', 'Coefficient'))
    for i in range(len(args.bins)):
        for j in range(len(lambdas)):
            fh.write(lambdas[j] + '\t' + str(args.bins[i]) + '\t{:.6f}'.format(mat[j,i]) +'\n')