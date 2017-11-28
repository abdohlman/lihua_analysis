import pandas as pd
import numpy as np

conditions = ['no1','no2','no3','no4','allsamples']
methods = ['DESeq2','DESeq2Block','EdgeR','EdgeRBlock']

indices = [ cond+'-'+method for cond in conditions for method in methods ]

table = pd.DataFrame(columns=['pval_up','pval_dn','FDR_up','FDR_dn'])

for cond in conditions:
	for method in methods:

		filename = '.'.join(['DBA_atac',cond,method,'total.csv'])

		df = pd.read_table(filename, sep=',')
		# print len(df[(df['p-value']< 0.05) & (df['Fold'] < 0.0)])
		# table.add(filename, axis='index')
		table.loc[cond+'-'+method, 'pval_up'] = len(df[(df['p-value']< 0.05) & (df['Fold'] > 1.0)])
		table.loc[cond+'-'+method, 'pval_dn'] = len(df[(df['p-value']< 0.05) & (df['Fold'] < -1.0)])
		table.loc[cond+'-'+method, 'FDR_up'] = len(df[(df['FDR']< 0.05) & (df['Fold'] > 1.0)])
		table.loc[cond+'-'+method, 'FDR_dn'] = len(df[(df['FDR']< 0.05) & (df['Fold'] < -1.0)])

		deps = df[ np.abs(df['Fold']) > 1.0 ]
		deps_bed = deps.iloc[:, [:2]]
		deps_bed.to_csv('bed/'+filename.replace('.csv','.bed'), sep='\t')
#		df = pd.read_table(filename)

table.to_csv('DEP_counts.txt',sep='\t')
print table