'''
5.22.20
Joao Ascensao

good one:
exp='E7'
e1='JL'
gene='yhfM'
rep=1

'''

import pandas as pd
import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats import multitest
import statsmodels.api as sm
from statsmodels.stats import multitest




dirs={
	'E_CD':'../../E_CD/',
	'E_F':'../../E_F/',
	'E_GHI':'../../E_GHI/',
	'E_MNO':'../../E_MNO/',
	'E_PQT':'../../E_PQT/',
	'E_SLR':'../../E_SLR/',
	'E_U':'../../E_U/',
	'E_VW':'../../E_VW/',
	'E_XY':'../../E_XY/',
}

exps = {
	'S':(['E_SLR', 'E_CD', 'E_F', 'E_GHI', 'E_MNO', 'E_PQT' 'E_VW', 'E_XY'],
		['S','C','FS','G','M','P','WS','X']),
	'L':(['E_SLR', 'E_CD', 'E_F', 'E_GHI', 'E_MNO', 'E_PQT', 'E_VW','E_VW', 'E_XY'],
		['L''D','FL','H','N','Q','WL','V','Y',]),
	'R':(['E_SLR', 'E_GHI', 'E_MNO', 'E_PQT', 'E_U'],
		['R','I','O','T','U'])}

expnames = {
	'S':'Mono',
	'L':'Mono',
	'R':'Mono 1',
	'C':'L in maj',
	'D':'S in maj 1',
	'FS':'Eco Eq 1',
	'FL':'Eco Eq 1',
	'G':'Ac exp',
	'H':'Ac exp',
	'I':'Ac exp',
	'M':'1:10 dil',
	'N':'1:10 dil',
	'O':'1:10 dil',
	'P':'Glu exp',
	'Q':'Glu exp',
	'T':'Glu exp',
	'U':'Mono 2',
	'WS':'S in maj 2',
	'WL':'S in maj 2',
	'V':'S in maj 3',
	'X':'Eco Eq 2',
	'Y':'Eco Eq 2',
}

df = pd.read_csv('../ecoli_net_B.csv').drop(columns='Unnamed: 0')
#df = pd.read_csv('../../regulonDB/genetic_network.csv').drop(columns='Unnamed: 0').rename(columns={'tf gene_ID':'gene_ID 1','reg gene_ID':'gene_ID 2'})
genes = pd.read_csv('../../../../rel606_genes_wrangled.csv').rename(columns={'gene_ID':'gene_ID 1', 'gene_symbol':'gene_symbol 1'})

genes['gene_ID 2'] = genes['gene_ID 1']
genes['gene_symbol 2'] = genes['gene_symbol 1']
genes = genes[['gene_ID 1','gene_ID 2','gene_symbol 1','gene_symbol 2']]
#print(genes)

# export statistics of cofitness network
# distribution of connections per gene

dfc = pd.concat((df[['gene_ID 1','ll']].rename(columns={'gene_ID 1':'gene_ID'}),df[['gene_ID 2','ll']].rename(columns={'gene_ID 2':'gene_ID'})))
#dfc = pd.concat((df[['gene_ID 1']].rename(columns={'gene_ID 1':'gene_ID'}),df[['gene_ID 2']].rename(columns={'gene_ID 2':'gene_ID'})))

edgenums=[]
gids=[]
for i,group in dfc.groupby(by='gene_ID'):
	edgenums.append(len(group)) #np.sum(group['ll'])
	gids.append(i)

edgedf = pd.DataFrame({'edgenum':edgenums,'gene_ID':gids})

blue='#69a2ff'
red='#f52257'

c=0
data=[]
for strain in exps:
	exps0,exps1 = exps[strain]
	for i,exp0 in enumerate(exps0):
		exp1 = exps1[i]
		fitness = pd.read_csv('../../{}/data/fitness/{}_fitness.csv'.format(dirs[exp0],exp1))
		df = fitness.merge(edgedf,how='inner',on='gene_ID')

		cond1=df['significant']==1
		condp=df['s']>0
		condn=df['s']<0

		conda=df['significant']==0
		condb=df['s']>-0.005
		condc=df['s']<0.005


		neutrals=np.array(df[conda & condb & condc]['edgenum']) #neutral
		beneficial=np.array(df[cond1 & condp]['edgenum']) #beneficial
		deleterious=np.array(df[cond1 & condn]['edgenum']) #deleterious

		beneficial_s=np.array(df[cond1 & condp]['s']) #beneficial
		deleterious_s=np.array(df[cond1 & condn]['s']) #deleterious
		bnorm = np.median(beneficial_s)
		dnorm = np.median(deleterious_s)

		
		
		wls_model_b = sm.OLS(np.concatenate((neutrals,beneficial),axis=None),
			sm.add_constant(np.concatenate((np.zeros_like(neutrals),beneficial_s)/bnorm,axis=None)))
		wls_model_d = sm.OLS(np.concatenate((neutrals,deleterious),axis=None),
			sm.add_constant(np.concatenate((np.zeros_like(neutrals),deleterious_s)/dnorm,axis=None)))
		results_b = wls_model_b.fit()
		results_d = wls_model_d.fit()



		data.append(pd.DataFrame({
			'side':'Beneficial',
			'beta':results_b.params[1],
			'std':results_b.HC0_se[1],
			'pval':results_b.pvalues[1],
			'strain':strain,
			'experiment':exp1,
			},index=[0]))
		data.append(pd.DataFrame({
			'side':'Deleterious',
			'beta':results_d.params[1],
			'std':results_d.HC0_se[1],
			'pval':results_d.pvalues[1],
			'strain':strain,
			'experiment':exp1,
			},index=[0]))




df_save=pd.concat(data,ignore_index=True)

df_save.to_csv('lm_degrees.csv')

