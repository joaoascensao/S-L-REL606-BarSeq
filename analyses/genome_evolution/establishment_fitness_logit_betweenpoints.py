'''

Fit model to compute correlations between fitness and log-odds of mutation establishment

'''

import pandas as pd
import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats import multitest

sns.set(style="whitegrid",font_scale=1.2)



muts = pd.read_csv('LTEE-Ecoli-data-clean.csv')
muts = muts[muts['intergenic']=='genic']

muts = muts[muts['mutation_category']!='snp_synonymous']

plucain = pd.read_csv('plucain_seq_clean2.csv')
plucain = plucain[plucain['Type']!='S-SNP'] # exclude synonymous snps


pstrains = {
	'S':'6.5S1',
	'L':'6.5L4',
}

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

timepoints = {
	'S':[20000,40000],
	'L':[20000,40000],
	'R':[5000,20000,40000],
}

SL_rel = {
	'S':['REL11830','REL11036'],
	'L':['REL11831','REL11035']
}

c=0
dsave=[]
for strain in exps:
	exps0,exps1 = exps[strain]
	for i,exp0 in enumerate(exps0):
		exp1 = exps1[i]
		fitness = pd.read_csv('../'+dirs[exp0]+'data/fitness/{}_fitness.csv'.format(exp1))
		
		df = muts.copy()
		data=[]
		gids=[]
		if strain=='R':
			for k,tt in enumerate(timepoints[strain]):
				d0 = df[df['time']==tt]
				d0 = d0[d0['mutator_status']=='non-mutator']
				fit0 = fitness.copy()
				fit0['time']=tt
				d0 = d0[['gene_ID','time']].drop_duplicates()
				d0['fix']=1
				if k!=0:
					d0.set_index('gene_ID',drop=True,inplace=True)
					d0 = d0.drop(index=gids,errors='ignore').reset_index()
				gids += list(d0['gene_ID'])
				data.append(
					fit0.merge(d0,how='outer',on=['gene_ID','time']).dropna(subset=['s'])
					)
		else:
			for k,ss in enumerate(SL_rel[strain]):
				d0 = df[df['strain']==ss]
				tt = np.mean(d0['time'])
				fit0 = fitness.copy()
				fit0['strain']=ss
				fit0['time']=tt
				d0['fix']=1
				if k!=0:
					#pass
					d0.set_index('gene_ID',drop=True,inplace=True)
					d0 = d0.drop(index=gids,errors='ignore').reset_index()
				gids += list(d0['gene_ID'])
				d1 = fit0.merge(d0,how='outer',on=['gene_ID','strain','time']).dropna(subset=['s'])#.set_index('gene_ID',drop=True)
				plucain0 = plucain[plucain['Clone']==pstrains[strain]]
				d1 = d1.drop(index=plucain0['gene_ID'],errors='ignore').reset_index()
				data.append(d1)
		
		d2 = pd.concat(data,ignore_index=True).fillna(value=0)

		dd=[]
		bb=[]
		aa=[]
		bb_std=[]
		dd_std=[]
		aa_std=[]
		pval_b=[]
		pval_d=[]
		pval_a=[]
		cts=[]
		if strain=='R':
			cts.append('0k')
		else:
			cts.append('6.5k')
		for tt in timepoints[strain]:
			dc = d2[d2['time']==tt]
			
			cond1=dc['significant']==1
			condp=dc['s']>0
			condn=dc['s']<0

			conda=dc['significant']==0
			condb=dc['s']>-0.005
			condc=dc['s']<0.005
			
			dp = dc[(cond1 & condp) | (conda & condb & condc)]
			dn = dc[(cond1 & condn) | (conda & condb & condc)]

			bnorm=np.median(np.array(dc[cond1 & condp]['s'])) #beneficial
			dnorm=np.median(np.array(dc[cond1 & condn]['s'])) #deleterious


			logreg_bene = sm.Logit(np.array(dp['fix']), sm.add_constant(np.array(dp['s'])/bnorm),missing='raise').fit()
			logreg_del = sm.Logit(np.array(dn['fix']), sm.add_constant(np.array(dn['s'])/dnorm),missing='raise').fit()
			print(logreg_bene.summary())
			print(logreg_bene.summary2())
			if (tt/1000).is_integer():
				cts.append(str(int(tt/1000))+'k')
			else:
				cts.append(str(tt/1000)+'k')
			dsave.append(pd.DataFrame({
				'side':'Beneficial',
				'beta':logreg_bene.params[1],
				'std':logreg_bene.bse[1],
				'pval':logreg_bene.pvalues[1],
				'timepoint':cts[-2]+'-'+cts[-1],
				'strain':strain,
				'experiment':exp1,
				},index=[0]))
			dsave.append(pd.DataFrame({
				'side':'Deleterious',
				'beta':logreg_del.params[1],
				'std':logreg_del.bse[1],
				'pval':logreg_del.pvalues[1],
				'timepoint':cts[-2]+'-'+cts[-1],
				'strain':strain,
				'experiment':exp1,
				},index=[0]))

def FDR_correction(pvals,alpha=0.05):
	rejected,pvals_corr,_,_=multitest.multipletests(pvals, alpha=alpha, method='fdr_bh', is_sorted=False)
	return [int(i) for i in rejected], pvals_corr


df_save=pd.concat(dsave,ignore_index=True)
s_signi,s_pval_corr=FDR_correction(df_save['pval'])
df_save['significant']=s_signi
df_save['corrected pval']=s_pval_corr

df_save.to_csv('estfitness_logit_betweenpoints.csv')

