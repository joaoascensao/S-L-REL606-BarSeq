'''

Fit model to compute correlations between fitness and log-changes in gene expression

'''

import pandas as pd
import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats import multitest
import statsmodels.api as sm
from statsmodels.stats import multitest



fit3 = pd.read_csv('fit3.csv')

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
		['R','I','O','T','U'])
}

contrasts={
		0:"X6kS.rel606",
		1:"X6kL.rel606",
		2:"X17kL.X6kL",
		3:"X17kS.X6kS",
		4:"X40kL.X6kL",
		5:"X40kS.X6kS",
		6:"X6kS.X6kL",
		7:"X17kS.rel606",
		8:"X17kL.rel606",
		9:"X40kS.rel606",
		10:"X40kL.rel606",
		}

# log-fold changes relevant for each genetic background
contrasts_per_exp = {
	'S':[0,3,5],
	'L':[1,2,4],
	'R':[0,1,7,8,9,10],
}

pval_name = 'p.value.'
coeff_name = 'coefficients.'
stdev_name = 'stdev.unscaled.'

blue='#69a2ff'
red='#f52257'

c=0
data=[]
for strain in exps:
	exps0,exps1 = exps[strain]
	for i,exp0 in enumerate(exps0):
		exp1 = exps1[i]
		fitness = pd.read_csv('../'+dirs[exp0]+'data/fitness/{}_fitness.csv'.format(exp1))
		df = fitness.merge(fit3,how='inner',on='gene_ID')
		

		for constrast_num in contrasts_per_exp[strain]:
			ct = contrasts[constrast_num]
			cond1=df['significant']==1
			condp=df['s']>0
			condn=df['s']<0

			conda=df['significant']==0
			condb=df['s']>-0.005
			condc=df['s']<0.005


			neutrals=np.array(df[conda & condb & condc][coeff_name+ct]) #neutral
			beneficial=np.array(df[cond1 & condp][coeff_name+ct]) #beneficial
			deleterious=np.array(df[cond1 & condn][coeff_name+ct]) #deleterious

			beneficial_s=np.array(df[cond1 & condp]['s']) #beneficial
			deleterious_s=np.array(df[cond1 & condn]['s']) #deleterious
			neutral_s=np.abs(np.array(df[conda & condb & condc]['s'])) #deleterious

			neutrals_w=np.array(df[conda & condb & condc][stdev_name+ct])**-2 #neutral
			beneficial_w=np.array(df[cond1 & condp][stdev_name+ct])**-2 #beneficial
			deleterious_w=np.array(df[cond1 & condn][stdev_name+ct])**-2 #deleterious

			bnorm=np.median(beneficial_s)
			dnorm=np.median(deleterious_s)

			
			wls_model_b = sm.WLS(np.concatenate((neutrals,beneficial),axis=None),
				sm.add_constant(np.concatenate((neutral_s,beneficial_s),axis=None)/bnorm),
				weights=np.concatenate((neutrals_w,beneficial_w),axis=None))
			wls_model_d = sm.WLS(np.concatenate((neutrals,deleterious),axis=None),
				sm.add_constant(np.concatenate((neutral_s,deleterious_s),axis=None)/dnorm),
				weights=np.concatenate((neutrals_w,deleterious_w),axis=None))
			results_b = wls_model_b.fit()
			results_d = wls_model_d.fit()


			data.append(pd.DataFrame({
				'side':'Beneficial',
				'beta':results_b.params[1],
				'std':results_b.HC0_se[1],
				'pval':results_b.pvalues[1],
				'contrast':ct.replace('X','').replace('rel606','R'),
				'strain':strain,
				'experiment':exp1,
				},index=[0]))
			data.append(pd.DataFrame({
				'side':'Deleterious',
				'beta':results_d.params[1],
				'std':results_d.HC0_se[1],
				'pval':results_d.pvalues[1],
				'contrast':ct.replace('X','').replace('rel606','R'),
				'strain':strain,
				'experiment':exp1,
				},index=[0]))



def FDR_correction(pvals,alpha=0.05):
	rejected,pvals_corr,_,_=multitest.multipletests(pvals, alpha=alpha, method='fdr_bh', is_sorted=False)
	return [int(i) for i in rejected], pvals_corr


df_save=pd.concat(data,ignore_index=True)
s_signi,s_pval_corr=FDR_correction(df_save['pval'])
df_save['significant']=s_signi
df_save['corrected pval']=s_pval_corr

df_save.to_csv('foldchange_fitness_lm_byfitness.csv')

