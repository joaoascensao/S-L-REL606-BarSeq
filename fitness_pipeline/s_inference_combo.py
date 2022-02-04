'''
MLE to infer strength of selection on knockouts, combining all replicates

python s_inference.py dir/count dir/out dir/technical_noise dir/outliers dir/meanfitness mincount

'''

import pandas as pd
import scipy as sp 
import numpy as np 
import sys
import time
import s_likelihoods as lh
import s_likelihoods_combo as lhc
from multiprocessing import Pool, TimeoutError
from itertools import compress
pd.set_option('mode.chained_assignment', None)

if len(sys.argv)!=9:
	print(sys.argv)
	raise Exception('Specify all files')


min_bc=4
mincount=np.int(sys.argv[8])
minabscount=3
RD_threshold=6

countdir=sys.argv[1]
outdir=sys.argv[2]
betadir=sys.argv[3]
outlierdir=sys.argv[4]
fitdir=sys.argv[5]
meanfitdir=sys.argv[6]
explist=pd.read_csv(sys.argv[7])


for experiment,group in explist.groupby(by='Experiment'):
	
	print('Estimating fitness for {}'.format(experiment))
	#s_data=[]
	counts={}
	shareddata={}
	itgdata={}
	tvs={}
	c=0
	for replicate in group['Replicate']:
		fit0 = pd.read_csv(fitdir+'/{}{}_fitness.csv'.format(experiment,replicate))
		if c==0:
			fit = fit0
		else:
			fit = fit0.merge(fit,how='outer',on=['gene_ID','type','gene_symbol','description'])
		c+=1
		counts0=pd.read_csv(countdir + '/{}{}_counts.csv'.format(experiment,replicate)).set_index('barcode',drop=True)
		total_counts=pd.read_csv(countdir + '/{}{}_Rtot.csv'.format(experiment,replicate)).drop(columns='Unnamed: 0')
		betas=pd.read_csv(betadir+ '/{}{}_betas.csv'.format(experiment,replicate))
		xbar=pd.read_csv(meanfitdir+ '/{}{}_xbar.csv'.format(experiment,replicate))
		outliers=pd.read_csv(outlierdir+'/{}{}_outliers.csv'.format(experiment,replicate))
		cond1=outliers['RD']>RD_threshold
		cond2=outliers['RD']==-2
		outliers=outliers[cond1 | cond2]
		counts0.drop(index=outliers['barcode'],inplace=True)

		#sampled_users = np.random.choice(counts0["gene_ID"].unique(), 30)
		#counts0 = counts0.groupby('gene_ID').filter(lambda x: x["gene_ID"].values[0] in sampled_users)
		
		shareddata[replicate]=lh.getshareddata(betas,xbar)
		tv = list(total_counts.columns)
		tvs[replicate]=tv
		d0 = tv[0]
		dl = tv[-1]
		counts0=counts0[counts0[d0]>minabscount]
		counts[replicate]=counts0

	def get_s_est(ggene):
		g,_ = ggene
		data_list={}
		totsize=0
		nrep=0
		for replicate in group['Replicate']:
			data_list0=[]
			counts0=counts[replicate]
			gene=counts0[counts0['gene_ID']==g[0]]

			totsize+=len(gene)
			nrep+=1
			
			gene = gene[tvs[replicate]]
			d0=tvs[replicate][0]
			
			if len(gene)<1:
				return None
			while np.min(gene[d0])<mincount and len(gene)>1:
				minloc=gene[d0].idxmin()
				minseries=gene.loc[minloc]
				gene.drop(inplace=True,index=minloc)
				minloc2=gene[d0].idxmin()
				gene.loc[minloc2]+=minseries
			
			if len(gene)<1:
				return None
			
			for _,row in gene.iterrows():
				d00=lh.getgenedata(row,betas)
				if len(d00['r'])>2:
					data_list0.append(d00)
			if len(data_list0)>0:
				data_list[replicate]=data_list0
			
		if totsize<min_bc*nrep:
			return None
		if np.sum([len(data_list[rep]) for rep in data_list])==0:
			return None

		s,s_std,s_pval=lhc.s_maxlikelihood_combo(data_list, shareddata)
		
		return pd.DataFrame({
				'gene_ID':g[0],
				'type':g[1],
				'gene_symbol':g[2],
				'description':g[3],
				's':s,
				's std':s_std,
				'pval':s_pval,
			},index=[0])
	
	fit.set_index(['gene_ID','type','gene_symbol','description'],inplace=True)
	
	batches = list(fit.iterrows())
	
	with Pool(processes=4) as pool:
		s_data = pool.map(get_s_est, batches)
	#s_data=list(s_data)
	s_data.append(None)
	cond = pd.notnull(s_data)
	s_data = list(compress(s_data, list(cond)))
	df_save=pd.concat(s_data,ignore_index=True)
	signi,pval_corr=lh.FDR_correction(df_save['pval'])
	df_save['significant']=signi
	df_save['corrected pval']=pval_corr
	df_save.to_csv(fitdir+'/{}_fitness.csv'.format(experiment))
	print('Wrote '+fitdir+'/{}_fitness.csv'.format(experiment))