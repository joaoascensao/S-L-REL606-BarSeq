'''
3.31.21
Joao Ascensao

MLE to infer strength of selection on knockouts, per biological replicate

python s_inference.py dir/count dir/out dir/technical_noise dir/outliers dir/meanfitness experiment,replicate mincount

'''

import pandas as pd
import scipy as sp 
import numpy as np 
import sys
import time
import s_likelihoods as lh
from multiprocessing import Pool, TimeoutError
from functools import partial
from itertools import compress
import os
pd.set_option('mode.chained_assignment', None)

if len(sys.argv)!=8:
	raise Exception('Specify all files')


min_bc=4
mincount=np.int(sys.argv[7])
minabscount=3
RD_threshold=6

countdir=sys.argv[1]
outdir=sys.argv[2]
betadir=sys.argv[3]
outlierdir=sys.argv[4]
meanfitdir=sys.argv[5]

exprep=sys.argv[6].split(',')
experiment=exprep[0]
replicate=exprep[1][0]
counts=pd.read_csv(countdir + '/{}{}_counts.csv'.format(experiment,replicate)).set_index('barcode',drop=True)
total_counts=pd.read_csv(countdir + '/{}{}_Rtot.csv'.format(experiment,replicate)).drop(columns='Unnamed: 0')
betas=pd.read_csv(betadir+ '/{}{}_betas.csv'.format(experiment,replicate))
xbar=pd.read_csv(meanfitdir+ '/{}{}_xbar.csv'.format(experiment,replicate))
outliers=pd.read_csv(outlierdir+'/{}{}_outliers.csv'.format(experiment,replicate))
cond1=outliers['RD']>RD_threshold
cond2=outliers['RD']==-2
outliers=outliers[cond1 | cond2]

counts.drop(index=outliers['barcode'],inplace=True)


#sampled_users = np.random.choice(counts0["gene_ID"].unique(), 30)
#counts = counts0.groupby('gene_ID').filter(lambda x: x["gene_ID"].values[0] in sampled_users)

print('Estimating fitness for {}{}'.format(experiment,replicate))

tv = list(total_counts.columns)
d0 = tv[0]
dl = tv[-1]

counts=counts[counts[d0]>minabscount]

#counts = counts[counts['gene_symbol']=='lacY']

shareddata=lh.getshareddata(betas,xbar)

start=time.time()

def get_s_est(ggene):
	g,gene = ggene
	if len(gene)<min_bc or pd.isnull(g[0]):
		return None

	gene = gene[tv]

	#merge initially low count barcodes into next lowest bc
	while np.min(gene[d0])<mincount and len(gene)>1:
		minloc=gene[d0].idxmin()
		minseries=gene.loc[minloc]
		gene.drop(inplace=True,index=minloc)
		minloc2=gene[d0].idxmin()
		gene.loc[minloc2]+=minseries


	if len(gene)<1:
		return None
	
	data_list=[]
	for _,row in gene.iterrows():
		d00=lh.getgenedata(row,betas)
		if len(d00['r'])>2:
			data_list.append(d00)
	if len(data_list)==0:
		return None
	s,s_std,s_pval,success=lh.s_maxlikelihood(data_list, shareddata)
	#print(g[2],s,s_std)
	if success:
		return pd.DataFrame({
				'gene_ID':g[0],
				'type':g[1],
				'gene_symbol':g[2],
				'description':g[3],
				's':s,
				's std':s_std,
				'pval':s_pval,
			},index=[0])
	else:
		return None

batches = list(counts.groupby(by=['gene_ID','type','gene_symbol','description']))
with Pool(processes=4) as pool:
	s_data = pool.map(get_s_est, batches)
s_data.append(None)
cond = pd.notnull(s_data)
s_data = list(compress(s_data, list(cond)))

df_save=pd.concat(s_data,ignore_index=True)
s_signi,s_pval_corr=lh.FDR_correction(df_save['pval'])
df_save['s significant']=s_signi
df_save['s corrected pval']=s_pval_corr
df_save.to_csv(outdir+'/{}{}_fitness.csv'.format(experiment,replicate))
print('Wrote '+outdir+'/{}{}_fitness.csv'.format(experiment,replicate))
print((time.time()-start)/3600)
