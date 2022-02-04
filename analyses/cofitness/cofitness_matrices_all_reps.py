'''

Compute cofitness matrices for each genetic background, only using either replicate "1" or "2" for each experiment

'''

import pandas as pd
import scipy as sp 
import numpy as np 
import sys
import matplotlib.pyplot as plt
import scipy.interpolate
import seaborn as sns
import time
from multiprocessing import Pool, TimeoutError
from statsmodels.stats import multitest
from functools import partial
from itertools import compress

strain=sys.argv[1]
rep=sys.argv[2]
print('Running {}{}'.format(strain,rep))


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

exps0 = {
	'S':(['E_SLR', 'E_CD', 'E_F', 'E_GHI', 'E_MNO', 'E_PQT' 'E_VW', 'E_XY'],
		['S','C','FS','G','M','P','WS','X']),
	'L':(['E_SLR', 'E_CD', 'E_F', 'E_GHI', 'E_MNO', 'E_PQT', 'E_VW','E_VW', 'E_XY'],
		['L''D','FL','H','N','Q','WL','V','Y',]),
	'R':(['E_SLR', 'E_GHI', 'E_MNO', 'E_PQT', 'E_U'],
		['R','I','O','T','U'])}

exps = exps0[strain][0]
exps1 = exps0[strain][1]

for i,exp in enumerate(exps):
	iname=exps1[i]
	dfi = pd.read_csv(dirs[exp]+'data/fitness/{}{}_fitness.csv'.format(iname,rep))
	dfi['isig'] = 1 - dfi['s significant']
	dfi=dfi[['gene_ID','s','s std','isig']].rename(columns={'s':'s'+iname, 's std':'s std'+iname, 'isig':'isig'+iname})
	if i==0:
		df = dfi
	else:
		df = df.merge(dfi,how='outer',on='gene_ID')

# exlude genes that are not significantly non-neutral in at least one experiment
for i,iname in enumerate(exps1):
	df['isig' + iname] = df['isig' + iname].fillna(1)
	cond0 = df['isig' + iname] == 0
	if i==0:
		cond = cond0
	else:
		cond = cond | cond0

alls = ['s'+exp for exp in exps1]
allstds = ['s std'+exp for exp in exps1]
df = df[cond]
df['notnull'] = df[alls].apply(lambda row: np.sum(pd.notnull(row)), axis=1)
df = df[df['notnull']>=4] # exclude genes that were measured in <4 experiments


gids = list(df['gene_ID'])
gid_pairs = []
for i1,gid1 in enumerate(gids):
	for i2,gid2 in enumerate(gids[i1+1:]):
		gid_pairs.append((gid1,gid2))



def wcov(x,y,mx,my,w):
	return np.sum(w * (x - mx) * (y - my))

def wperson(x,y,w):
	mx, my = (np.sum(i * w) / np.sum(w) for i in [x, y])
	return wcov(x,y,mx,my,w)/np.sqrt(wcov(x,x,mx,mx,w) * wcov(y,y,my,my,w))

def unison_shuffled_copies(a, b):
    assert len(a) == len(b)
    p = np.random.permutation(len(a))
    return a[p], b[p]

def permutation_test_single(s10,std10,s20,std20,n):
	s1 = np.random.normal(loc=s10,scale=std10)
	s2 = np.random.normal(loc=s20,scale=std20)
	s1,std1=unison_shuffled_copies(s1,std10)
	s2,std2=unison_shuffled_copies(s2,std20)
	w = 1/(std1**2 + std2**2)
	return wperson(s1,s2,w)

def permutation_test_ntimes(s1,std1,s2,std2,n,corr0,nperm=300):
	corrs=[]
	for i in range(nperm):
		corrs.append(permutation_test_single(s1,std1,s2,std2,n))
	corrs = np.sort(corrs)
	ps = np.linspace(1/nperm,1,nperm)
	if corr0>np.max(corrs) or corr0<np.min(corrs):
		return 0
	else:
		ii = np.asarray(corrs>=corr0).nonzero()
		x  = corr0
		x1 = corrs[ii[0][0]]
		x2 = corrs[ii[0][0]-1]
		y1 = ps[ii[0][0]]
		y2 = ps[ii[0][0]-1]
		p = (y1*(x2 - x) + y2*(x - x1))/(x2 - x1)
		if corr0>0:
			return 1-p
		else:
			return p

start=time.time()

def getcofitness(gid_pairs):
	gid1, gid2 = gid_pairs
	gene1 = df[df['gene_ID']==gid1]
	gene2 = df[df['gene_ID']==gid2]
	s1 = np.array(gene1[alls])[0]
	std1 = np.array(gene1[allstds])[0]
	s2 = np.array(gene2[alls])[0]
	std2 = np.array(gene2[allstds])[0]
	cond1 = pd.notnull(s1)
	cond2 = pd.notnull(s2)
	s1 = s1[cond1 & cond2]
	std1 = std1[cond1 & cond2]
	s2 = s2[cond1 & cond2]
	std2 = std2[cond1 & cond2]
	n=len(s1)
	w = 1/(std1**2 + std2**2)
	try:
		corr = wperson(s1,s2,w)

		if np.abs(corr)<0.5:
			pval=1
		else:
			pval=permutation_test_ntimes(s1,std1,s2,std2,n,corr)

		return pd.DataFrame({
				'gene_ID 1':gid1,
				'gene_ID 2':gid2,
				'corr':corr,
				'pval':pval,
			},index=[0])
	except:
		return None


with Pool(processes=4) as pool:
	cofitness_data = pool.map(getcofitness, gid_pairs)
cofitness_data.append(None)
cond = pd.notnull(cofitness_data)
cofitness_data = list(compress(cofitness_data, list(cond)))

print(time.time()-start)

def FDR_correction(pvals,alpha=0.05):
	rejected,pvals_corr,_,_=multitest.multipletests(pvals, alpha=alpha, method='fdr_bh', is_sorted=False)
	return [int(i) for i in rejected], pvals_corr

df_save=pd.concat(cofitness_data,ignore_index=True)
s_signi,s_pval_corr=FDR_correction(df_save['pval'])
df_save['significant']=s_signi
df_save['corrected pval']=s_pval_corr
df_save.to_csv('{}{}_cofitness.csv'.format(strain,rep))