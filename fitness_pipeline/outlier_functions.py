'''
Functions for outlier detection

'''
import pandas as pd
import scipy as sp 
import numpy as np
import random
from itertools import compress
from multiprocessing import Pool, TimeoutError
from functools import partial
import s_likelihoods as lh
from scipy import interpolate
import os

# trimmed weighted mean
def trimmed_wmean(x0,std,trim=0.3):
	x = np.sort(x0)
	n = np.int(np.round(len(x)*trim))
	inv_var = std[n:len(x)-n]**-2
	return np.sum(x[n:len(x)-n]*inv_var)/np.sum(inv_var)

# weighted median
def weighted_median(x,std):
	xa=np.argsort(x)
	w = np.array(std)**-2
	w = w[xa]/np.sum(w)
	if len(x)==2:
		return x[0]*w[0] + x[1]*w[1]
	wc=np.insert(np.cumsum(w),0,0)
	cdf = interpolate.interp1d([(wc[i]+wc[i+1])/2 for i in range(len(wc)-1)], x[xa],kind='linear')
	return cdf(0.5)
	
# compute log-likelihood ratio deviation of each barcode from rest of gene using a jacknife approach
def gene_boot_sampling(i,data_list,shareddata):
	dsample = random.sample(data_list,int(np.ceil(len(data_list)/2)))
	sl=np.array([dd['s'] for dd in dsample])
	stds=np.array([dd['std'] for dd in dsample])
	
	if np.all(np.isfinite(stds)):
		if len(dsample)<10:
			sm = weighted_median(sl,stds)
		else:
			sm=trimmed_wmean(sl,stds)
	else:
		sm = np.median(sl)
	
	LRs=[]
	for data in data_list:
		nll_null = lh.add_likelihoods(sm,[data],shareddata)
		LR = data['nll'] - nll_null
		LRs.append(LR)
	LRs=np.array(LRs)/np.median(LRs)
	return list(LRs)

# compute resistant diagnostic
def resistant_diagnostics(data_list,shareddata,njack=200):

	batches = list(range(njack))
	boot_sampling_p=partial(gene_boot_sampling, data_list=data_list,shareddata=shareddata)
	with Pool(processes=4) as pool:
		LRs = pool.map(boot_sampling_p, batches)
	
	LRs = [LR for LR in LRs if isinstance(LR,list)]

	if len(LRs)<50:
		return None, False
	u = []
	for i in range(len(data_list)):
		ui = np.max([l[i] for l in LRs])
		u.append(ui)
	
	unorm=np.median(u)
	RDs = {}
	for i,data in enumerate(data_list):
		RDs[data['bc']] = u[i]/unorm
	return RDs, True

