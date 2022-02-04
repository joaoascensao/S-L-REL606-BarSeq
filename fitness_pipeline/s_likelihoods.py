'''
9.13.20
Joao Ascensao

Functions for knockout fitness likelihoods & helper functions

'''
import pandas as pd
import scipy as sp 
import numpy as np 
from scipy import optimize
from scipy.optimize import minimize
from scipy import interpolate
import scipy.special
import copy
import re
import random
from statsmodels.stats import multitest
import time
from itertools import compress
from multiprocessing import Pool, TimeoutError
from functools import partial
import matplotlib.pyplot as plt
#import sys
#import statsmodels.api as sm

largenum = 1e9

# clean up all data for a gene
def getgenedata(countdf,betas):
	gene_data=betas.copy()
	r1=[]

	for _,row in gene_data.iterrows():
		d1=str(row['Day label'])
		if re.match(r'^-?\d+(?:\.\d+)?$', d1) is not None:
			d10=np.int(np.float(d1))
			d1=str(d10)
		r1.append(np.float(countdf[d1]))
	r2=copy.deepcopy(r1)
	
	#get rid of timepoints after the bc goes extinct
	if r1[-1]==0:
		x=np.array(r1)
		xa=np.flip(np.short(x==0))
		xb=np.flip(np.cumprod(xa))
		wh = np.where(xb==1)
		xb[wh[0][0]]=0
		r1 = x[xb==0]

	
	gene_data_dic = gene_data[['beta + Ne R','R']].to_dict(orient='series')
	gene_data_dic['r']=r1
	gene_data_dic2={}
	for key in gene_data_dic:
		gene_data_dic2[key] = np.longdouble(np.array(gene_data_dic[key]))

	k = gene_data_dic2['beta + Ne R']
	r = gene_data_dic2['r']
	R = gene_data_dic2['R']


	f0guess = r[0]/R[0]
	fl = (r[0] - 2.75*r[0]*np.sqrt(k[0]))/R[0]
	fu = (r[0] + 2.75*r[0]*np.sqrt(k[0]))/R[0]
	if fl<=0:
		fl=1e-10
	if fu>1:
		fu=1
	elif fu<=0:
		fu=50/R[0]
	
	gene_data_dic2['f0range']=(f0guess,fl,fu)

	return {'f0range':gene_data_dic2['f0range'], 'r':gene_data_dic2['r']}

# clean up data that is shared between genes, i.e. variance parameters and mean fitness changes
def getshareddata(betas,xbar_df):
	gene_data=betas.copy()
	xbar=[]
	xbar_df['Day label'] = xbar_df['Day label'].astype(str)
	for _,row in gene_data.iterrows():
		d1=str(row['Day label'])
		if re.match(r'^-?\d+(?:\.\d+)?$', d1) is not None:
			d10=np.int(np.float(d1))
			d1=str(d10)

		
		if row['Day']==0:
			xbar.append(0)
		else:
			xbar.append(np.float(xbar_df[xbar_df['Day label']==d1]['xbar'].mean()))
	gene_data['xbar']=xbar
	gene_data_dic = gene_data[['beta + Ne R','R','gens','xbar']].to_dict(orient='series')
	
	gene_data_dic2={}
	for key in gene_data_dic:
		gene_data_dic2[key] = np.longdouble(np.array(gene_data_dic[key]))
	return gene_data_dic2

# stirlings approx
def gamma_stirlings(x):
	return 0.5*np.log(2*np.pi) + (x - 0.5)*np.log(x) - x + 1/(12*x)

# likelihood for a single barcode, single timepoint
def og_likelihood0(s,gs):
	g1,g2,g3,b1,b2,kappa,r=gs
	if b1: #stirlings approx
		dg1 = gamma_stirlings(g1)
	else:
		dg1 = np.longdouble(sp.special.loggamma(np.float64(g1)))
	
	if b2: #stirlings approx
		dg2 = gamma_stirlings(g2)
	else:
		dg2 = np.longdouble(sp.special.loggamma(np.float64(g2)))
	dg3 = gamma_stirlings(g3)
	return dg1 - dg2 - np.log(kappa)*g2 + r*np.log(kappa - 1) - r*np.log(kappa) - dg3

# add all log likelihoods across time for a single barcode, then integrate over nuisance intercept parameter
def int_likelihood(s,data,fl,fu):
	r,kappa,m0=data
	f0 = np.linspace(fl,fu,10**3)

	ll=0
	for i,ri in enumerate(r):
		m = f0*m0[i]
		g1=ri + m/(kappa[i]-1)
		g2=m/(kappa[i]-1)
		g3=ri+1
		
		b1=np.all(g1>1)
		b2=np.all(g2>1)

		gs=(g1,g2,g3,b1,b2,kappa[i],ri)
		ll0=og_likelihood0(s,gs)
		#print(ll0)
		ll+=ll0
	
	cond=np.isfinite(ll)
	f0 = f0[cond]
	ll = ll[cond]
	if len(ll)<=10:
		return np.nan
	mm=np.max(ll)	
	li = np.log(np.trapz(np.exp(ll-mm),x=f0)) + mm
	return li

# add likelihoods over all barcodes within a gene
def add_likelihoods(s,data_list,shareddata):
	lls=0
	s=np.longdouble(s)
	k = shareddata['beta + Ne R']
	R = shareddata['R']
	xbar = shareddata['xbar']
	t = shareddata['gens']
	#print(shareddata,data_list)
	for data in data_list:
		r = data['r']
		_,fl,fu = data['f0range']
		m0=R*np.exp((s-xbar)*t)
		

		dd = (r,k,m0)
		#lls0 = taylor_int_likelihood(s,dd,f0guess,fl,fu)
		lls0 = int_likelihood(s,dd,fl,fu)
		if not np.isfinite(lls0):
			return largenum
		lls+=lls0
	if not np.isfinite(lls):
		return largenum
	return np.float(lls*-1)

# get bracket for a minimization of negative log-likelihood
def get_bracket(data_list,shareddata,objective,delta=0.03,nmax=5,combo=False):
	eps=1e-9
	s = 0
	c = 0
	if combo:
		for replicate in data_list:
			for _,data in enumerate(data_list[replicate]):
				r = data['r']
				R = shareddata[replicate]['R']
				gens=shareddata[replicate]['gens']
				s0=(np.log(r[2]/R[2] + eps) - np.log(r[0]/R[0] + eps))/gens[2]
				if np.abs(s0)>0.5:
					s0=0
				s+=s0
				c+=1
	else:
		for _,data in enumerate(data_list):
			r = data['r']
			R = shareddata['R']
			gens=shareddata['gens']
			s0=(np.log(r[2]/R[2] + eps) - np.log(r[0]/R[0] + eps))/gens[2]
			if np.abs(s0)>0.5:
				s0=0
			s+=s0
			c+=1
	sm=s/c
	sf=objective(sm)
	u = sm + delta
	uf = objective(u)
	c=0
	while uf<sf*1.02:
		u = u + delta
		uf = objective(u)
		if c==nmax:
			return np.nan, np.nan, False
		c+=1
	
	l = sm - delta
	lf = objective(l)
	c=0
	while lf<sf*1.02:
		l = l - delta
		lf = objective(l)
		if c==nmax:
			return np.nan, np.nan, False
		c+=1
	if not np.all(np.isfinite([uf,sf,lf])):
		return np.nan, np.nan, False
	if not np.all(np.isfinite([u,sm,l])):
		return np.nan, np.nan, False
	if sm>u or sm<l:
		return np.nan, np.nan, False
	if sf==largenum or lf==largenum or uf==largenum:
		return np.nan, np.nan, False
	
	return (l,sm,u), (lf,sf,uf), True

def get_boot_bracket(bracket,objective,delta=0.02,nmax=5):
	sm=bracket[1]
	sf=objective(sm)
	u = bracket[2]
	uf = objective(u)
	c=0
	while uf<sf*1.02:
		u = u + delta
		uf = objective(u)
		if c==nmax:
			return np.nan, np.nan, False
		c+=1
	
	l = bracket[0]
	lf = objective(l)
	c=0
	while lf<sf*1.02:
		l = l - delta
		lf = objective(l)
		if c==nmax:
			return np.nan, np.nan, False
		c+=1
	if not np.all(np.isfinite([uf,sf,lf])):
		return np.nan, np.nan, False
	if not np.all(np.isfinite([u,sm,l])):
		return np.nan, np.nan, False
	if sm>u or sm<l:
		return np.nan, np.nan, False
	if sf==largenum or lf==largenum or uf==largenum:
		return np.nan, np.nan, False
	
	return (l,sm,u), (lf,sf,uf), True

# get standard error from likelihood via numerical approx of second derivative of log-likelihood
def std_from_likelihood(optim_s, optim_ll, get_ll, h=1e-5):
	deriv2 = (-get_ll(optim_s + 2*h) + 16*get_ll(optim_s + h) - 30*optim_ll + 16*get_ll(optim_s - h) - get_ll(optim_s - 2*h))/(12*(h**2))
	std = 1/np.sqrt(deriv2)
	return std

# p-value calculation
def posterior_pval(s,std,obj):
	lowbound=s-6*std
	upbound=s+6*std
	if lowbound>0:
		lowbound=-1e-9
	if upbound<0:
		upbound=1e-9
	svec = np.concatenate((np.linspace(lowbound, s-2*std-0.001,50),np.linspace(s-2*std,s+2*std,10**2),np.linspace(s+2*std+0.001,upbound,50)),axis=None)
	logpost = np.array([obj(si)*-1 for si in svec])
	
	log_post_int= interpolate.interp1d(svec,logpost,kind='cubic')
	svec1 = np.linspace(min(svec),max(svec),5*10**4)
	ll0 = log_post_int(0)
	log_post2 = log_post_int(svec1)

	lam = np.nan_to_num(ll0-log_post2)
	mm=np.max(log_post2)
	post2=np.nan_to_num(np.exp(log_post2-mm))
	post3=post2/np.sum(post2)

	pval=np.sum(post3[lam>0])
	return pval

# get maximum likelihood estimate of fitness
def s_maxlikelihood(data_list, shareddata, boot=False, bounds=(-1.5,0.8), boot_bracket=False):
	#start=time.time()
	success=False
	objective = lambda s: add_likelihoods(s, data_list, shareddata)

	if not boot:
		bracket,bracket_eval,bracket_success=get_bracket(data_list,shareddata,objective)
	else:
		bracket,bracket_eval,bracket_success=get_boot_bracket(boot_bracket,objective)
	
	if bracket_success:
		s,sfun,success = minimize_parabolic_interp(objective, bracket, bracket_eval)
	if not success:
		res = sp.optimize.minimize_scalar(objective,bounds=bounds,method='Bounded',options={'xatol': 1e-4})
		s=res.x
		sfun=res.fun
		success=res.success
	if not success:
		res = sp.optimize.dual_annealing(objective,bounds=bounds)
		s=res.x
		sfun=res.fun
		success=res.success
		
	if success and sfun!=largenum and s>(bounds[0]+0.01) and s<(bounds[1]-0.01):
		if boot:
			return s, sfun, True
		s_std=std_from_likelihood(s, sfun, objective)
		try:
			s_pval = posterior_pval(s,s_std,objective)
		except:
			s_pval = 1
		return s, s_std, s_pval, True
	else:
		if boot:
			return None, None, False
		return np.nan, np.nan, np.nan, False

# BH FDR correaction
def FDR_correction(pvals,alpha=0.05):
	rejected,pvals_corr,_,_=multitest.multipletests(pvals, alpha=alpha, method='fdr_bh', is_sorted=False)
	return [int(i) for i in rejected], pvals_corr

# more efficient minimization via parabolic interpolation
def minimize_parabolic_interp(objective, bracket, bracket_eval, maxiter=500, xtol=1e-4):
	eps=1e-10
	nobj=lambda x: objective(x)*-1
	x = copy.deepcopy(bracket)
	y = tuple([b*-1 for b in bracket_eval])
	for _ in range(maxiter):
		
		if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
			return np.nan, np.nan, False
		den = (x[1] - x[0])*(y[1] - y[2]) - (x[1] - x[2])*(y[1] - y[0])
		if den==0:
			return np.nan, np.nan, False
		x3 = x[1] - 0.5*(((x[1] - x[0])**2)*(y[1] - y[2]) - ((x[1] - x[2])**2)*(y[1] - y[0]))/den
		y3 = nobj(x3)
		
		if not np.isfinite(x3) or not np.isfinite(y3):
			return np.nan, np.nan, False

		if np.abs((x3-x[1])/(x[1] + eps))<xtol:
			return x3, y3*-1, True
		
		if x3>x[1]:
			x = (x[1],x3,x[2])
			y = (y[1],y3,y[2])
		else:
			x = (x[0],x3,x[1])
			y = (y[0],y3,y[1])
	return np.nan, np.nan, False
