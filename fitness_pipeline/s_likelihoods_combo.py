'''
9.13.20
Joao Ascensao

Functions for likelihoods & helper functions


'''
import pandas as pd
import scipy as sp 
import numpy as np 
from scipy import optimize
from scipy.optimize import minimize
from scipy import interpolate
import scipy.special
import matplotlib.pyplot as plt
import copy
import re
import random
from statsmodels.stats import multitest
import time
from itertools import compress
from multiprocessing import Pool, TimeoutError
from functools import partial
import s_likelihoods as lh

largenum = 1e9


def add_likelihoods_acrossreps(s, data_list, shareddata):
	nll=0
	for replicate in data_list:
		nll+=lh.add_likelihoods(s, data_list[replicate], shareddata[replicate])
	return nll

def s_maxlikelihood_combo(data_list, shareddata, bounds=(-1.5,0.8)):
	#start=time.time()
	success=False
	objective = lambda s: add_likelihoods_acrossreps(s, data_list, shareddata)

	bracket,bracket_eval,bracket_success=lh.get_bracket(data_list,shareddata,objective,combo=True)
	
	if bracket_success:
		s,sfun,success = lh.minimize_parabolic_interp(objective, bracket, bracket_eval)
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
		s_std=lh.std_from_likelihood(s, sfun, objective)
		s_pval = lh.posterior_pval(s,s_std,objective)
		return s, s_std, s_pval
	else:
		return np.nan, np.nan, np.nan
