'''

Bootstrapping procedure to calculate support for each branch

'''

import pandas as pd
import scipy as sp 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
import WeightedCorr as wc
import scipy.spatial.distance
import random

darkgrey='#40464d'
sns.set(style="whitegrid",font_scale=1.4)

red='#d91e3a'
blue='#2086e6'
darkgrey='#40464d'


strain_names={
		'R':'REL606',
		'S':'6.5k S',
		'L':'6.5k L',
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


cmaps = {
	'R':'Blues',
	'S':'Oranges',
	'L':'Greens',
}

def get_linkdic(link,n):
	linkdic={}
	v={}
	vb={}
	p={}
	for k,row in enumerate(link[:-1]):
		linkdic[len(link)+k+1] = []
		a = int(row[0])
		b = int(row[1])
		for c in [a,b]:
			if c > len(link):
				linkdic[len(link)+k+1] += linkdic[c]
			else:
				linkdic[len(link)+k+1].append(c)
		ss = linkdic[len(link)+k+1]
		linkdic[len(link)+k+1] = set(ss)
		x = np.zeros(n)
		for s in ss:
			x[s] = 1
		v[len(link)+k+1] = x
		vb[len(link)+k+1] = np.abs(1-x)
		p[len(link)+k+1] = min(np.sum(x),np.sum(np.abs(1-x)))
	return linkdic, v, vb, p

def unison_resample_gids(gids,gid_boot,a,b,c):
    p=[gids.index(gid) for gid in gid_boot]
    return a[p], b[p], c[p]

def bootstrap(exps,strain,linkbundle,fitmap,intgid,nboot=5000):
	hamm = lambda a: np.count_nonzero(a[0]!=a[1])
	linkdic0, v0, vb0, p0 = linkbundle
	n=len(exps[strain][0])
	bootres={}
	for key in linkdic0:
		bootres[key]=0

	for i in range(nboot):
		corrmat = np.zeros((n,n))
		gid_boot = random.choices(intgid,k=len(intgid))
		for jA,exp0A in enumerate(exps[strain][0]):
			for jB in range(jA,len(exps[strain][0])):
				if jA!=jB:
					a0,b0,w0,gids=fitmap[(jA,jB)]
					a,b,w=unison_resample_gids(gids,gid_boot,a0,b0,w0)
					rho = wc.wpearson(a,b,w)
					if rho<0:
						rho=0
					corrmat[jA,jB]=rho
					corrmat[jB,jA]=rho
		distmat=(1-corrmat)*(1-np.diag(np.ones(n)))
		pp = sp.spatial.distance.squareform(distmat.clip(min=0))
		link = scipy.cluster.hierarchy.linkage(pp,method='ward')
		linkbundle_boot = get_linkdic(link,n)
		linkdic_boot, v_boot, _, _ = linkbundle_boot 
		for key1 in linkdic0:
			if p0[key1]<=2:
				for key2 in linkdic_boot:
					if linkdic0[key1] == linkdic_boot[key2]:
						bootres[key1] += 1
						break
			else:
				delta=[]
				for key2 in linkdic_boot:
					dd = min( hamm((v0[key1],v_boot[key2])), hamm((vb0[key1],v_boot[key2])) )
					delta.append(dd)
				bootres[key1] += 1 - np.min(delta)/(p0[key1] - 1)

	for key in bootres:
		bootres[key] = bootres[key]/nboot
	return bootres



minstd=0.01
links={}
for i,strain in enumerate(strain_names):
	if strain=='R':
		continue
	data=[]
	n=len(exps[strain][0])
	corrmat = np.zeros((n,n))
	fitmap={}
	gids=[]
	for jA,exp0A in enumerate(exps[strain][0]):
		exp1A=exps[strain][1][jA]
		fitA = pd.read_csv(dirs[exp0A] + 'data/fitness/{}_fitness.csv'.format(exp1A)).dropna(subset=['s','s std']).rename(columns={'s':'sA','s std':'stdA'})
		fitA = fitA[fitA['stdA']<minstd]
		for jB in range(jA,len(exps[strain][0])):
			if jA==jB:
				fitA = pd.read_csv(dirs[exp0A] + 'data/fitness/{}1_fitness.csv'.format(exp1A)).dropna(subset=['s','s std']).rename(columns={'s':'sA','s std':'stdA'})
				fitB = pd.read_csv(dirs[exp0A] + 'data/fitness/{}2_fitness.csv'.format(exp1A)).dropna(subset=['s','s std']).rename(columns={'s':'sB','s std':'stdB'})
				fit = fitA.merge(fitB,on='gene_ID',how='inner')
				fit['w'] = 1/(fit['stdA']**2 + fit['stdB']**2)
				rho = wc.wpearson(np.array(fit['sA']), np.array(fit['sB']), np.array(fit['w']))
				#corrmat[jA,jB]=rho
				#pass
			else:
				exp0B=exps[strain][0][jB]
				exp1B=exps[strain][1][jB]
				fitB = pd.read_csv(dirs[exp0B] + 'data/fitness/{}_fitness.csv'.format(exp1B)).dropna(subset=['s','s std']).rename(columns={'s':'sB','s std':'stdB'})
				fitB = fitB[fitB['stdB']<minstd]
				fit = fitA.merge(fitB,on='gene_ID',how='inner').reset_index()
				fit['w'] = 1/(fit['stdA']**2 + fit['stdB']**2)
				rho = wc.wpearson(np.array(fit['sA']), np.array(fit['sB']), np.array(fit['w']))
				if rho<0 and rho>-0.0006:
					rho=0
				corrmat[jA,jB]=rho
				corrmat[jB,jA]=rho
				gid=list(fit['gene_ID'])
				gids.append(set(gid))
				fitmap[(jA,jB)] = (np.array(fit['sA']), np.array(fit['sB']), np.array(fit['w']),gid)#,np.array(fit['stdA']),np.array(fit['stdB'])
				

	for l,gid in enumerate(gids):
		if l==0:
			intgid = gid
		else:
			intgid = intgid & gid

	intgid = list(intgid)
	print(len(intgid))

	names = {ll:expnames[o] for ll,o in enumerate(exps[strain][1])}
	distmat=(1-corrmat)*(1-np.diag(np.ones(n)))
	pp = sp.spatial.distance.squareform(distmat.clip(min=0))
	link = scipy.cluster.hierarchy.linkage(pp,method='ward')

	print(names)
	print(link)
	linkbundle = get_linkdic(link,len(names))
	bootres = bootstrap(exps,strain,linkbundle,fitmap,intgid)
	print(bootres)



	




