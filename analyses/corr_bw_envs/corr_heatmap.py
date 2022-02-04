'''

Calculate correlation of fitness effects between environments and spit out heatmap (Fig 3)

'''

import pandas as pd
import scipy as sp 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
import WeightedCorr as wc
import scipy.spatial.distance

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

minstd=0.01
links={}
for i,strain in enumerate(strain_names):
	
	data=[]
	n=len(exps[strain][0])
	corrmat = np.zeros((n,n))

	for jA,exp0A in enumerate(exps[strain][0]):
		exp1A=exps[strain][1][jA]
		fitA = pd.read_csv(dirs[exp0A] + 'data/fitness/{}_fitness.csv'.format(exp1A)).dropna(subset=['s','s std']).rename(columns={'s':'sA','s std':'stdA'})
		fitA = fitA[fitA['stdA']<minstd]
		for jB in range(jA,len(exps[strain][0])):
			if jA!=jB:
				exp0B=exps[strain][0][jB]
				exp1B=exps[strain][1][jB]
				fitB = pd.read_csv(dirs[exp0B] + 'data/fitness/{}_fitness.csv'.format(exp1B)).dropna(subset=['s','s std']).rename(columns={'s':'sB','s std':'stdB'})
				fitB = fitB[fitB['stdB']<minstd]
				fit = fitA.merge(fitB,on='gene_ID',how='inner')
				fit['w'] = 1/(fit['stdA']**2 + fit['stdB']**2)
				rho = wc.wpearson(np.array(fit['sA']), np.array(fit['sB']), np.array(fit['w']))
				if rho<0 and rho>-0.0006:
					rho=0
				corrmat[jA,jB]=rho
				corrmat[jB,jA]=rho

	names = [expnames[o] for o in exps[strain][1]]
	distmat=(1-corrmat)*(1-np.diag(np.ones(n)))
	pp = sp.spatial.distance.squareform(distmat.clip(min=0))
	link = scipy.cluster.hierarchy.linkage(pp,optimal_ordering=True,method='ward')

	plt.figure()
	np.fill_diagonal(corrmat,None)
	corrmat_df = pd.DataFrame(corrmat, index=names, columns=names)
	cg=sns.clustermap(corrmat_df,row_linkage=link,col_linkage=link,cmap=cmaps[strain],
		**{'cbar':False,'linewidths':1.5,'vmin':0,'vmax':0.65,'linecolor':darkgrey,'annot':True})
	plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
	plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=25)
	links[strain]=link


	plt.tight_layout()
	plt.savefig('corr_heatmap_{}.png'.format(strain),dpi=600)



# cross-correlations

combos=[('R','S'),('R','L'),('S','L')]
sns.set(style="whitegrid",font_scale=1.1)
for i,combo in enumerate(combos):
	strain1,strain2 = combo
	n1=len(exps[strain1][0])
	n2=len(exps[strain2][0])
	corrmat = np.zeros((n1,n2))
	for jA,exp0A in enumerate(exps[strain1][0]):
		exp1A=exps[strain1][1][jA]
		fitA = pd.read_csv(dirs[exp0A] + 'data/fitness/{}_fitness.csv'.format(exp1A)).rename(columns={'s':'sA','s std':'stdA'})
		fitA = fitA[fitA['stdA']<minstd]
		for jB,exp0B in enumerate(exps[strain2][0]):
			exp0B=exps[strain][0][jB]
			exp1B=exps[strain][1][jB]
			fitB = pd.read_csv(dirs[exp0B] + 'data/fitness/{}_fitness.csv'.format(exp1B)).rename(columns={'s':'sB','s std':'stdB'})
			fitB = fitB[fitB['stdB']<minstd]
			fit = fitA.merge(fitB,on='gene_ID',how='inner')
			fit['w'] = 1/(fit['stdA']**2 + fit['stdB']**2)
			#rho = wc.WeightedCorr(x=fit['sA'], y=fit['sB'], w=fit['w'] )(method='pearson')
			rho = wc.wpearson(np.array(fit['sA']), np.array(fit['sB']), np.array(fit['w']))
			corrmat[jA,jB]=rho
	names1 = [expnames[o] for o in exps[strain1][1]]
	names2 = [expnames[o] for o in exps[strain2][1]]
	corrmat_df = pd.DataFrame(corrmat, index=names1, columns=names2)
	plt.figure()
	cg=sns.clustermap(corrmat_df,row_linkage=links[strain1],col_linkage=links[strain2],cmap="pink_r",
		**{'cbar':False,'linewidths':1.5,'vmin':0,'vmax':0.4,'linecolor':darkgrey,'annot':True})
	plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
	plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=25)

	plt.tight_layout()
	plt.savefig('cross_corr_heatmap_{}_{}.png'.format(strain1,strain2),dpi=600)


