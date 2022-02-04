import pandas as pd
import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import random
import time
import copy
sns.set(style="whitegrid",font_scale=1.2)

strains = ['R','S','L']
gn=pd.read_csv('../ecoli_net_B.csv')
nw={}
for strain in strains:
	nw[strain] = pd.read_csv('../../{}_cofitness.csv'.format(strain))[['gene_ID 1','gene_ID 2','corr']].rename(columns={'corr':'corr{}'.format(strain)})

stderr= lambda th: (1-th[0]**2)/np.sqrt(th[1]-2)


data=[]
strainpairs=[('R','S'),('R','L'),('S','L')]
for strain1,strain2 in strainpairs:
	df = nw[strain1].merge(nw[strain2],on=['gene_ID 1','gene_ID 2'],how='inner')
	df['id']=list(range(len(df)))
	#rho,p=sp.stats.pearsonr(df0['corr{}'.format(strain1)],df0['corr{}'.format(strain2)])
	print(strain1,strain2)
	df0 = df.merge(gn,on=['gene_ID 1','gene_ID 2'],how='inner')

	df1 = df0[df0['ll']>5.5]
	rho,p=sp.stats.pearsonr(df1['corr{}'.format(strain1)],df1['corr{}'.format(strain2)])
	data.append(pd.DataFrame({
		'strains':strain1+', '+strain2,
		'rho':rho,
		'stderr':stderr([rho,len(df1)]),
		'cat':3,
		},index=[0]))
	
	df1 = df0[df0['ll']<3]
	rho,p=sp.stats.pearsonr(df1['corr{}'.format(strain1)],df1['corr{}'.format(strain2)])
	data.append(pd.DataFrame({
		'strains':strain1+', '+strain2,
		'rho':rho,
		'stderr':stderr([rho,len(df1)]),
		'cat':1,
		},index=[0]))


	df1 = df0[(df0['ll']>3) & (df0['ll']<5.5)]
	rho,p=sp.stats.pearsonr(df1['corr{}'.format(strain1)],df1['corr{}'.format(strain2)])
	data.append(pd.DataFrame({
		'strains':strain1+', '+strain2,
		'rho':rho,
		'stderr':stderr([rho,len(df1)]),
		'cat':2,
		},index=[0]))

	df2 = df.set_index('id').drop(index=df0['id'])
	rho,p=sp.stats.pearsonr(df2['corr{}'.format(strain1)],df2['corr{}'.format(strain2)])
	data.append(pd.DataFrame({
		'strains':strain1+', '+strain2,
		'rho':rho,
		'stderr':stderr([rho,len(df2)]),
		'cat':0,
		},index=[0]))

dft = pd.concat(data)
print(dft)


darkgrey='#40464d'
plt.figure()
g=sns.barplot(data=dft,x='strains',y='rho',hue='cat', palette='summer_r')#
g.legend(labels=['0',r'$x<3$',r'$3<x<5.5$',r'$x>5.5$'],title=r'EcoliNet Score, $x$',handles=[g.patches[0],g.patches[3],g.patches[6],g.patches[9]])
x_coords = [p.get_x() + 0.5*p.get_width() for p in g.patches]
y_coords = [p.get_height() for p in g.patches]

errors = pd.DataFrame({'x':x_coords,'rho':y_coords}).merge(dft[['rho','stderr']],on='rho',how='inner')
plt.errorbar(errors['x'], errors['rho'], yerr=errors['stderr'], fmt=' ',zorder=np.inf,ecolor=darkgrey)

plt.ylabel(r'Correlation between cofitness')
plt.xlabel('Strain Pairs')
plt.tight_layout()
plt.savefig('corr_cofit.png',dpi=600)
