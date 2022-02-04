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
operons= pd.read_csv('../../../Good_etal_2017/door_operon_list.txt',sep='\t').rename(columns={'Synonym':'gene_ID 1'})
operons['gene_ID 2']=operons['gene_ID 1']
nw={}
for strain in strains:
	nw[strain] = pd.read_csv('../../{}_cofitness.csv'.format(strain))[['gene_ID 1','gene_ID 2','corr']].rename(columns={'corr':'corr{}'.format(strain)})

stderr= lambda th: (1-th[0]**2)/np.sqrt(th[1]-2)


data=[]
strainpairs=[('R','S'),('R','L'),('S','L')]

def sameoperonfunc(row):
	if np.int(row['OperonID 1'])==np.int(row['OperonID 2']):
		return True
	else:
		return False

for strain1,strain2 in strainpairs:
	df = nw[strain1].merge(nw[strain2],on=['gene_ID 1','gene_ID 2'],how='inner')
	df['id']=list(range(len(df)))
	#rho,p=sp.stats.pearsonr(df0['corr{}'.format(strain1)],df0['corr{}'.format(strain2)])
	print(strain1,strain2)
	df0 = df.merge(operons[['gene_ID 1','OperonID']].rename(columns={'OperonID':'OperonID 1'}),on=['gene_ID 1'],how='inner')
	df0 = df0.merge(operons[['gene_ID 2','OperonID']].rename(columns={'OperonID':'OperonID 2'}),on=['gene_ID 2'],how='inner')
	'''
	sameoperon=[]
	for pp,row in df0.iterrows():
		if np.int(row['OperonID 1'])==np.int(row['OperonID 2']):
			sameoperon.append(True)
		else:
			sameoperon.append(False)

	df0['Same Operon']=sameoperon
	'''
	df0['Same Operon'] = df0[['OperonID 1','OperonID 2']].apply(lambda row: sameoperonfunc(row), axis=1)
	df1 = df0[df0['Same Operon']==True]
	rho,p=sp.stats.pearsonr(df1['corr{}'.format(strain1)],df1['corr{}'.format(strain2)])
	data.append(pd.DataFrame({
		'strains':strain1+', '+strain2,
		'rho':rho,
		'stderr':stderr([rho,len(df1)]),
		'Same Operon':True,
		},index=[0]))
	
	df1 = df0[df0['Same Operon']==False]
	rho,p=sp.stats.pearsonr(df1['corr{}'.format(strain1)],df1['corr{}'.format(strain2)])
	data.append(pd.DataFrame({
		'strains':strain1+', '+strain2,
		'rho':rho,
		'stderr':stderr([rho,len(df1)]),
		'Same Operon':False,
		},index=[0]))


red='#d91e3a'
blue='#2086e6'
dft = pd.concat(data)
print(dft)


darkgrey='#40464d'
plt.figure()
g=sns.barplot(data=dft,x='strains',y='rho',hue='Same Operon', palette=[blue,red])#
x_coords = [p.get_x() + 0.5*p.get_width() for p in g.patches]
y_coords = [p.get_height() for p in g.patches]

errors = pd.DataFrame({'x':x_coords,'rho':y_coords}).merge(dft[['rho','stderr']],on='rho',how='inner')
plt.errorbar(errors['x'], errors['rho'], yerr=errors['stderr'], fmt=' ',zorder=np.inf,ecolor=darkgrey)

plt.ylabel(r'Correlation between cofitness')
plt.xlabel('Strain Pairs')
plt.tight_layout()
plt.savefig('corr_operons.png',dpi=600)
