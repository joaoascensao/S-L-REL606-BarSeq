import pandas as pd
import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import networkx as nx
from networkx.algorithms import community
import copy
import random

darkgrey='#40464d'
sns.set(style="whitegrid",font_scale=1.2,rc= {'patch.edgecolor': darkgrey})

df = pd.read_csv('cc_weightedpaths.csv')

#for i,group in df.groupby(by=['strain 1','strain 2','distance']):
#	print(i)
#	print(np.sum(group['connected'])/len(group))

df.dropna(subset=['ll'])

nmax=8
lls=list(np.array(df['ll']))
cats=[]
spa = (np.max(lls) - np.min(lls))/nmax
ra = np.min(lls) + np.array(range(nmax+1))*spa
midpoints = [(ra[i]+ra[i+1])/2 for i in range(len(ra)-1)]

print(midpoints)

d2=[]
ra2=ra[:-1]
for lli in lls:
	ii = np.where(lli>=ra2)[0][-1]
	d2.append(midpoints[ii])

df['d2']=d2
#df['Distance'] = [np.int(d) for d in list(np.array(df[df['distance']<=4]['distance']))]
#print(df['Distance'])
df['x'] = df['strain 1'] + df['strain 2']
df['x'] = df['x'].apply(lambda x: x[0]+', '+x[1])
print(df['x'])
plt.figure()

g=sns.barplot(data=df,x='x',y='connected',hue='d2',ci=68.26, palette='summer_r')#pink_r
g.figure.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=min(lls), vmax=max(lls), clip=False), cmap='summer_r'),label='EcoliNet Score')
plt.xlabel('Strain 1, Strain 2')
plt.ylabel('')
plt.title('Probability two genes are in the same cluster in strain 2,\nif they were together in strain 1')
plt.legend('',frameon=False)
plt.tight_layout()
plt.savefig('pathlengths_weighted1.png',dpi=500)
