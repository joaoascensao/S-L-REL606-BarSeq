import pandas as pd
import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from networkx.algorithms import community
import copy
import random
sns.set(style="whitegrid",font_scale=1.2)

df = pd.read_csv('gnpathlengths_community_changes.csv')

for i,group in df.groupby(by=['strain 1','strain 2','distance']):
	print(i)
	print(np.sum(group['connected'])/len(group))

df.dropna(subset=['distance'],inplace=True)
df = df[df['distance']<=4]
df['Distance'] = [np.int(d) for d in list(np.array(df['distance']))]
print(df['Distance'])
df['x'] = df['strain 1'] + df['strain 2']
df['x'] = df['x'].apply(lambda x: x[0]+', '+x[1])
print(df['x'])
plt.figure()
sns.barplot(data=df,x='x',y='connected',hue='Distance',ci=68.26)
plt.xlabel('Strain 1, Strain 2')
plt.ylabel('')
plt.title('Probability two genes are in the same cluster in strain 2,\nif they were together in strain 1')
plt.tight_layout()
plt.savefig('gnpathlengths_community_changes.png',dpi=500)
