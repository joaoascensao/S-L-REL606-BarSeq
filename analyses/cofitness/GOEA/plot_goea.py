
import pandas as pd
import scipy as sp 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
import sys
sns.set(style="whitegrid",font_scale=1)

strains=['S','L','R']


for exp in strains:
	for side in range(10):
		for w in [0,1]:
			try:
				df=pd.read_excel('res/GOEnrichment_{}_{}.xlsx'.format(exp,w))
			except:
				continue
			cc='navy'

			df=df[df['NS']=='BP']
			

			df=df[df['study_count']>0]
			df=df[df['name']!='biological_process'].sort_values("p_fdr_bh")

			print(df)

			nc=min((20,len(df)))
			print(nc)

			if len(df)==0:
				continue

			plt.figure()
			sns.barplot(y='name',x='study_count',data=df.head(20),palette=sns.light_palette(cc, reverse=True,n_colors=nc),ci=None)
			plt.ylabel('')
			plt.xlabel('# Genes')
			plt.tight_layout()
			plt.savefig('plots/GOEA_{}_{}.png'.format(exp,w+1),dpi=600)