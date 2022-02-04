import pandas as pd
import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import random
import networkx as nx
from networkx.algorithms import shortest_paths
import time
sns.set(style="whitegrid",font_scale=1.2)

strains = ['S','L','R']

gn = pd.read_csv('ecoli_net_B.csv')
pathlen_df = pd.read_csv('shortest_paths_ecolinet.csv').drop(columns=['Unnamed: 0'])
#pathlen_dic = {}
#for i,row in pathlen_df.iterrows():
#	pathlen_dic[(str(row['source']),str(row['target']))] = np.int(np.float(row['distance']))
#	pathlen_dic[(str(row['target']),str(row['source']))] = np.int(np.float(row['distance']))
pathlen_dic=pathlen_df.set_index(['source','target']).to_dict('index')
G = nx.from_pandas_edgelist(gn, source='gene_ID 1', target='gene_ID 2')

all_communities={}
needed_genes=[]
genes_by_strain={}
for strain in strains:
	df = pd.read_csv('../{}_fluid_communities.csv'.format(strain))
	ngroups = df.groupby(by='cluster').ngroups
	all_communities[strain]=[]
	gbi=[]
	for n in range(ngroups):
		gi=list(df[df['cluster']==n]['gene_ID'])
		all_communities[strain].append(gi)
		needed_genes+=gi
		gbi+=gi
	genes_by_strain[strain]=list(set(gbi))





strain1_v=[]
strain2_v=[]
gAv=[]
gBv=[]
conn_v=[]
dist_v=[]
comm1_v=[]
comm2_v=[]
comm2b_v=[]



for i,geneA in enumerate(geneAs):
	start=time.time()
	for strain1 in strains:
		for strain2 in strains:
			if strain1!=strain2:
				# find community where geneA is 
				for i1,comm1 in enumerate(all_communities[strain1]):
					if geneA in comm1:
						for i2,comm2 in enumerate(all_communities[strain2]):
							if geneA in comm2:
								# choose some other random gene in comm1
								geneC = random.choice(comm1)
								cc=0
								while (not (geneC in genes_by_strain[strain2])) or (not (geneC in gg)):
									geneC = random.choice(comm1)
									cc+=1
									if cc==5:
										break
								if cc==5:
									break
								if geneC in comm2:
									connect = 1
									i2b=i2
								else:
									connect = 0
									i2b=None
									for i2b,comm2 in enumerate(all_communities[strain2]):
										if geneC in comm2:
											break
								try:
									distance = pathlen_dic[(geneA,geneC)]['distance']
								except:
									distance=np.nan

								strain1_v.append(strain1)
								strain2_v.append(strain2)
								gAv.append(geneA)
								gBv.append(geneC)
								conn_v.append(connect)
								dist_v.append(distance)
								comm1_v.append(i1)
								comm2_v.append(i2)
								comm2b_v.append(i2b)
								break
						break
	#print(time.time()-start)



df_save = pd.DataFrame({
				'strain 1':strain1_v,
				'strain 2':strain2_v,
				'gene A':gAv,
				'gene B':gBv,
				'connected':conn_v,
				'distance':dist_v,
				'cluster 1':comm1_v,
				'cluster 2':comm2_v,
				'cluster 2B':comm2b_v,
			}).drop_duplicates()
df_save.to_csv('gnpathlengths_community_changes.csv')
	
	