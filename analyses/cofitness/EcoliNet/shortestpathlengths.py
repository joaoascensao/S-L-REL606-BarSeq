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


gn = pd.read_csv('ecoli_net_B.csv')

G = nx.from_pandas_edgelist(gn, source='gene_ID 1', target='gene_ID 2')
shortest_paths = nx.shortest_path_length(G)


pathlen_df = []
sources=[]
targets=[]
ds=[]
for item in shortest_paths:
	source = item[0]
	for target in item[1]:
		if source!=target:
			sources.append(source)
			targets.append(target)
			ds.append(item[1][target])

pathlen_df=pd.DataFrame({
		'source':sources,
		'target':targets,
		'distance':ds,
		})
pathlen_df.to_csv('shortest_paths_ecolinet.csv')