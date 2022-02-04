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
gn['weight'] = 1/gn['ll']
G1 = nx.from_pandas_edgelist(gn, source='gene_ID 1', target='gene_ID 2')

strains=['R','S','L']
for strain in strains:
	df = pd.read_csv('../{}_cofitness_sig.csv'.format(strain))
	G2 = nx.from_pandas_edgelist(df, source='gene_ID 1', target='gene_ID 2')
	res= nx.graph_edit_distance(G1,G2,node_del_cost=0)
	print(strain,res)
	del G2