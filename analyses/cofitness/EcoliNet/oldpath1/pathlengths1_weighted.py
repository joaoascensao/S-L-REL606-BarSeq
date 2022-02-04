import pandas as pd
import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import random
import networkx as nx
import time
sns.set(style="whitegrid",font_scale=1.2)


gn = pd.read_csv('ecoli_net_B.csv').drop(columns=['Unnamed: 0'])
cc = pd.read_csv('gnpathlengths_community_changes.csv').drop(columns=['Unnamed: 0'])
sp = pd.read_csv('shortest_paths_ecolinet.csv').rename(columns={'source':'gene_ID 1','target':'gene_ID 2'}).drop(columns=['Unnamed: 0'])
sp = sp[sp['distance']==1]
cc = cc[cc['distance']==1]
df = sp.merge(gn,on=['gene_ID 1','gene_ID 2'],how='right')

d1 = cc.merge(df.rename(columns={'gene_ID 1':'gene A','gene_ID 2':'gene B'}), how='inner',on=['gene A','gene B'])
d2 = cc.merge(df.rename(columns={'gene_ID 2':'gene A','gene_ID 1':'gene B'}), how='inner',on=['gene A','gene B'])

pd.concat([d1,d2],ignore_index=True).to_csv('cc_weightedpaths.csv')