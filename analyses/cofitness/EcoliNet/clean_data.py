import pandas as pd
import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="whitegrid",font_scale=1.2)

k2b = pd.read_csv('genes_K2B.csv').rename(columns={'K Gene Name':'Symbol'})
ecoli_net = pd.read_csv('EcoliNet.v1.txt',sep='\t')
k12 = pd.read_csv('K12_gene_result.txt',sep='\t')
k12.dropna(subset=['Aliases'],inplace=True)
k12['orf']=k12['Aliases'].apply(lambda row: row.split(',')[0])

d1 = k2b.merge(k12[['orf','Symbol']],on='Symbol',how='inner')
d2 = ecoli_net.merge(d1.rename(columns={'orf':'orf1','gene_ID':'gene_ID 1'}),how='left',on='orf1')
d3 = d2.merge(d1.rename(columns={'orf':'orf2','gene_ID':'gene_ID 2'}),how='left',on='orf2')


d4=d3[['gene_ID 1','gene_ID 2','ll']]
d4.dropna(inplace=True)
d4.to_csv('ecoli_net_B.csv')