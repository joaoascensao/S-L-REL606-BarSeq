'''
7.1.19
Joao Ascensao

Perform a Gene Ontology Enrichment analysis on a set of genes that differ in their fitness effect
between S L and R, listGeneEffects_compare objects
Only include files with a minimum number of rows, minrows
'''
import pandas as pd
import scipy as sp 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import convert_Ecoli_geneIDs
import GOEA
from matplotlib.backends.backend_pdf import PdfPages

sns.set(style="whitegrid")



strains=['S','L','R']

symboldic=convert_Ecoli_geneIDs.uniprot_to_ncbi()

genes = pd.read_csv('../../../rel606_genes_wrangled.csv')

for strain in strains:
	df=pd.read_csv('../{}_fluid_communities.csv'.format(strain)).merge(genes,on='gene_ID',how='left')
	ngroups = df.groupby(by='cluster').ngroups
	for cat in range(ngroups):
		comp = df[df['cluster']==cat]
		k12_ids=[]
		for x in comp['gene_symbol']:
			if not pd.isnull(x) and x in symboldic:
				k12_ids.append(symboldic[x])
		print(k12_ids)
		GOEA.run_GOEA(k12_ids,strain,cat)



