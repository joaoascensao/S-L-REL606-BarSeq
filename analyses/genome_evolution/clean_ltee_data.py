'''

Clean clonal mutation data from Tenaillon et al. (2016), downloaded from https://barricklab.org/shiny/LTEE-Ecoli/

'''

import pandas as pd
import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns



muts = pd.read_csv('LTEE-Ecoli-data.csv')

data = []

def zero_padding(gidn):
	nf=5
	gid = np.str(gidn)
	if len(gid)==5:
		return 'ECB_'+gid
	else:
		return 'ECB_'+(5-len(gid))*'0'+gid

def row2df(row):
	if 'intergenic' in str(row['gene_position']):
		intergenic='intergenic'
	else:
		intergenic='genic'
	return pd.DataFrame({
			'population':row['population'],
			'time':row['time'],
			'strain':row['strain'],
			'clone':row['clone'],
			'mutator_status':row['mutator_status'],
			'type':row['type'],
			'gene_position':row['gene_position'],
			'gene_ID':row['locus_tag'],
			'mutation_category':row['mutation_category'],
			'snp_type':row['snp_type'],
			'html_mutation':row['html_mutation'],
			'start_position':row['start_position'],
			'end_position':row['end_position'],
			'intergenic':intergenic,
		},index=[0])


for i,row in muts.iterrows():
	lt=str(row['locus_tag'])
	if '::' in lt:
		if 't' in lt or 'r' in lt:
			continue
		ids = lt.replace('[','').replace(']','').replace('ECB_','').split('::')
		if not (len(ids[0])==5 and len(ids[1])==5):
			continue
		for gidn in list(range(np.int(ids[0]),np.int(ids[1])+1)):
			row0=row.copy()
			row0['locus_tag']=zero_padding(gidn)
			data.append(row2df(row0))
	else:
		data.append(row2df(row))


				
dff = pd.concat(data,ignore_index=True)
dff['mid'] = ['m'+str(i) for i in range(len(dff))]
dff.to_csv('LTEE-Ecoli-data-clean.csv')