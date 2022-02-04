'''
6.30.19
Joao Ascensao

Converts REL606 gene IDs to MG1655 IDs using uniprot's ecoli.txt https://www.uniprot.org/docs/ecoli

'''

import pandas as pd
import scipy as sp 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages
sns.set(style="whitegrid")

dire='./'

def uniprot_conv():
	conv=pd.read_csv(dire+'ecoli.txt')
	conv['b'], conv['JW'] = conv['Ordered locus name'].str.split(';', 1).str
	conv.set_index('b',inplace=True)


	symbollist=conv['symbols'].str.split(';').tolist()
	indices=conv.index
	symboldic={}


	for i,symbols in enumerate(symbollist):
		if isinstance(symbols, list):
			for symbol in symbols:
				symboldic.update({symbol:indices[i]})
	return symboldic

def ncbi_conv():
	#b ordered locus name : ncbi
	ncbi=pd.read_csv(dire+'gene_result_NCBI.txt',sep='\t').set_index('GeneID')


	aliaslist=ncbi['Aliases'].str.split(', ').tolist()
	indices=ncbi.index
	aliasesdic={}

	for i,aliases in enumerate(aliaslist):
		if isinstance(aliases, list):
			for alias in aliases:
				aliasesdic.update({alias:indices[i]})
	return aliasesdic

def uniprot_to_ncbi():
	uni=uniprot_conv()
	ncbi=ncbi_conv()
	return {symbol:ncbi[loc] for symbol, loc in uni.items() if loc in ncbi}
		

def ncbi_to_rel606ID(ids,name):
	rel606=pd.read_csv(dire+'rel606_genes_wrangled.csv')
	rel606=rel606[rel606['type']=='cds']
	rel606=rel606.set_index('gene_symbol')
	if name=='Early LTEE mutations':
		return [rel606['gene_ID'][id_] for id_ in ids if id_ in rel606['gene_ID']]

	aliasesdic=uniprot_to_ncbi()
	outp=[]
	for id_ in ids:
		for symbol,ncbi in aliasesdic.items():
			if id_==ncbi and symbol in rel606.index:
				outp.append(rel606['gene_ID'][symbol])
				#outp.append(symbol)
				break
	return outp


