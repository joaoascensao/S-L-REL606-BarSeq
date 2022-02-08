'''

Process & clean fit2.csv (output of geo.R)

'''

import pandas as pd
import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats import multitest



sns.set(style="whitegrid",font_scale=1.2)

genomes = pd.read_csv('genomes/all_genomes.csv',sep='\t')
rel606 = pd.read_csv('../rel606_genes_wrangled.csv')
fit2 = pd.read_csv('fit2.csv')

tax_priority = {
	199310:3,
	386585:2,
	511145:1,
	155864:4,
}

rel606_gids=[]
rel606_symbols=[]
strains=[]
for i,row in fit2.iterrows():
	success = False
	gid = row['genes.Gene.ID']
	if pd.notnull(gid):
		gid=str(gid)
		if '///' in gid:
			gids0=gid.split('///')
			taxs=[]
			for gid0 in gids0: #arrange by strain
				taxid = genomes[genomes['GeneID']==np.int(gid0)]['tax_id']
				if len(taxid)==1:
					taxid = int(taxid)
					taxs.append(tax_priority[taxid])
				elif len(taxid)==0:
					taxs.append(10)
				else:
					taxid = int(taxid.iloc[0])
					taxs.append(tax_priority[taxid])
			gids1 = [x for _, x in sorted(zip(taxs, gids0))]
		else:
			gids1 = [gid]

		for gid1 in gids1:
			symbol=genomes[genomes['GeneID']==np.int(gid1)]['Symbol']
			if not symbol.empty:
				symbol = str(symbol.iloc[0])
				potrow = rel606[rel606['gene_symbol']==symbol]
				if len(potrow)==1:
					rel606_gids.append(str(potrow['gene_ID'].iloc[0]))
					strains.append(tax_priority[int(genomes[genomes['GeneID']==np.int(gid1)]['tax_id'].iloc[0])])
					rel606_symbols.append(symbol)
					success = True
					break
			aliases = genomes[genomes['GeneID']==np.int(gid1)]['Aliases']
			if not aliases.empty and len(aliases)>1:
				aliases = str(aliases.iloc[0]).split(', ')
				for alias in aliases:
					potrow = rel606[rel606['gene_symbol']==alias]
					if len(potrow)==1:
						rel606_gids.append(str(potrow['gene_ID'].iloc[0]))
						strains.append(tax_priority[int(genomes[genomes['GeneID']==np.int(gid1)]['tax_id'].iloc[0])])
						rel606_symbols.append(alias)
						success = True
						break
				if success:
					break
	if not success:
		rel606_gids.append(np.nan)
		rel606_symbols.append(np.nan)
		strains.append(np.nan)


fit2 = fit2[["coefficients.X6kS.rel606","coefficients.X6kL.rel606","coefficients.X17kL.X6kL","coefficients.X17kS.X6kS","coefficients.X40kL.X6kL","coefficients.X40kS.X6kS","coefficients.X6kS.X6kL","coefficients.X17kS.rel606","coefficients.X17kL.rel606","coefficients.X40kS.rel606","coefficients.X40kL.rel606","stdev.unscaled.X6kS.rel606","stdev.unscaled.X6kL.rel606","stdev.unscaled.X17kL.X6kL","stdev.unscaled.X17kS.X6kS","stdev.unscaled.X40kL.X6kL","stdev.unscaled.X40kS.X6kS","stdev.unscaled.X6kS.X6kL","stdev.unscaled.X17kS.rel606","stdev.unscaled.X17kL.rel606","stdev.unscaled.X40kS.rel606","stdev.unscaled.X40kL.rel606",
"genes.Gene.ID","genes.GO.Function","genes.GO.Process","genes.GO.Component","genes.GO.Function.ID","genes.GO.Process.ID","genes.GO.Component.ID","p.value.X6kS.rel606","p.value.X6kL.rel606","p.value.X17kL.X6kL","p.value.X17kS.X6kS","p.value.X40kL.X6kL","p.value.X40kS.X6kS","p.value.X6kS.X6kL","p.value.X17kS.rel606","p.value.X17kL.rel606","p.value.X40kS.rel606","p.value.X40kL.rel606"]]
fit2['gene_ID']=rel606_gids
fit2['gene_symbol']=rel606_symbols
fit2['taxid']=strains

# adjust p-values for fdr correction
def FDR_correction(pvals,alpha=0.05):
	rejected,pvals_corr,_,_=multitest.multipletests(pvals, alpha=alpha, method='fdr_bh', is_sorted=False)
	return [int(i) for i in rejected], pvals_corr

pvals=["p.value.X6kS.rel606","p.value.X6kL.rel606","p.value.X17kL.X6kL","p.value.X17kS.X6kS","p.value.X40kL.X6kL",
	"p.value.X40kS.X6kS","p.value.X6kS.X6kL","p.value.X17kS.rel606","p.value.X17kL.rel606","p.value.X40kS.rel606","p.value.X40kL.rel606"]
for pval in pvals:
	_,pval_corr=FDR_correction(np.array(fit2[pval]))
	fit2[pval]=pval_corr


fit2 = fit2[pd.notnull(fit2['gene_ID'])]
fit3=[]
for i,group in fit2.groupby(by='gene_ID'):
	if len(group)==1:
		fit3.append(group)
	else:
		row0 = group[group['taxid']==1]
		if len(row0)==1:
			fit3.append(row0)



pd.concat(fit3,ignore_index=True).to_csv('fit3.csv')


