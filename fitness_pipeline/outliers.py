'''
Identify outlier barcodes
1) For intergenic (putatively neutral) barcodes, use a simple likelihood ratio test
2) For barcodes within a gene, use the resistant diagnostic from Rousseeuw (1985)

python s_inference.py dir/count dir/out dir/technical_noise dir/meanfitness experiment,replicate 

'''

import pandas as pd
import scipy as sp 
import numpy as np 
import sys
import time
import scipy.stats
import s_likelihoods as lh
import outlier_functions as of
from itertools import compress
pd.set_option('mode.chained_assignment', None)

if len(sys.argv)!=6:
	raise Exception('Specify all files')

starttime=time.time()
maxtime=100*60*60

min_bc=4
mincount=10
pval_min_itg=0.05
RD_threshold=2

countdir=sys.argv[1]
outdir=sys.argv[2]
betadir=sys.argv[3]
meanfitdir=sys.argv[4]

exprep=sys.argv[5].split(',')
experiment=exprep[0]
replicate=exprep[1][0]
counts=pd.read_csv(countdir + '/{}{}_counts.csv'.format(experiment,replicate)).set_index('barcode',drop=True)
total_counts=pd.read_csv(countdir + '/{}{}_Rtot.csv'.format(experiment,replicate)).drop(columns='Unnamed: 0')
betas=pd.read_csv(betadir+ '/{}{}_betas.csv'.format(experiment,replicate))
xbar=pd.read_csv(meanfitdir+ '/{}{}_xbar.csv'.format(experiment,replicate))




print('Finding outliers for {}{}'.format(experiment,replicate))

tv = list(total_counts.columns)
d0 = tv[0]
dl = tv[-1]

## get intergenic barcodes
cond1=pd.isnull(counts['gene_ID'])
cond2=pd.notnull(counts['pos'])
intergenic_bcs=counts[cond1 & cond2]

shareddata=lh.getshareddata(betas,xbar)

outliers=[]
RDlist=[]

# get outliers of intergenic barcodes by using a likelihood ratio test, compare to chi squared
for bc,row in intergenic_bcs.sample(frac=1).iterrows():
	if row[d0]>mincount:
		d00=lh.getgenedata(row[tv],betas)
		
		if len(d00['r'])<=2:
			outliers.append(bc) #get rid of barcodes that go extinct immediately
			RDlist.append(-1)
			continue
		sn,nll,success = lh.s_maxlikelihood([d00], shareddata, boot=True, boot_bracket=(-0.02,0,0.02))
		if success:
			nll_null = lh.add_likelihoods(0,[d00],shareddata)
			LR=2*(nll_null - nll)
			pval = 1-sp.stats.chi2.cdf(LR,1)
			if pval<pval_min_itg:
				outliers.append(bc)
				RDlist.append(-1)


# detect outliers within genes using the resistant diagnostic
for g,gene in counts.groupby(by=['gene_ID','type','gene_symbol','description']):
	if len(gene)>=min_bc and pd.notnull(g[0]):
		data_list=[]
		for bc,row in gene.iterrows():
			if row[d0]>mincount:
				d00=lh.getgenedata(row,betas)
				if len(d00['r'])<=2:
					#print(row)
					outliers.append(bc) #get rid of barcodes that go extinct immediately
					RDlist.append(-2)
					continue
				s,nll,success = lh.s_maxlikelihood([d00], shareddata, boot=True, boot_bracket=(-0.03,0,0.03))
				if success:
					s_std=lh.std_from_likelihood(s, nll, lambda s: lh.add_likelihoods(s, [d00], shareddata))
					d00['s']=s
					d00['nll']=nll
					d00['bc']=bc
					d00['std']=s_std
					data_list.append(d00)
		if len(data_list)<min_bc:
			continue
		RDs, success = of.resistant_diagnostics(data_list,shareddata)
		if success:
			out2 = [bc for bc in RDs if RDs[bc]>RD_threshold]
			outliers+=out2
			RD2 = [RDs[bc] for bc in RDs if RDs[bc]>RD_threshold]
			RDlist+=RD2
		if (time.time()-starttime)>maxtime:
			break

with open(outdir+'/{}{}_outliers.csv'.format(experiment,replicate), 'w') as f:
    f.write("barcode,RD\n")
    for i,item in enumerate(outliers):
        f.write("{},{}\n".format(item,RDlist[i]))

print('Identified {} outliers'.format(len(outliers)))
print('Wrote outliers to '+outdir+'/{}{}_outliers.csv'.format(experiment,replicate))

