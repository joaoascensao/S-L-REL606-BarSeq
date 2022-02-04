'''
5.22.20
Joao Ascensao

Inference of mean fitness
Maximum likelihood over chunks then take the median
With ~50 barcodes per chunk, barcodes with 10<counts<300

python xbar_inference.py dir/count dir/out dir/kt_errors metadata.csv experiment,replicate

'''
import pandas as pd
import numpy as np 
import sys


if len(sys.argv)!=6:
	raise Exception('Specify all files')

nchunks=100
minabscount=5

countdir=sys.argv[1]
outdir=sys.argv[2]
ktdir=sys.argv[3]
meta=pd.read_csv(sys.argv[4])


exprep=sys.argv[5].split(',')
experiment=exprep[0]
replicate=exprep[1][0]
counts=pd.read_csv(countdir + '/{}{}_counts.csv'.format(experiment,replicate))
total_counts=pd.read_csv(countdir + '/{}{}_Rtot.csv'.format(experiment,replicate)).drop(columns='Unnamed: 0')
betas=pd.read_csv(ktdir+ '/{}{}_betas.csv'.format(experiment,replicate))
tv = list(total_counts.columns)

d0 = tv[0]

counts=counts[counts[d0]>minabscount]


print('Inferring mean fitness dynamics for {}{}'.format(experiment,replicate))


# get number of day 0 technical reps
if len(experiment)==1:
	cond1=meta['Sample']==experiment
else:
	cond1=meta['Sample']==experiment[:-1]

cond2=meta['Replicate']==np.int(replicate)
cond3=meta['Day']==0
mm=meta[cond1 & (cond2 | cond3)]
mm2=meta[cond1 & cond3]
nd0=len(mm2)

if nd0==2: # if there are two technical replicates for initial time point
	#inverse variance weighted mean to get initial frequency
	f1 = counts['0_1']/np.float(total_counts['0_1'])
	f2 = counts['0_2']/np.float(total_counts['0_2'])
	beta1 = np.float(betas[betas['Day label']=='0_1']['beta + Ne R'])
	beta2 = np.float(betas[betas['Day label']=='0_2']['beta + Ne R'])
	counts['f0'] = (f1/beta1 + f2/beta2)/(1/beta1 + 1/beta2)
elif nd0==1:
	counts['f0'] = counts[tv[0]]/np.float(total_counts[tv[0]])
else:
	raise Exception('More than two d0 technical reps')

# select intergenic barcodes and split into chunks
cond1=pd.isnull(counts['gene_ID'])
cond2=pd.notnull(counts['pos'])
cond3=np.all(counts[tv]<300,axis=1)
counts_int=counts[cond1 & cond2 & cond3]
nc=np.floor(len(counts_int)/nchunks)
chunks=np.array_split(counts_int.sample(frac=1), nc)


days = list(betas[betas['Day']>0]['Day label'])
gens = list(betas[betas['Day']>0]['gens'])
days = [str(day) for day in days]

xbarv=[]
for gene0 in chunks:
	gene = gene0.sum()
	xbar0=[]
	for i,day in enumerate(days):
		xbar00 = np.float( (np.log(gene['f0']) - np.log(gene[day]/total_counts[day]))/gens[i])
		xbar0.append(xbar00)
	xbarv.append(xbar0)

	
xbarv=np.vstack(xbarv)
xbar=np.nanmedian(xbarv,axis=0)
a=np.sum(np.isfinite(xbarv),axis=0)
MAD = np.nanmedian(np.abs(xbarv - xbar),axis=0)/np.sqrt(a)

df_out = pd.DataFrame({
	'Day label':days,
	'xbar':xbar,
	'stderr':MAD/0.67449,
	'gens':gens,
})

df_out.to_csv(outdir+'/{}{}_xbar.csv'.format(experiment,replicate))
print('Wrote'+outdir+'/{}{}_xbar.csv'.format(experiment,replicate))

