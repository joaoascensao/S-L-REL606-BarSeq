'''
Estimates noise parameters from BarSeq data
   - Only uses intergenic barcodes with initial reads between 40 and 200 count
   - Use MAD to calculate kappa
   - Need to rescale barcode numbers from first day

python kt_errors.py dir/count dir/out metadata.csv experiment,replicate

'''
import pandas as pd
import scipy as sp 
import numpy as np 
import sys
import statsmodels.api as sm
import random


if len(sys.argv)!=5:
	raise Exception('Specify all files')

nboot=500 #number of bootstrapped samples to get kappa std

countdir=sys.argv[1]
outdir=sys.argv[2]
meta=pd.read_csv(sys.argv[3])


exprep=sys.argv[4].split(',')
experiment=exprep[0]
replicate=exprep[1][0]
counts=pd.read_csv(countdir + '/{}{}_counts.csv'.format(experiment,replicate))
Rtot=pd.read_csv(countdir+ '/{}{}_Rtot.csv'.format(experiment,replicate))
print('Calculating noise model error parameters for {}{}'.format(experiment,replicate))

#select intergenic barcodes
cond1=pd.isnull(counts['gene_ID'])
cond2=pd.notnull(counts['pos'])
counts_int=counts[cond1 & cond2]

#select relevant metadata
if len(experiment)==1:
	cond1=meta['Sample']==experiment
else:
	cond1=meta['Sample']==experiment[:-1]

cond2=meta['Replicate']==np.int(replicate)
cond3=meta['Day']==0
mm=meta[cond1 & (cond2 | cond3)]
#number of day 0 technical reps
mm2=meta[cond1 & cond3]
nd0=len(mm2)

finalday=max(mm['Day']) #final day
ndl = len(mm[mm['Day']==finalday]) #number of final day technical replicates

# What are the adjacent days? 
# this is all very inefficient, sorry
daypairs=[]
if nd0>1 and ndl==1:
	days=np.sort(list(set(mm['Day'])))
	for i,d1 in enumerate(days[:-1]):
		for d2 in days[i+1:]:
			if d1==0:
				for rep in range(nd0):
					gen1=0
					gen2=np.float(mm[mm['Day']==d2]['Generations'])
					daypairs.append((str(d1)+'_'+str(rep+1),str(d2),gen1,gen2,int(d2)-int(d1),d1,d2))
			else:
				gen1=np.float(mm[mm['Day']==d1]['Generations'])
				gen2=np.float(mm[mm['Day']==d2]['Generations'])
				daypairs.append((str(d1),str(d2),gen1,gen2,int(d2)-int(d1),d1,d2))
elif nd0>1 and ndl>1:
	days=np.sort(list(set(mm['Day'])))
	for i,d1 in enumerate(days[:-1]):
		for d2 in days[i+1:]:
			if d1==0:
				if d2==finalday:
					continue
				for rep in range(nd0):
					gen1=0
					gen2=np.float(mm[mm['Day']==d2]['Generations'])
					daypairs.append((str(d1)+'_'+str(rep+1),str(d2),gen1,gen2,int(d2)-int(d1),d1,d2))
			elif d2==finalday:
				if d1==0:
					continue
				for rep in range(ndl):
					gen1=np.float(mm[mm['Day']==d1]['Generations'])
					gen2=np.float(mm[mm['Day']==d2]['Generations'].mean()) #works but not the best way to do this
					daypairs.append((str(d1),str(d2)+'_'+str(rep+1),gen1,gen2,int(d2)-int(d1),d1,d2))
			else:
				gen1=np.float(mm[mm['Day']==d1]['Generations'])
				gen2=np.float(mm[mm['Day']==d2]['Generations'])
				daypairs.append((str(d1),str(d2),gen1,gen2,int(d2)-int(d1),d1,d2))
else:
	days=np.sort(mm['Day'])
	for i,d1 in enumerate(days[:-1]):
		for d2 in days[i+1:]:
			gen1=np.float(mm[mm['Day']==d1]['Generations'])
			gen2=np.float(mm[mm['Day']==d2]['Generations'])
			daypairs.append((str(d1),str(d2),gen1,gen2,int(d2)-int(d1),d1,d2))

# Calculate kappa for each day pair
kdata={
	'Day label 1':[],
	'Day label 2':[],
	'Day 1':[],
	'Day 2':[],
	'kappa':[],
	'kappa std':[],
	'R1':[],
	'R2':[],
	'ntransfers':[],
	'gens 1':[],
	'gens 2':[],
}

for dp in daypairs:
	dl1,dl2,gens1,gens2,dt,d1,d2=dp
	R1=np.int(Rtot[dl1])
	R2=np.int(Rtot[dl2])
	cd=counts_int[[dl1,dl2]]
	cond1=cd[dl1]>50 
	cond2=cd[dl1]<500
	cd=cd[cond1 & cond2].sample(frac=1).reset_index(drop=True)

	kva=np.array(np.sqrt(cd[dl1]/R1) - np.sqrt(cd[dl2]/R2))

	kmed = np.median(kva)
	kappa_t = (np.median(np.abs(kva - kmed))/0.67449)**2
	kv=list(kva)
	boot_samples=[]
	for j in range(nboot):
		kv_boot=np.array(random.choices(kv,k=len(kv)))
		kmed_boot = np.median(kv_boot)
		kappa_boot= (np.median(np.abs(kv_boot - kmed_boot))/0.67449)**2
		boot_samples.append(kappa_boot)
	kdata['Day label 1'].append(dl1)
	kdata['Day label 2'].append(dl2)
	kdata['Day 1'].append(d1)
	kdata['Day 2'].append(d2)
	kdata['kappa'].append(kappa_t)
	kdata['kappa std'].append(np.std(boot_samples))
	kdata['R1'].append(R1)
	kdata['R2'].append(R2)
	kdata['ntransfers'].append(dt)
	kdata['gens 1'].append(gens1)
	kdata['gens 2'].append(gens2)


k_df=pd.DataFrame(kdata)
k_df.to_csv(outdir+'/{}{}_kappa.csv'.format(experiment,replicate))
print('Wrote '+outdir+'/{}{}_kappa.csv'.format(experiment,replicate))