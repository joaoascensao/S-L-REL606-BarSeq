'''
6.15.19
Joao Ascensao
Uses merged codes/close files along with pool file (from DesignRandomPool.pl)
Creates a table of barcodes with # reads at each time point, insertion position, gene, type of insertion (cds,ncrna,intergenic)
Need to specify n time points in days
Assume barcodes map to 'rcbarcode' in pool file
Also drops barcodes with less than n_cut reads at time 0 
Returns: output/filename.csv, a table of barcode counts per time point, with tnseq information
--
improved speed 6/27/19
--
2.4.20
Modifications to accomodate S/L pooling
--
5.18.20
Generalized

mixed:
python combine_multicode_tn.py n_cut output/dir rel606_genes_wrangled.csv pool/dir dir/error_corrected_codes metadata.csv experiment replicate {S/L/R}

single:
python combine_multicode_tn.py n_cut output/dir rel606_genes_wrangled.csv pool/dir dir/error_corrected_codes metadata.csv experiment,replicate
--
11.30.20

Can now handle technical replicates in first and last data points

'''

import pandas as pd
import scipy as sp 
import numpy as np 
import sys

if len(sys.argv)<8:
	raise Exception('Specify all arguments')

if len(sys.argv)==10:
	mixed=True
	experiment=sys.argv[7]
	replicate=sys.argv[8]
	SLR=sys.argv[9]
elif len(sys.argv)==8:
	mixed=False
	exprep=sys.argv[7].split(',')
	experiment=exprep[0]
	replicate=exprep[1][0]


excludeEnds=0.05#what percent of the gene do I exclude from counting? 5% each end default

n_cut=int(sys.argv[1])

genes=pd.read_csv(sys.argv[3])
pooldir=sys.argv[4]
countsdir=sys.argv[5]#dir/error_corrected_codes

meta=pd.read_csv(sys.argv[6])

cond1=meta['Sample']==experiment
cond2=meta['Replicate']==np.int(replicate)
cond3=meta['Day']==0
mm=meta[cond1 & (cond2 | cond3)]

if mixed:
	print('Combining code counts and Tn pool for {}{}{}'.format(experiment,SLR,replicate))
else:
	print('Combining code counts and Tn pool for {}{}'.format(experiment,replicate))


if mixed:
	poolfile = pooldir + '/' + SLR + '.pool'
else:
	tnpoolname=mm['TnPool'].iloc[0]
	poolfile = pooldir + '/' + tnpoolname + '.pool'
pool=pd.read_csv(poolfile,sep='\t').rename(columns={'barcode':'barcode0'}).rename(columns={'rcbarcode':'barcode'}).set_index('barcode')

#number of day 0 technical reps
mm2=meta[cond1 & cond3]
nd0=len(mm2)

finalday=max(mm['Day']) #final day
ndl = len(mm[mm['Day']==finalday]) #number of final day technical replicates

c=0
ll=1
Rtot={}
for i,row in mm.iterrows():
	num_primer=int(row['Primer'])
	primer0='IT0'+str(int(num_primer))
	if int(num_primer)<10:
		#primer='IT00'+str(int(num_primer))
		primer='IT0'+str(int(num_primer))
	else:
		primer='IT0'+str(int(num_primer))
	if mixed:
		filename=countsdir+'/{}_{}.merged2.codes'.format(primer0,SLR)
	else:
		filename=countsdir+'/{}.merged2.codes'.format(primer0)
	day=np.int(row['Day'])

	if day==0 and nd0>1:
		day='0_'+str(row['Replicate'])
	if day==finalday and ndl>1:
		day=str(finalday)+'_'+str(ll)
		ll+=1

	if c==0:
		counts=pd.read_csv(filename).set_index('barcode')
		counts.rename(columns={primer:day},inplace=True)
	else:
		df=pd.read_csv(filename).set_index('barcode')
		df.rename(columns={primer:day},inplace=True)
		counts=counts.merge(df,how='outer',on='barcode')
	#counts=counts.astype({day:np.int})
	Rtot[day]=np.int(np.sum(counts[day]))
	c+=1


counts.fillna(0,inplace=True)

counts=counts.astype(np.int)

if n_cut>0:
	if nd0==1:
		counts=counts[counts[0]>(n_cut-1)]
	else:
		for j in range(nd0):
			rep0=str(j+1)
			counts=counts[counts['0_'+rep0]>(n_cut-1)]

merged=counts.merge(pool[['strand','pos']],how='left',on='barcode')

merged['gene_ID']=np.nan
merged['type']=np.nan
merged['gene_symbol']=np.nan
merged['description']=np.nan

pos_s=merged['pos']

for index, row in genes.iterrows():
	l=(row['stop']-row['start'])*excludeEnds
	bb=pos_s.between(row['start']+l,row['stop']-l)
	if np.sum(bb)>0:
		merged.loc[bb,'gene_ID']=row['gene_ID']
		merged.loc[bb,'type']=row['type']
		merged.loc[bb,'gene_symbol']=row['gene_symbol']
		merged.loc[bb,'description']=row['description']

if mixed:
	merged.to_csv(sys.argv[2]+'/{}{}{}_counts.csv'.format(experiment,SLR,replicate))
	pd.DataFrame(Rtot,index=[0]).to_csv(sys.argv[2]+'/{}{}{}_Rtot.csv'.format(experiment,SLR,replicate))
	print('Wrote {}{}{}_counts.csv'.format(experiment,SLR,replicate))
else:
	merged.to_csv(sys.argv[2]+'/{}{}_counts.csv'.format(experiment,replicate))
	pd.DataFrame(Rtot,index=[0]).to_csv(sys.argv[2]+'/{}{}_Rtot.csv'.format(experiment,replicate))
	print('Wrote {}{}_counts.csv'.format(experiment,replicate))

