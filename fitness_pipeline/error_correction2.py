'''
5.19.20
Joao Ascensao

Calculates Levenshtein distances between barcodes, merges barcodes with low edit distance together (<=3)
Only consider merging low abundance barcodes (<3) into high abundance barcodes (>20)
Only merge barcodes that uniquely map at a given distance threshold
Add in heuristics to speed up computation by reducing number of comparisons:
   - 

python error_correction2.py xx.codes output/directory primer SLR[optional]

'''
import pandas as pd
import scipy as sp 
import numpy as np 
import sys
import Levenshtein
import time
from multiprocessing import Pool, TimeoutError
from functools import partial
import swifter

dcut=4

if len(sys.argv)<4:
	raise Exception('Specify all files')

if len(sys.argv)==5:
	mixed=True
else:
	mixed=False

codes=pd.read_csv(sys.argv[1]).set_index('barcode')
outdir=sys.argv[2]
num_primer=sys.argv[3]
primer0='IT0'+num_primer
if int(num_primer)<10:
	primer='IT00'+num_primer
else:
	primer='IT0'+num_primer


minority=codes[codes[primer]<5].reset_index()
majority=codes[codes[primer]>70].reset_index()

dna=['A','T','G','C']

five_mers={
	'a':(1,6),
	'b':(6,11),
	'c':(13,18),
}
five_mers_keys=[key for key in five_mers]
def getsummaries(df):
	df['first']=df['barcode'].apply(lambda x: x[0:4])
	#df['last']=df['barcode'].apply(lambda x: x[-1])
	for neu in dna:
		df[neu]=df['barcode'].apply(lambda x: x.count(neu))
	for key in five_mers:
		b,l=five_mers[key]
		df[key]=df['barcode'].apply(lambda x: x[b:l])
	return df


minority=getsummaries(minority)
majority=getsummaries(majority)

majority['barcode2']=majority['barcode'].copy(deep=True)
minority.set_index('barcode',inplace=True)
majority.set_index('barcode',inplace=True)

print('Starting error correction 2 for primer {}'.format(primer0))

def minLeven(minbc,mt):
	minrow=minority.loc[minbc]
	mtf=mt['first']
	mif=minrow['first']
	cf=pd.eval("mtf==mif")
	mt2 = mt[cf]

	res=mt2['barcode2'].swifter.progress_bar(False).apply(lambda majbc: Levenshtein.distance(minbc,majbc))
	#mt.loc[:,'dist']=doLeven(list(mt['barcode2']),minbc)
	mindist=res.min()
	if mindist<=dcut:
		if len(res[res==mindist])==1:
			indmin=res.idxmin()
			parent=str(mt2.loc[indmin,'barcode2'])

			#print(c,parent,minbc,mindist,np.round(time.time()-start),codes.loc[parent,primer],len(mt))
			#codes.loc[parent,primer] += np.int(minrow[primer])
			#todrop.append(minbc)

			return parent, np.int(minrow[primer]), minbc
	return None, None, None




batches=list(minority.index)

with Pool(processes=4) as pool:
	maj,incr,mi = zip(*pool.map(partial(minLeven,mt=majority),batches))

mil = [mm for mm in mi if pd.notnull(mm)]


codes.drop(inplace=True,index=mil)
for i,m in enumerate(maj):
	if pd.notnull(m):
		codes.loc[m,primer] += incr[i]

if mixed:
	SLR=sys.argv[4]
	codes.to_csv('{}/{}_{}.merged2.codes'.format(outdir,primer0,SLR))
	print('Wrote {}/{}_{}.merged2.codes'.format(outdir,primer0,SLR))
else:
	codes.to_csv('{}/{}.merged2.codes'.format(outdir,primer0))
	print('Wrote {}/{}.merged2.codes'.format(outdir,primer0))
