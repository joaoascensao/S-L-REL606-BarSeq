'''
6.14.19
Joao Ascensao
Merges codes and close files (outputs of MultiCodes.pl)
For off-by-one pairs, if there are n>5 type 1 barcodes and n<5 type 2 barcodes,
merge type 2 barcodes into type 1
Usage: python merge_counts_close.py xx.codes xx.close output/directory
Returns a modified codes file, xx.merged.codes

---
2.4.20
Modify to handle S&L libraries pooled together in same experiment
1. Remove BCs found in both libraries
2. Sort close by minority/majority. Group the minority barcodes. Throw away all that map to >1 majority bc
3. Merge close/counts
4. Split into S and L

---
5.18.20
Generalized to allow for any combination of two libraries mixed together. Provide both the pool directory location and pool names

Usage: python merge_counts_close.py xx.codes xx.close output/directory pool/directory pool1 pool2
'''
import pandas as pd
import scipy as sp 
import numpy as np 
import sys

if len(sys.argv)!=7:
	print(sys.argv)
	raise Exception('Specify codes and close files')


codes=pd.read_csv(sys.argv[1],sep='\t').set_index('barcode')
close=pd.read_csv(sys.argv[2],sep='\t')


primer=sys.argv[1].split('.')[-2].split('/')[-1]

print('Starting error correction 1 for primer {}'.format(primer))

genpool=lambda xx: pd.read_csv(xx,sep='\t').rename(columns={'barcode':'barcode0'}).rename(columns={'rcbarcode':'barcode'}).set_index('barcode')
pooldir = sys.argv[4]
name1=sys.argv[5]
name2=sys.argv[6]
pool1=genpool(pooldir + '/' + name1 + '.pool')
pool2=genpool(pooldir + '/' + name2 + '.pool')

#barcodes in both the S and L pools
shared_barcodes=pool1.merge(pool2,how='inner',left_index=True,right_index=True).index

#get rid of shared barcodes in codes
codes.drop(index=shared_barcodes,errors='ignore',inplace=True)

#make new data frame for close for majority/minority barcodes. fill iteratively

cc=[]
for i,row in close.iterrows():
	if row['count1']>5 and row['count2']<5:
		data={
		'majoritycode':row['code1'],
		'minoritycode':row['code2'],
		'majoritycount':row['count1'],
		'minoritycount':row['count2'],
		}
	elif row['count1']<5 and row['count2']>5:
		data={
		'majoritycode':row['code2'],
		'minoritycode':row['code1'],
		'majoritycount':row['count2'],
		'minoritycount':row['count1'],
		}
	else:
		continue
	cc.append(pd.DataFrame(data,index=[0]))

close2=pd.concat(cc,ignore_index=True)


#drop rows where minority barcodes that map to more than 1 majority barcode
gc=close2.groupby(by='minoritycode')
close_filter=gc.filter(lambda x: len(x)<2)

#merge minority into majority (in codes), drop minority bcs
close_filter_sum=close_filter[['majoritycode','minoritycount']].groupby(by='majoritycode').sum()

primer2=list(codes)[0]
codes.drop(index=close_filter['minoritycode'],inplace=True,errors='ignore')

mergecodes=codes.merge(close_filter_sum,how='left',left_index=True,right_index=True).fillna(0)
mergecodes['count']=mergecodes[primer2]+mergecodes['minoritycount']


counts1=pool1.merge(mergecodes,how='inner',left_index=True,right_index=True)[['count']].rename(columns={'count':primer})
counts2=pool2.merge(mergecodes,how='inner',left_index=True,right_index=True)[['count']].rename(columns={'count':primer})

counts1.to_csv('{}/{}_{}.merged.codes'.format(sys.argv[3],primer,name1))
counts2.to_csv('{}/{}_{}.merged.codes'.format(sys.argv[3],primer,name2))
print('Wrote {}/{}.merged.codes'.format(sys.argv[3],primer))