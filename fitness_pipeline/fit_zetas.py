'''
Using previously calculated kappas, calculate zetas (technical noise), Ne (effective population size), and total variance per time point

python fit_zetas.py dir/kappa dir/count dir/out metadata.csv experiment,replicate

'''
import pandas as pd
import scipy as sp 
import numpy as np 
import scipy.optimize
import sys
import random
import re

if len(sys.argv)!=6:
	raise Exception('Specify all files')

Ne0=1e7 #initial guess for Ne

kappadir=sys.argv[1]
countdir=sys.argv[2]
outdir=sys.argv[3]
meta=pd.read_csv(sys.argv[4])
exprep=sys.argv[5].split(',')
experiment=exprep[0]
replicate=exprep[1][0]
Rtot=pd.read_csv(countdir+ '/{}{}_Rtot.csv'.format(experiment,replicate))
kappa=pd.read_csv(kappadir+ '/{}{}_kappa.csv'.format(experiment,replicate))

#select relevant metadata
if len(experiment)==1:
	cond1=meta['Sample']==experiment
else:
	cond1=meta['Sample']==experiment[:-1]

cond2=meta['Replicate']==np.int(replicate)
cond3=meta['Day']==0
mm=meta[cond1 & (cond2 | cond3)]

print('Fitting beta and Ne for {}{}'.format(experiment,replicate))

# objective function to minimize
def obj_single(kappa,kappa_std,ntransfers,beta1,beta2,Ne):
    kappa_est = beta1 + beta2 + np.float(ntransfers)/(4*Ne)
    return ((kappa_est - np.float(kappa))/np.float(kappa_std))**2

# add up objective functions across all kappas
def add_obj(theta,data,tv):
    mse=0
    Ne=np.exp(theta[-1])
    for i,row in data.iterrows():
        l1=str(row['Day label 1'])
        l2=str(row['Day label 2'])
        if not re.match(r'^-?\d+(?:\.\d+)?$', l1) is None:
            l1=np.str(np.int(np.float(l1)))
        if not re.match(r'^-?\d+(?:\.\d+)?$', l2) is None:
            l2=np.str(np.int(np.float(l2)))
        
        beta1=np.exp(theta[tv.index(l1)])
        beta2=np.exp(theta[tv.index(l2)])
        mse+=obj_single(row['kappa'],row['kappa std'],row['ntransfers'],beta1,beta2,Ne)
    return mse



tv=list(kappa['Day label 1']) + list(kappa['Day label 2'])
tv=list(set([str(t) for t in tv]))

bounds=[]
theta0=[]
for i,dl in enumerate(tv):
    R = np.float(Rtot[dl])
    bounds.append((np.log(1/(4*R)),np.log(100/(4*R))))
    theta0.append(np.log(5/(4*R)))

bounds.append((np.log(1e3),np.log(1e9)))
theta0.append(np.log(Ne0))

# minimize objective function, get zetas and Ne
res=sp.optimize.minimize(lambda theta: add_obj(theta,kappa,tv),theta0, method='L-BFGS-B',bounds=bounds)

theta = np.exp(res.x)
Ne = theta[-1]

bdata={
	'Day label':[],
	'Day':[],
    'beta':[],
	'beta + Ne R':[],
	'R':[],
	'gens':[],
    'Ne':[],
}

for i,dl in enumerate(tv):
    R = np.float(Rtot[dl])
    beta = theta[i]
    if re.match(r'^-?\d+(?:\.\d+)?$', dl) is None:
        day = np.int(np.float(dl.split('_')[0]))
    else:
        day = np.int(np.float(dl))
    
    gens=np.float(mm[mm['Day']==day]['Generations'].mean())

    # original name for zeta was beta
    bdata['Day label'].append(dl)
    bdata['Day'].append(day)
    bdata['beta'].append(beta)
    bdata['beta + Ne R'].append((4*beta + (1/Ne))*R)
    bdata['R'].append(R)
    bdata['gens'].append(gens)
    bdata['Ne'].append(Ne)

b_df=pd.DataFrame(bdata).sort_values('Day')
b_df.to_csv(outdir+'/{}{}_betas.csv'.format(experiment,replicate))
print('Wrote '+outdir+'/{}{}_betas.csv'.format(experiment,replicate))