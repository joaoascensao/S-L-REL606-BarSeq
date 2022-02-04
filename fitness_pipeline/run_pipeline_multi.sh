#!/bin/bash
# pipeline for BarSeq experiments with two libraries per flask


OUTDIR='../E_VW/'
METADATA='../E_VW/W_meta.csv'
EXPLIST='../E_VW/W_exps.csv'

SEQDIR=''

POOL1='S'
POOL2='L'

sh makedirs.sh $OUTDIR
sh run_MultiCodes.sh $METADATA $OUTDIR $SEQDIR
sh multi_error_correction.sh $METADATA $OUTDIR $POOL1 $POOL2
sh combine_multicode_tn_multi.sh $METADATA $OUTDIR # fix .sh to include $POOL1 and $POOL2 and automatically extract experiment letter

sh kt_errors.sh $EXPLIST $METADATA $OUTDIR
sh fit_betas.sh $EXPLIST $METADATA $OUTDIR
sh xbar_inference.sh $EXPLIST $METADATA $OUTDIR
sh outliers.sh $EXPLIST $METADATA $OUTDIR
sh s_inference.sh $EXPLIST $METADATA $OUTDIR

