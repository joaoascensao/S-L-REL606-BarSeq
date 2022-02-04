#!/bin/bash
# pipeline for BarSeq experiments with only one library per flask

OUTDIR='../E_VW/'
METADATA='../E_VW/V_meta.csv'
EXPLIST='../E_VW/V_exps.csv'

SEQDIR=''


sh makedirs.sh $OUTDIR
sh run_MultiCodes.sh $METADATA $OUTDIR $SEQDIR
sh single_error_correction.sh $METADATA $OUTDIR
sh combine_multicode_tn_single.sh $EXPLIST $METADATA $OUTDIR

sh kt_errors.sh $EXPLIST $METADATA $OUTDIR
sh fit_betas.sh $EXPLIST $METADATA $OUTDIR
sh xbar_inference.sh $EXPLIST $METADATA $OUTDIR
sh outliers.sh $EXPLIST $METADATA $OUTDIR
sh s_inference.sh $EXPLIST $METADATA $OUTDIR
