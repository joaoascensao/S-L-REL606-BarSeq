#!/bin/bash
# sh makedirs.sh outdir
outdir=$1

mkdir -p $outdir
mkdir -p $outdir/data
mkdir -p $outdir/data/MultiCode_tables/
mkdir -p $outdir/data/errorcorrectedcodes0/
mkdir -p $outdir/data/errorcorrectedcodes/
mkdir -p $outdir/data/bc_counts/
rm -r $outdir/data/kt_errors/
rm -r $outdir/data/meanfitness/
mkdir -p $outdir/data/meanfitness/
mkdir -p $outdir/data/kt_errors/
mkdir -p $outdir/data/outliers/
mkdir -p $outdir/data/fitness/


echo 'Made directories in '$1