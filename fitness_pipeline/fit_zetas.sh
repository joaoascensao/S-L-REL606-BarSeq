#!/bin/bash
# kt_errors.sh experiment_list.csv metadata.csv /outdir 

for i in $(awk -F',' 'FNR>1' $1)
do
	python fit_betas.py $3/data/kt_errors $3/data/bc_counts $3/data/kt_errors $2 $i
done