#!/bin/bash
# xbar_inference.sh experiment_list.csv metadata.csv /outdir 

for i in $(awk -F',' 'FNR>1' $1)
do
	python xbar_inference.py $3/data/bc_counts $3/data/meanfitness $3/data/kt_errors $2 $i
done