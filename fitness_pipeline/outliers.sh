#!/bin/bash
# outliers.sh experiment_list.csv metadata.csv /outdir 

for i in $(awk -F',' 'FNR>1' $1)
do
	python outliers.py $3/data/bc_counts $3/data/outliers $3/data/kt_errors $3/data/meanfitness $i #> $3/data/outliers/$i.out
done