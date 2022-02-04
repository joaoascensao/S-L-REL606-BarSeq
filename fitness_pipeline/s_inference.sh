#!/bin/bash

# s_inference.sh experiment_list.csv metadata.csv /outdir 

for i in $(awk -F',' 'FNR>1' $1)
do
	python s_inference.py $3/data/bc_counts $3/data/fitness $3/data/kt_errors $3/data/outliers $3/data/meanfitness $i 80 #> $3/data/fitness/$i.out
done

python s_inference_combo.py $3/data/bc_counts $3/data/fitness $3/data/kt_errors $3/data/outliers $3/data/fitness $3/data/meanfitness $1 80
