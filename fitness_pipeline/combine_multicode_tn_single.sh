#!/bin/bash
# combine_multicode_tn_single.sh experiment_list.csv metadata.csv /outdir 
n_cut=0 #minimum number of initial counts for a barcode
for i in $(awk -F',' 'FNR>1' $1)
do
	python combine_multicode_tn.py $n_cut $3/data/bc_counts ../rel606_genes_wrangled.csv ../TnPools $3/data/errorcorrectedcodes $2 $i
done

