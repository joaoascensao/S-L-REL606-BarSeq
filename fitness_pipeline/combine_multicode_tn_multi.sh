#!/bin/bash
# combine_multicode_tn_multi.sh metadata.csv /outdir
n_cut=0 #minimum number of initial counts for a barcode
for exp in W #experiments
do
	for rep in 1 2 #replicates
	do
		for SLR in S L
		do
			python combine_multicode_tn.py $n_cut $2/data/bc_counts ../rel606_genes_wrangled.csv ../TnPools $2/data/errorcorrectedcodes $1 $exp $rep $SLR
		done
	done
done
