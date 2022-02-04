#!/bin/bash
# multi_error_correction.sh metadata.csv /outdir pool1 pool2
for i in $(awk -F',' 'FNR>1 { print $4 }' $1)
do
	python multi_error_correction.py $2/data/MultiCode_tables/IT0${i}.codes $2/data/MultiCode_tables/IT0${i}.close $2/data/errorcorrectedcodes0/ ../TnPools $3 $4
	for j in $3 $4
	do
		python error_correction2.py $2/data/errorcorrectedcodes0/IT0${i}_${j}.merged.codes $2/data/errorcorrectedcodes $i $j
	done
done
