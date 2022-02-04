#!/bin/bash
# run_MultiCodes.sh metadata.csv /outdir /rawdatadir
for i in $(awk -F',' 'FNR>1 { print $4 }' $1)
do
	ip=$(printf '%02d' $i)
	x=`find $3 -name "*IT0${ip}_*.fastq"  -type f`
	perl feba_pipeline/bin/MultiCodes.pl -out $2/data/MultiCode_tables/IT0${i} -index IT0${ip} -bs3 < $x
done