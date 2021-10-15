#! /bin/bash

bed=$1
factor_file=$2
annot=$3
log=$4
tmp=$5
output=$6

total=$(python workflow/scripts/single/find-unique.py $log)
sort -k1,1 -k2,2n $bed > $tmp
factor=$(cat $factor_file)

bedtools intersect -a $annot -b $tmp -sorted -c | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7*'"$factor"'/'"$total"',$7*1000000000/($3-$2)/'"$total"'}' > $output
