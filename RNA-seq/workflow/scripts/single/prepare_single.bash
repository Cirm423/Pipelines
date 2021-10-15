#! /bin/bash

treat=$1
id_name=$2
factor_file=$3
output=$4

total=`awk -F"\t" '{if(NR!=1) x+=$5} END{print int(x)}' $treat `
factor=$(cat $factor_file)

paste $id_name $treat | awk -F"\t" 'BEGIN{OFS="\t"} {if(NR!=1&&$7!=0) {print $2,$3,$5,$6,$7*'"$factor"'/'"$total"',$8,$9} else if(NR!=1) {print $2,$3,$5,$6,1*'"$factor"'/'"$total"',$8,$9} }' > $output