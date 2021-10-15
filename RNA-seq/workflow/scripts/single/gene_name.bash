#! /bin/bash

ref=$1
out=$2

echo -e "Gene_ID\tGene_Name" > $out
awk '{if($3=="gene") {split($10,a,"\"");split($14,b,"\"");print a[2]"\t"b[2]}}' $ref | \
sort -k1,1 >> $out