#!/bin/bash

IN_DIR=$1
BED_FILES="ls $1*.bed"
OUT_DIR=$2
LEN=$3

for f in $BED_FILES
do
    fname=$(basename ${f%.bed})
    bedGraphToBigWig $f $LEN $OUT_DIR$fname.bw
done

#to know it's done
touch $OUT_DIR/convert_done.txt
