#!/bin/bash

BDG_FILE=$1
OUT_DIR=$2

for chr in `cut -f 1 $BDG_FILE | sort | uniq`; do
    echo $chr >> $OUT_DIR/chroms.txt   
    grep -w $chr $BDG_FILE > $OUT_DIR/$chr.bed
done


#to know it's done
touch $OUT_DIR/split_done.txt
