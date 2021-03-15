#!/bin/bash
#calculates the PCR Bottleneck coefficient (PBC) for SAMPLES 
#(i.e. BAMS, i.e. READS)
#ref: https://genome.ucsc.edu/ENCODE/qualityMetrics.html search "PBC"
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -i|--input)
    INPUT="$2"
    shift # past argument
    ;;
    -o|--output) #OUTPUT: summary
    OUTPUT="$2"
    shift
    ;;
    --ignored)
    IGNORED=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

#I know it's one ugly command but...it's ripped from chilin/library/dc.py
#GENERATE the histogram of locations and counts
bedtools bamtobed -i $INPUT | awk '{{l[$1"\\t"$2"\\t"$3"\\t"$6]+=1}} END {{for(i in l) print l[i]}}' | awk '{{n[$1]+=1}} END {{for (i in n) print i"\t"n[i]}}' | sort -k1n - > $OUTPUT

##summarize the HISTOGRAM calculate N1, Nd, and PBC
#pbc=$(awk '{{ if (NR==1) {{N1=$2}} Nd+=$2 }} END {{print N1"\t"Nd"\t"N1/Nd}}' $OUTPUT)
#echo $pbc >> $OUTPUT

