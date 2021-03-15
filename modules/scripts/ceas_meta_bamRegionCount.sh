#!/bin/bash
#calculates the PCR Bottleneck coefficient (PBC) for SAMPLES 
#(i.e. BAMS, i.e. READS)
#ref: https://genome.ucsc.edu/ENCODE/qualityMetrics.html search "PBC"
PVAL=1E-9

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -i|--input) 
    INPUT="$2"
    shift
    ;;
    -p|--pval)
    PVAL="$2" #OPTIONAL ARG b/c we set it up top!
    shift # past argument
    ;;
    -b|--bed) #bed file
    BED="$2"
    shift
    ;;
    -o|--output) 
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

count=$(bedtools intersect -f $PVAL -wa -u -abam $INPUT -b $BED -bed | wc -l)
total=$(samtools flagstat $INPUT | head -1 | cut -d" " -f1)
echo $count,$total > $OUTPUT
