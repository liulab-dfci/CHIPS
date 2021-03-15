#!/bin/bash
#calculates FRiP, outputs # reads under peak and # total reads to STOUD
#ref: http://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash 
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -a|--bam)
    BAM="$2"
    shift # past argument
    ;;
    -b|--bed)
    BED="$2"
    shift # past argument
    ;;
    -p)
    PVAL="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

#echo $BAM $BED $PVAL
fr=$(bedtools intersect -f $PVAL -wa -u -abam $BAM -b $BED -bed | wc -l)
total=$(samtools flagstat $BAM | head -1 | cut -d" " -f1) 
echo -e "#reads_under_peaks\t"$fr
echo -e "#total_reads\t"$total
