#!/bin/bash
#outputs to STDOUT: CHROM\tREADS\t% of total MAPPED reads
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -a|--bam)
    BAM="$2"
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

mappedRow="5q;d" #expect the total mapped reads to be on the 5th line -SED synt

total=$(samtools flagstat $BAM | sed $mappedRow | cut -d" " -f1) 
#samtools idxstats $BAM | awk -v total="$total" '{print $1 "\\t" $3 "\\t" $3/total*100}'
samtools idxstats $BAM | awk -v total="$total" '{printf "%s\t%s\t%.1f\%\n", $1, $3, $3/total*100}'
echo -e "#total_reads\t"$total
