#!/bin/bash
#subsamples bam to NUM

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -i|--input) 
    INPUT="$2"
    shift
    ;;
    -n|--num)
    NUM="$2" 
    shift # past argument
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

total=$(samtools view -c $INPUT)
frac=$(echo "$NUM / $total" | bc -l)
#CHECK totals
if [ "$frac" '>' 1.0 ]; then
    frac=1.0
fi

samtools view -b -s $frac $INPUT > $OUTPUT
