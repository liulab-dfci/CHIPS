#!/bin/bash
#ref: http://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash 
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -i|--input)
    INPUT="$2"
    shift # past argument
    ;;
    -n)
    NUM="$2"
    shift # past argument
    ;;
    -o|--output)
    OUTPUT="$2"
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

#echo $INPUT $NUM $OUTPUT

count=$(samtools view -Sc $INPUT)
frac=$(echo "$NUM/$count" | bc -l)
#echo $count $frac
if [ "$count" -le "$NUM" ]; then
    #do nothing- just convert sam to bam
    echo "FRiP: read count < $NUM, take ALL reads"
    #samtools view -b $INPUT > $OUTPUT
    cp $INPUT $OUTPUT
else
    echo "FRiP: sampling $NUM reads"
    samtools view -b -s $frac $INPUT> $OUTPUT
fi
