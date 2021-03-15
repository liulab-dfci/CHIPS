#!/bin/bash
#A script that will take a HOMER (motif program) annotatePeaks.pl output
#and generate two output files--a TSV and CSV equivalent file with:
#PeakID (....) simplified to just PeakID

Input=$1
Out_tsv=$2
Out_csv=$3

#Create header--THIS doesn't work! it's not TSV when it's saved as a bash-var
#Header=`head -n 1 $1 | awk -F"\t" '{print "PeakID\t" substr($0, index($0, $2))}'`
#echo $Header

#TSV output
#echo $Header > $Out_tsv
head -n 1 $Input | awk -F"\t" '{print "PeakID\t" substr($0, index($0, $2))}' > $Out_tsv
tail -n +2 $Input >> $Out_tsv

#CSV output
cat $Out_tsv | tr "\\t" "," > $Out_csv
