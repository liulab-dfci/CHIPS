#!/usr/bin/env python
"""Script to collect the fastqc statistics from across all samples. 
outputs to stdout:
Sample,MedianQuality,MedianGC
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FILE_1] -f [FILE_2] ...-f [FILE_N]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list of fastqc/{sample}.csv files")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.files:
        optparser.print_help()
        sys.exit(-1)

    print(",".join(["Sample","MedianQuality","MedianGC"]))

    for f in sorted(options.files):
        #UGH: this script is ugly!!
        #TRY to infer the SAMPLE NAMES--SAMPLE.virusseq.ReadsPerGene.out.tab
        sampleID = f.strip().split("/")[-1].split('.')[0]
        if sampleID.endswith('_stats'):
            sampleID = sampleID.replace("_stats","")

        f = open(f)
        #first line quality; second GC
        qual = int(f.readline().strip().split(",")[1])
        gc = int(f.readline().strip().split(",")[1])
        print(",".join([sampleID,str(qual),str(gc)]))
        f.close()

if __name__=='__main__':
    main()


