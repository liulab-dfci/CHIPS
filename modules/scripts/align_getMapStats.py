#!/usr/bin/env python
"""Script to collect the mapping statistics from across all samples. 
outputs to stdout:
Sample,Mapped,Total,Percentage
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FPKM FILE_1] -f [FPKM FILE_2] ...-f [FPKM FILE_N]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list of sample_mapping.txt files (note: these are snakemake temp files)")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.files:
        optparser.print_help()
        sys.exit(-1)

    print(",".join(["Sample","Total","Mapped","UniquelyMapped"]))

    for f in options.files:
        #UGH: this script is ugly!!
        #TRY to infer the SAMPLE NAMES--SAMPLE.virusseq.ReadsPerGene.out.tab
        sampleID = f.strip().split("/")[-1].split('.')[0]
        #ALSO remove the suffix '_mapping' from the sampleID name
        if sampleID.endswith('_mapping'):
            sampleID = sampleID.replace("_mapping","")

        f = open(f)
        total = int(f.readline().strip().split()[0])
        #skip 3 lines
        l = f.readline()
        l = f.readline()
        l = f.readline()
        mapped = int(f.readline().strip().split()[0])
        #skip 8 lines
        l = f.readline()
        l = f.readline()
        l = f.readline()
        l = f.readline()
        l = f.readline()
        l = f.readline()
        l = f.readline()
        l = f.readline()
        uniq_mapped = int(f.readline().strip())
        print(",".join([sampleID,str(total),str(mapped),str(uniq_mapped)]))#"%.2f" % (float(mapped)/total *100)]))

if __name__=='__main__':
    main()


