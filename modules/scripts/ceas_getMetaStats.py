#!/usr/bin/env python
"""Script to collect the peak distribution stats
Outputs: 
Run,Total,Promoter,Exon,Intron,Intergenic
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FILE_1] -f [FILE_2] ...-f [FILE_N] -o [output]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list of _sorted_peaks.narrowPeak files")
    optparser.add_option("-o", "--output", help="output file")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    out = open(options.output,"w")
    out.write("%s\n" % ",".join(["Run","Total","Promoter","Exon","Intron","Intergenic"]))

    for f in options.files:
        #TRY to infer the RUN NAMES
        runID = ".".join(f.strip().split("/")[-1].split('.')[:-1])
        if runID.endswith('_summary'):
            runID = runID.replace("_summary","")

        f = open(f)
        #NOTE: these files are python dictionaries
        tmp = eval(f.readline().strip())
        prom = int(tmp['Promoter'])
        exon = int(tmp['Exon'])
        intron = int(tmp['Intron'])
        intergen = int(tmp['Intergenic'])
        tot = sum([prom,exon,intron,intergen])

        #convert to %

        out.write("%s\n" % ",".join([runID, str(tot), str(prom), str(exon),\
                                         str(intron), str(intergen)]))
        f.close()
    out.close()

if __name__=='__main__':
    main()


