#!/usr/bin/env python
"""Script to collect and calcuate the DHS overlap of the top5k peaks
Outputs: 
Run,Total,DHS
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
    out.write("%s\n" % ",".join(["Run","Total","DHS"]))

    for f in options.files:
        #TRY to infer the RUN NAMES
        runID = ".".join(f.strip().split("/")[-1].split('.')[:-1])
        if runID.endswith('_DHS_stats'):
            runID = runID.replace("_DHS_stats","")

        f = open(f)
        total = int(f.readline().strip().split()[0])
        dhs = int(f.readline().strip().split()[0])
        out.write("%s\n" % ",".join([runID,str(total),str(dhs)]))
        f.close()
    out.close()

if __name__=='__main__':
    main()


