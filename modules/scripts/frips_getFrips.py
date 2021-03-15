#!/usr/bin/env python
"""Script to collect and calcuate the frips statistics from across all runs. 
Outputs: 
Run,Total,ReadsInPeaks,FRiP
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FILE_1] -o [output]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="Frip.txt")
    optparser.add_option("-o", "--output", help="output file")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    out = open(options.output,"w")
    out.write("%s\n" % ",".join(["Run","Total","ReadsInPeaks","FRiP"]))

    for f in options.files:
        #TRY to infer the RUN NAMES
        runID = ".".join(f.strip().split("/")[-1].split('.')[:-1])
        if runID.endswith('_frip'):
            runID = runID.replace("_frip","")

        f = open(f)
        readsInPeaks = int(f.readline().strip().split("\t")[1])
        totalReads = int(f.readline().strip().split("\t")[1])
        frip = "%.1f" % (float(readsInPeaks)/totalReads *100.0)

        out.write("%s\n" % ",".join([runID,str(totalReads),str(readsInPeaks),str(frip)]))
        f.close()
    out.close()

if __name__=='__main__':
    main()


