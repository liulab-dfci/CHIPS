#!/usr/bin/env python3
"""GALI BAI
Script to calculate PBC by dividing N1/Nd
Prints out a csv file with header- "Sample", "PBC"
"""

import os
import sys
from collections import defaultdict
from optparse import OptionParser
import numpy as np

def main():
    usage = "USAGE: %prog -p [input pbc stats] -o [output csv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-p", "--pbc", help="pbc.csv file in analysis/frips/")
    optparser.add_option("-o", "--output", help="output .plotly file")
    (options, args) = optparser.parse_args(sys.argv)

    #Setup dictionaries for storing data
    pbc_out=defaultdict(list)

    ##format pbc stats
    fhd = open(options.pbc, "rt")
    next(fhd)
    for line in fhd:
        line = line.strip().split(",")
        n = [int(line[1])/int(line[2])]
        pbc_out[line[0]].extend(n)
    fhd.close()

    #write to the output file
    out = open(options.output, "w")
    out.write("%s\n" % ",".join(["Sample", "PBC"]))
    for k in pbc_out.keys():
        l= [str(k)]
        l.extend(pbc_out[k])
        m = ",".join(map(str,l))
        out.write("%s\n" % m)

    out.close()

if __name__ == '__main__':
    main()
