#!/usr/bin/env python3
"""GALI BAI
Script to concatenate mapping and pbc information into read level summary table
Prints out a csv file with header- "Sample", "Total(M)", "Mapped(M)", "Uniqely Mapped(M)", "PBC"
"""

import os
import sys
from collections import defaultdict
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -m [input mapping stats] -p [input pbc stats] -o [output csv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-m", "--map", help="mapping.csv file in analysis/align/")
    optparser.add_option("-p", "--pbc", help="pbc.csv file in analysis/frips/")
    optparser.add_option("-o", "--output", help="output .plotly file")
    (options, args) = optparser.parse_args(sys.argv)

    #Setup dictionaries for storing data
    mapping_out=defaultdict(list)
    pbc_out=defaultdict(list)

    ##format mapping stats
    fhd = open(options.map, "rt")
    next(fhd)
    for line in fhd:
        line = line.strip().split(",")
        m = [round(int(line[1])/1000000),round(int(line[2])/1000000),round(int(line[3])/1000000)]
        mapping_out[line[0]].extend(m)
    fhd.close()

    ##format pbc stats
    fhd = open(options.pbc, "rt")
    next(fhd)
    for line in fhd:
        line = line.strip().split(",")
        n = [round(int(line[1])/int(line[2]),2)]
        pbc_out[line[0]].extend(n)
    fhd.close()

    #write to the output file
    out = open(options.output, "w")
    #write header
    out.write("%s\n" % ",".join(["Sample", "Total (M)", "Mapped (M)", "Uniquely Mapped (M)", "PBC"]))
    #concatenate each line by matching sample name
    for k in mapping_out.keys():
        l= [str(k)]
        v= mapping_out[k] + pbc_out[k]
        l.extend(v)
        m = ",".join(map(str,l))
        out.write("%s\n" % m)

    out.close()

if __name__ == '__main__':
    main()
