#!/usr/bin/env python
"""Script to parse out the read length fastqc_data.txt (single)
NOTE: just extract that read length and dump to file
"""

import os
import sys
from optparse import OptionParser
import numpy as np

def main():
    usage = "USAGE: %prog -f [path to fastqc_data.txt] -o [path to output file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="path to fastqc_data.txt file")
    optparser.add_option("-o", "--output", help="output file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.file: 
        optparser.print_help()
        sys.exit(-1)

    with open(options.file) as f:
        for l in f:
            if l.strip().startswith("Sequence length"):
                ReadLength = l.strip().split()[2].split('-')
                if len(ReadLength) == 1:
                    outputlenth = ReadLength[0]
                else:
                    intReadLength = [int(i) for i in ReadLength]
                    outputlenth = int(round(np.median(intReadLength)))
    with open(options.output, "w") as o:
        o.write(str(outputlenth))

if __name__=='__main__':
    main()