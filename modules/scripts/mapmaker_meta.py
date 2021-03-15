#!/usr/bin/env python
"""Script to generate a mapmaker metasheet based on the chips run names
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FILE_1] -f [FILE_2] ...-f [FILE_N] -o [output]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-n", "--names", action="append", help="run names")
    optparser.add_option("-o", "--output", help="output file")

    (options, args) = optparser.parse_args(sys.argv)

    #if not options.files or not options.output:
    if not options.names or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #OUTPUT
    out = open(options.output,"w")
    #header
    out.write("Sample,Type,comp_B_vs_A\n")
    for n in options.names:
        out.write("%s\n" % ",".join([n,"",""]))
    out.close()

if __name__=='__main__':
    main()


