#!/usr/bin/env python
"""Script to parse out the Per sequence quality from fastqc_data.txt (single)
NOTE: just extract that section and dump to file
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [path to fastqc_data.txt] -s [Section to parse out, e.g. Per sequence quality]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="path to fastqc_data.txt file")
    optparser.add_option("-s", "--section", help="Section to parse out, e.g. Per sequence quality")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.file or not options.section:
        optparser.print_help()
        sys.exit(-1)

    f = open(options.file)
    flag = False
    for l in f:
        #look for the section/module
        if l.strip().startswith(">>%s" % options.section):
            flag = True
            continue
        if l.strip().startswith(">>END_MODULE"):
            flag = False

        if flag:
            print(l.strip())

if __name__=='__main__':
    main()


