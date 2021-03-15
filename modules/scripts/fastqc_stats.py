#!/usr/bin/env python
"""Script to calc the media Per sequence quality and Per seq GC content
OUTPUT: csv file (to stdout)
"""

import os
import sys
from optparse import OptionParser

def readToTuple(f_path):
    """Reads in a two-col file (tab-delim) and returns a list of tuples"""
    f = open(f_path)
    ls = []
    for l in f:
        if l.startswith("#"):
            continue
        ls.append(tuple(l.strip().split("\t")))
    return ls

def calcMedian(list_o_tuples):
    """Given a list of tuples (A, B), where A = category, and B = counts,
    returns A that represents the median count value"""
    #calc total
    ct = 0
    for (a, b) in list_o_tuples:
        ct += float(b)

    med = ct / 2
    #find A
    ct = 0
    for (i, (a, b)) in enumerate(list_o_tuples):
        ct += float(b)

        if ct > med:
            break

    #print (i, a, b)
    return a


def main():
    usage = "USAGE: %prog -a [_perSeqQual.txt] -b [_perSeqGC.txt]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-a", "--file1", help="path to _perSeqQual.txt file")
    optparser.add_option("-b", "--file2", help="path to _perSeqGC.txt file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.file1 or not options.file2:
        optparser.print_help()
        sys.exit(-1)


    ql = readToTuple(options.file1)
    medQual = calcMedian(ql)
    
    gc = readToTuple(options.file2)
    medGC = calcMedian(gc)

    print("%s,%s" % ("MedianQual", medQual))
    print("%s,%s" % ("MedianGC", medGC))

if __name__=='__main__':
    main()
