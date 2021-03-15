#!/usr/bin/env python

"""Given a set of directories, reads in the available bamRegionStat files, 
i.e. {sample}.exon {sample}.promoter {sample}.DHS and reports the count
"""
import os
import sys
from optparse import OptionParser
import json

def getCount(file_path):
    """returns the first line of the file as a tuple"""
    if os.path.exists(file_path):
        f = open(file_path)
        tmp = tuple(f.readline().strip().split(","))
        f.close()
        return tmp

    return ("--", False)

def main():
    usage = "USAGE: %prog -p [PROMOTER] -d [DHS] -e [EXON] -o [OUTPUT]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-d", "--dhs", help="DHS from bam")
    optparser.add_option("-p", "--prom", help="prom from bam")
    optparser.add_option("-e", "--exon", help="exon from bam")
    optparser.add_option("-o", "--output", help="output json file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.dhs or not options.prom or not options.exon or not options.output:
        optparser.print_help()
        sys.exit(-1)

    DHS = getCount(options.dhs)
    prom = getCount(options.prom)
    exon = getCount(options.exon)
    #UGLY way to get the total:
    if DHS[1]:
        total = DHS[1]
    elif prom[1]:
        total = prom[1]
    elif exon[1]:
        total = exon[1]
    else:
        total = 1
    
    meta_bam = {"Exon": int(exon[0]), "Promoter": int(prom[0]), "DHS": int(DHS[0]), "Total": int(total)}

    with open(options.output,"w") as out:
        json.dump(meta_bam, out)
    

if __name__=='__main__':
    main()
