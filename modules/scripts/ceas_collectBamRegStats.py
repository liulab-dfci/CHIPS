#!/usr/bin/env python

"""Given a set of directories, reads in the available bamRegionStat files, 
i.e. {sample}.exon {sample}.promoter {sample}.DHS and reports the count
"""
import os
import sys
from optparse import OptionParser

def getCount(file_path):
    """returns the first line of the file as a tuple"""
    if os.path.exists(file_path):
        f = open(file_path)
        tmp = tuple(f.readline().strip().split(","))
        f.close()
        return tmp

    return ("--", False)

def main():
    usage = "USAGE: %prog -d [directory_1] -d [dir_2] ...-d [dir_N]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-d", "--dirs", action="append", help="list of bam region dirs")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.dirs:
        optparser.print_help()
        sys.exit(-1)

    print(",".join(["Sample","DHS","Promoter","Exon","Total"]))

    for d in options.dirs:
        if os.path.isdir(d):
            sampleID = d.strip().split("/")[-1]
            DHS = getCount(os.path.join(d, "%s.DHS" % sampleID))
            prom = getCount(os.path.join(d, "%s.promoters" % sampleID))
            exon = getCount(os.path.join(d, "%s.exons" % sampleID))

            #UGLY way to get the total:
            if DHS[1]:
                total = DHS[1]
            elif prom[1]:
                total = prom[1]
            elif exon[1]:
                total = exon[1]
            else:
                total = '1'

            print(",".join([sampleID,DHS[0],prom[0],exon[0],total]))

if __name__=='__main__':
    main()
