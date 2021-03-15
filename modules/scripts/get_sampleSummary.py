#!/usr/bin/env python
"""Generates the sampleSummary.csv file"""
import os
import sys
from optparse import OptionParser

def humanReadable(n):
    """Given an int, e.g. 52071886 returns a human readable string 5.2M
    ref: http://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size
    """
    #DROPPING the human readable part and just moving to millions
    # for unit in ['','K','M','B','T','P','E','Z']:
    #     if abs(n) < 1000.0:
    #         return "%3.1f%s" % (n, unit)
    #     n /= 1000.0
    # return "%.1f%s" % (num, 'Y')

    return "%.1f" % (n/1000000.0)
    
    
def parseCSV(csv_file):
    """parses a CSV file, using the first line as a header (of column names)
    ASSUMES: sampleIDs are the first column
    returns a dictionary of dictionarys, e.g.
    {sampleID: {col1: val1, col2: val2, ...}}
    """
    ret = {}
    f = open(csv_file)
    cols = f.readline().strip().split(",")
    
    #read the rest
    for l in f:
        tmp = l.strip().split(',')
        #enter new record - don't enter the first col, which is sampleID
        ret[tmp[0]] = dict(zip(cols[1:], tmp[1:]))
    f.close()
    return ret

def addStat(d, csv, fieldToGet, fieldToStore, readable=False):
    """ADD fieldToGet (value from csv) to d AS fieldToStore
    this fn has side-effects
    readable = flag to make an number humanReadable (see fn above)
    """
    for s in list(d.keys()):
        if s in csv:
            val = humanReadable(int(csv[s][fieldToGet])) if readable \
                else csv[s][fieldToGet]
            d[s][fieldToStore] = val
        else:
            d[s][fieldToStore] = "NA"
    

def main():
    usage = "USAGE: %prog -f [fastqc.csv] -m [mapping.csv] -p [pbc.csv]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--fastqc", help="fastqc.csv file")
    optparser.add_option("-m", "--mapping", help="mapping.csv file")
    optparser.add_option("-p", "--pbc", help="pbc.csv file")
    (options, args) = optparser.parse_args(sys.argv)

    if (not options.fastqc or not options.mapping or not options.pbc):
        optparser.print_help()
        sys.exit(-1)

    #print(",".join(["Sample","EstFragment","StdDev"]))

    #handle fastqc - get MedianQuality, store as FastQC
    tmp = parseCSV(options.fastqc)
    samples = sorted(list(tmp.keys()))

    #initialize stats table--set keys to samples
    stats = {sampleID:{} for sampleID in samples}
    
    #get MedianQuality and set it as FastQC
    addStat(stats, tmp, 'MedianQuality', 'FastQC')

    #HANDLE mapping.csv
    tmp = parseCSV(options.mapping)
    addStat(stats, tmp, 'Total', 'TotalReads(M)', True)
    addStat(stats, tmp, 'Mapped', 'MappedReads(M)', True)
    addStat(stats, tmp, 'UniquelyMapped', 'UniqMappedReads(M)', True)

    #HANDLE pbc.csv
    tmp = parseCSV(options.pbc)
    addStat(stats, tmp, 'Nd', 'UniqLoc4M')
    addStat(stats, tmp, 'N1', 'UniqLoc1read4M')
    #PBC = N1/Nd
    for s in samples:
        stats[s]['PBC'] = \
            int(stats[s]['UniqLoc1read4M']) / float(stats[s]['UniqLoc4M'])*100
        #FORMAT: convert to 2-decimal percentage
        stats[s]['PBC'] = "%.2f" % stats[s]['PBC']

    #MAKE UniqLoc4M and UniqLoc1read4M  human readable
    for s in samples:
        stats[s]['UniqLoc4M'] = humanReadable(int(stats[s]['UniqLoc4M']))
        stats[s]['UniqLoc1read4M'] = humanReadable(int(stats[s]['UniqLoc1read4M']))

    #print(stats)

    #OUTPUT- fields defines the column order
    fields = ['FastQC', 'TotalReads(M)', 'MappedReads(M)','UniqMappedReads(M)',
              'UniqLoc4M', 'UniqLoc1read4M', 'PBC']
    print(",".join(['Sample'] + fields))

    for s in samples:
        vals = [stats[s][f] for f in fields]
        print(",".join([s] + vals))

if __name__=='__main__':
    main()
