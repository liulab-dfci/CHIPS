#!/usr/bin/env python

"""Given a set of summary files, this combines them into a single file 
suitable for the report
"""
import os
import sys
from optparse import OptionParser

def humanReadable(n):
    """Given an int, e.g. 52071886 returns a human readable string 5.2M
    ref: http://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size
    """
    for unit in ['','K','M','B','T','P','E','Z']:
        if abs(n) < 1000.0:
            return "%3.1f%s" % (n, unit)
        n /= 1000.0
    return "%.1f%s" % (num, 'Y')
    
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
    usage = "USAGE: %prog -p [peakStats.csv] -f [frips.csv] ... -o [output]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-p", "--peaks", help="peakStats.csv file")
    optparser.add_option("-f", "--frips", help="frips.csv file")
    optparser.add_option("-d", "--dhs", help="dhs.csv file")
    optparser.add_option("-m", "--meta", help="meta.csv file")
    optparser.add_option("-o", "--output", help="output file")

    (options, args) = optparser.parse_args(sys.argv)

    if (not options.peaks or not options.frips or not options.output):
        optparser.print_help()
        sys.exit(-1)

    tmp = parseCSV(options.peaks)
    runs = sorted(list(tmp.keys()))

    #initialize stats table--set keys to samples
    stats = {runID:{} for runID in runs}
    
    #get MedianQuality and set it as FastQC
    addStat(stats, tmp, 'Total', 'TotalPeaks')
    addStat(stats, tmp, '10FC', 'FC>10')
    addStat(stats, tmp, '20FC', 'FC>20')

    #HANDLE mapping.csv
    tmp = parseCSV(options.frips)
    addStat(stats, tmp, 'FRiP', 'FRiP')
    #addStat(stats, tmp, 'Mapped', 'MappedReads', True)

    #HANDLE DHS
    tmp = parseCSV(options.dhs)
    for r in runs:
        #HANDLE zero-counts:
        dhs = float(tmp[r]['DHS'])
        tot = int(tmp[r]['Total'])
        ratio = dhs/tot if tot else 'NA'
        stats[r]['DHS_peaks'] = "%s" % int(dhs)
        stats[r]['DHS_%'] = "%.2f" % (ratio*100) if ratio != 'NA' else "NA"

    #HANDLE META
    tmp = parseCSV(options.meta)
    for r in runs:
        tot = int(tmp[r]['Total'])
        stats[r]['Promoter'] = "%.1f" % (float(tmp[r]['Promoter'])/tot *100) if tot else 'NA'
        stats[r]['Exon'] = "%.1f" % (float(tmp[r]['Exon'])/tot *100) if tot else 'NA'
        stats[r]['Intron'] = "%.1f" % (float(tmp[r]['Intron'])/tot *100) if tot else 'NA'
        stats[r]['Intergenic'] = "%.1f" % (float(tmp[r]['Intergenic'])/tot *100) if tot else 'NA'

    #OUTPUT- fields defines the column order
    fields = ['TotalPeaks', 'FC>10', 'FC>20', 'FRiP','DHS_peaks', 'DHS_%', 
              'Promoter','Exon','Intron','Intergenic']
    out = open(options.output,"w")
    out.write("%s\n" % ",".join(['Run'] + fields))

    for s in runs:
        vals = [stats[s][f] for f in fields]
        out.write("%s\n" % ",".join([s] + vals))

if __name__=='__main__':
    main()
