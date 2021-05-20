#!/usr/bin/env python3
"""GALI BAI
Script to concatenate peak stats, frip score, mapped stats, and dhs stats into peak level summary table
Prints out a csv file with header- Run,Total,10FC,20FC,FRiP,% prom,% exons,% introns,% inter,% DHS
"""

import os
import sys
from collections import defaultdict
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -p [input peaks stats] -f [input frip stats] -m [input peak annotation stats] -d [input dhs stats] -o [output csv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-p", "--peaks", help="peakStats.csv file in analysis/peaks/")
    optparser.add_option("-f", "--frip", help="frips.csv file in analysis/frips/")
    optparser.add_option("-m", "--meta", help="meta.csv file in analysis/ceas/")
    optparser.add_option("-d", "--dhs", help="dhs.csv file in analysis/ceas/")
    optparser.add_option("-s", "--sum", help="output peak level summary file")
    optparser.add_option("-r", "--out_frip", help="formatted frip score summary table")
    optparser.add_option("-a", "--out_annotation", help="formatted peak annotation table")
    optparser.add_option("-o", "--out_dhs", help="formatted dhs summary table")
    (options, args) = optparser.parse_args(sys.argv)

    #Setup dictionaries for storing data
    peaks_out=defaultdict(list)
    frips_out=defaultdict(list)
    meta_out=defaultdict(list)
    meta_nor_out=defaultdict(list)
    dhs_out=defaultdict(list)

    ##format peaks stats
    fhd = open(options.peaks, "rt")
    next(fhd)
    for line in fhd:
        line = line.strip().split(",")
        l = [int(line[1]), int(line[2]), int(line[3])]
        peaks_out[line[0]].extend(l)
    fhd.close()

    ##format frip stats
    fhd = open(options.frip, "rt")
    next(fhd)
    for line in fhd:
        line = line.strip().split(",")
        l = [line[3]]
        frips_out[line[0]].extend(l)
    fhd.close()

    ##format peaks annotations stats
    fhd = open(options.meta, "rt")
    next(fhd)
    for line in fhd:
        line = line.strip().split(",")
        if (int(line[1]) != 0): #Check for div by 0
            r = [round(100*(int(line[2])/int(line[1])), 2), round(100*(int(line[3])/int(line[1])), 2), round(100*(int(line[4])/int(line[1])), 2), round(100*(int(line[5])/int(line[1])) ,2)]
            not_r = [100*(int(line[2])/int(line[1])), 100*(int(line[3])/int(line[1])), 100*(int(line[4])/int(line[1])), 100*(int(line[5])/int(line[1]))]
        else:
            r = [0,0,0,0]
            not_r = [0,0,0,0]
        meta_out[line[0]].extend(r)
        meta_nor_out[line[0]].extend(not_r)
    fhd.close()

    ##format dhs stats
    fhd = open(options.dhs, "rt")
    next(fhd)
    for line in fhd:
        line = line.strip().split(",")
        if (int(line[1]) != 0): #Check for div by 0
            l = [round(100*(int(line[2])/int(line[1])), 2)]
        else:
            l = [0]

        dhs_out[line[0]].extend(l)
    fhd.close()


    #write peaks level summary output file
    out = open(options.sum, "w")
    #write header
    out.write("%s\n" % ",".join(["Run", "Total", "10FC", "20FC", "FRiP", "% prom", "% exons", "% introns", "% inter", "% DHS"]))
    #concatenate each line by matching sample name
    for k in peaks_out.keys():
        l= [str(k)]
        v= peaks_out[k] + frips_out[k] + meta_out[k] + dhs_out[k]
        l.extend(v)
        m = ",".join(map(str,l))
        out.write("%s\n" % m)

    out.close()

    out = open(options.out_frip, "w")
    out.write("%s\n" % ",".join(["Sample", "FRiP"]))
    for k in peaks_out.keys():
        s = [str(k)]
        v = frips_out[k]
        s.extend(v)
        m = ",".join(map(str,s))
        out.write("%s\n" % m)
    out.close()

    out = open(options.out_annotation, "w")
    out.write("%s\n" % ",".join(["Sample", "% peaks in promoters", "% peaks in exons", "% peaks in introns", "% peaks in intergenic regions"]))
    for k in meta_nor_out.keys():
        s = [str(k)]
        v = meta_nor_out[k]
        s.extend(v)
        m = ",".join(map(str,s))
        out.write("%s\n" % m)
    out.close()

    out = open(options.out_dhs, "w")
    out.write("%s\n" % ",".join(["Sample", "% peaks in DHS"]))
    for k in dhs_out.keys():
        s = [str(k)]
        v = dhs_out[k]
        s.extend(v)
        m = ",".join(map(str,s))
        out.write("%s\n" % m)
    out.close()

if __name__ == '__main__':
    main()
