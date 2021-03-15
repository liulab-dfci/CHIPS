#!/usr/bin/env python

import os
import sys
from optparse import OptionParser
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json

def main():
    usage = "USAGE: %prog -f [PEAKS FILE] -o [OUTPUT FIGURE]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="_sorted_peaks.narrowPeaks files")
    optparser.add_option("-c", "--ceas", help="_summary.txt files")
    optparser.add_option("-o", "--output", help="output figure path")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.file:
        optparser.print_help()
        sys.exit(-1)

    f = options.file
    #TRY to infer the RUN NAMES
    runID = ".".join(f.strip().split("/")[-1].split('.')[:-1])
    if runID.endswith('_sorted_peaks'):
        runID = runID.replace("_sorted_peaks","")
    f = open(f)
    #start the counts
    tot = fc_10 = fc_20 = 0
    for l in f:
        tmp = l.strip().split("\t") 
        #note FC is 7th col
        fc = float(tmp[6])
        if fc >= 20.0:
            fc_20 += 1
        if fc >= 10.0:
            fc_10 += 1
        tot += 1
    f.close()

    c = options.ceas
    with open(c,"r") as ceas_meta:
        ceas_meta=(ceas_meta.read().replace("\'","\""))
        ceas = json.loads(ceas_meta)
        prom = ceas['Promoter']
        exon = ceas['Exon']
        intr = ceas['Intron']
        inte = ceas['Intergenic']

    data = pd.Series([tot,fc_10,fc_20,prom, exon, intr, inte],
                      index=['Total Peaks','10 Fold Change','20 Fold Change','Promoter', 'Exon', 'Intron', 'Intergenic'])

    sns.set(style="whitegrid",font_scale=1.2)
    f, ax= plt.subplots(figsize = (10, 9))
    ax.set_title('%s Reads in Peaks' % runID)
    sns.barplot(x=data.values,y=data.index,palette="Blues_d",)
    f.savefig(options.output,dpi=100,bbox_inches='tight')     

if __name__=='__main__':
    main()

    