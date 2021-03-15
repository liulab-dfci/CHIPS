#!/usr/bin/env python

import os
import sys
from optparse import OptionParser
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

def main():
    usage = "USAGE: %prog -f [FRIP FILE_1] -o [OUTPUT FIGURE]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="FRIP files")
    optparser.add_option("-o", "--output", help="output figure path")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.file:
        optparser.print_help()
        sys.exit(-1)

    f = options.file
    #UGH: this script is ugly!!
    #TRY to infer the SAMPLE NAMES--SAMPLE.virusseq.ReadsPerGene.out.tab
    sampleID = f.strip().split("/")[-1].split('.')[0]
    #ALSO remove the suffix '_mapping' from the sampleID name
    if sampleID.endswith('.rep1_frip.txt'):
        sampleID = sampleID.replace(".rep1_frip.txt","")

    f = open(f)
    ReadsInPeaks = int(f.readline().strip().split()[1])
    total = int(f.readline().strip().split()[1])
    f.close()
    
    data = pd.Series([total,ReadsInPeaks],
                      index=['Total Reads', 'Reads in Peaks'])

    sns.set(style="whitegrid",font_scale=1.2)
    f, ax= plt.subplots(figsize = (10, 5))
    ax.set_title('%s Reads in Peaks' % sampleID)
    sns.barplot(x=data.values,y=data.index,palette="Blues_d")
    f.savefig(options.output,dpi=100,bbox_inches='tight')       

if __name__=='__main__':
    main()

    