#!/usr/bin/env python

import os
import sys
from optparse import OptionParser
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    usage = "USAGE: %prog -f [FPKM FILE_1] -o [OUTPUT FIGURE]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--file", help="sample_mapping.txt files")
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
    if sampleID.endswith('_mapping'):
        sampleID = sampleID.replace("_mapping","")

    f = open(f)
    total = int(f.readline().strip().split()[0])
    #skip 3 lines
    l = f.readline()
    l = f.readline()
    l = f.readline()
    mapped = int(f.readline().strip().split()[0])
    #skip 8 lines
    l = f.readline()
    l = f.readline()
    l = f.readline()
    l = f.readline()
    l = f.readline()
    l = f.readline()
    l = f.readline()
    l = f.readline()
    uniq_mapped = int(f.readline().strip())
    f.close()
    
    data = pd.Series([total,mapped,uniq_mapped],
                      index=['Total Reads', 'Mapped Reads', 'Uniquely Mapped Reads'])

    sns.set(style="whitegrid",font_scale=1.2)
    f, ax= plt.subplots(figsize = (9, 5))
    ax.set_title('%s Mapped Rate' % sampleID)
    sns.barplot(x=data.values,y=data.index,palette="Blues_d")
    f.savefig(options.output, dpi=100,bbox_inches='tight')       

if __name__=='__main__':
    main()

    