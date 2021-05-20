#!/usr/bin/env python3
"""GALI BAI
Script to generate coordinates extended and tss only bed files using refGene.bed.
Prints out two bed files: extend.bed and tss.bed
"""

import os
import sys
from collections import defaultdict
from optparse import OptionParser
import pandas as pd
import numpy as np

def main():
    usage = "USAGE: %prog -i [assembly_refGene.bed] -u [extend upstream] -d [extend downsream] -e [coord extended bed file] -t [tss bed file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-i", "--input", help="refGene.bed file in reference folder")
    optparser.add_option("-u", "--upstream", type="int", help="extend the left border of gene coordinates to this bp")
    optparser.add_option("-d", "--downstream", type="int", help="extend the right border of gene coordinates to this bp")
    optparser.add_option("-e", "--extend", help="coordinates extended refGene.bed file")
    optparser.add_option("-t", "--tss", help="bed file with only tss coordinates of refGene.bed")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.input:
        optparser.print_help()
        sys.exit(-1)

    def pygenomeTacks_bed(ref_bed, extend_bed, tss_bed, upstream, downstream):
        #deduplicate based on gene symbol
        df_bed= pd.read_csv(ref_bed, sep = '\t',header=0, index_col=False).drop_duplicates(subset=['symbol'], keep='first')
        #extend the region of GenomeTrack view based on users' input
        df_bed["left_extend"] = df_bed["start"] - upstream
        df_bed["right_extend"] = df_bed["end"] + downstream
        extend_list = []
        tss_list = []
        for i in df_bed.index.tolist():
            chromosome = df_bed.loc[i, "chromosome"]
            start = df_bed.loc[i, "start"]
            end = df_bed.loc[i, "end"]
            if df_bed.loc[i, "left_extend"] < 0:
                istart = 0
            else:
                istart = df_bed.loc[i, "left_extend"]
            iend = df_bed.loc[i, "right_extend"]
            coordinate = '%s:%s-%s' % (chromosome, istart, iend)
            product_accession = df_bed.loc[i, "product_accession"]
            strand = df_bed.loc[i,"strand"]
            symbol = df_bed.loc[i, "symbol"]
            TSS = df_bed.loc[i, "TSS"]
            try:
                TSS_end = df_bed.loc[(i+1),"TSS"]
            except:
                TSS_end = df_bed.loc[i, "TSS"] + 1
            extend_list.append([chromosome, str(start), str(end), str(symbol), str(coordinate), str(product_accession), str(strand), str(TSS)])
            if TSS < TSS_end:
                tss_list.append([chromosome, str(TSS), str(TSS_end), str(coordinate), str(product_accession), str(strand), str(symbol), str(TSS)])
        print(extend_list[0:10])
        print(tss_list[0:10])
        f1 = open(extend_bed, "w")
        f1.write("\n".join(list(map(lambda x: "\t".join(x), extend_list))))
        f1.close()
        f2 = open(tss_bed, "w")
        f2.write("\n".join(list(map(lambda x: "\t".join(x), tss_list))))
        f2.close()

    pygenomeTacks_bed(options.input, options.extend, options.tss, options.upstream, options.downstream)

if __name__ == '__main__':
    main()
