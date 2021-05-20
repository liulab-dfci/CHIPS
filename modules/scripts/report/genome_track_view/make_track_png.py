#!/usr/bin/env python3
"""GALI BAI
Script to generate a list of genome track view plot for user defined gene list.
Prints out a png files.
"""
import os
import sys
from collections import defaultdict
from optparse import OptionParser
import pandas as pd
import subprocess

def main():
    usage = "USAGE: %prog -i [tracks_all_vlines.ini] -e [extended.bed] -g [list of genes] -o [list of genome track view plots]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-i", "--input", help="vlines addded track files")
    optparser.add_option("-e", "--extend", help="coordinates extended refGene.bed file")
    optparser.add_option("-g", "--genes", action="append", help="list of genes to plot in genome track view")
    optparser.add_option("-o", "--output", action="append", help="list of genome track view plots")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.input:
        optparser.print_help()
        sys.exit(-1)

    lookup_coords = pd.read_csv(options.extend, sep = '\t', header=None, index_col=3).iloc[:,-4]
    #print(options.genes)
    #print(options.output)
    for list_num, gene in enumerate(options.genes):
        #print(list_num)
        #print(gene)
        if gene in pd.read_csv(options.extend, sep = '\t',header=None, index_col=None).iloc[:,3].values:
            region_plot = lookup_coords[gene]
            print(gene)
            print(region_plot)
            os.system("pyGenomeTracks --tracks {input} --region {region} --trackLabelFraction 0.2 --width 38 --dpi 130 -o {output}".format(input = options.input, region = region_plot, output = options.output[list_num]))
            #print(tmp)
            #cmd.append(tmp)
        else:
            print(gene + ' not found')
    #call= " && ".join(cmd)
    #print(call)
    #os.system(call)

if __name__ == '__main__':
    main()
