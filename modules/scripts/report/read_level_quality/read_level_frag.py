#!/usr/bin/env python3
"""GALI BAI
Script to draw histogram on each fragment file and concatenate them.
Prints out a csv file with fragment length as rawname.
"""

import os
import sys
import pandas as pd
from collections import defaultdict
from optparse import OptionParser
import numpy as np

def main():
    usage = "USAGE: %prog -c [input contamination table] -o [output csv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--frag", action="append", help="a list of _frag.txt files in analysis/frag/")
    optparser.add_option("-o", "--output", help="output .plotly file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.frag:
        optparser.print_help()
        sys.exit(-1)

    frag_out=defaultdict(list)

    for file in options.frag:
        frag = pd.read_csv(file, header=None)
        fpath = os.path.basename(file)
        fname = os.path.splitext(fpath)[0]
        sname = fname.split("_frags")[0]
        counts, bins = np.histogram(frag, bins=np.linspace(1,1000,200), density=True)
        frag_df = pd.DataFrame(counts).reset_index()
        frag_df.columns = ['Bin','Density']
        frag_df['Fragment size in bp'] = frag_df['Bin']*5
        frag_df = frag_df.set_index('Fragment size in bp')
        frag_out[sname] =frag_df['Density']

    dfs = pd.concat(frag_out, axis=1)
    dfs.to_csv(options.output)

if __name__ == '__main__':
    main()
