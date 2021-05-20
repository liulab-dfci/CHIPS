#!/usr/bin/env python3
"""GALI BAI
Script to sum value in myco contamination into one column
Prints out a csv file with header- "Sample", "S_cerevisiae", "dm3", "e_coli", "myco"
"""

import os
import sys
import pandas as pd
from collections import defaultdict
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -c [input contamination table] -o [output csv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-c", "--contam", help="contamination.csv file in analysis/contam/")
    optparser.add_option("-o", "--output", help="output .plotly file")
    (options, args) = optparser.parse_args(sys.argv)

    ##adds up columns start with "myco_"
    df = pd.read_csv(options.contam,index_col=0)
    myco_cols = [col for col in df.columns if 'myco_' in col]
    df['myco'] = df[myco_cols].sum(axis=1)
    cols_to_include = [col2 for col2 in df.columns if 'myco_' not in col2]
    df[cols_to_include].to_csv(options.output)


if __name__ == '__main__':
    main()
