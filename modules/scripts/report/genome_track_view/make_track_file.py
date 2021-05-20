#!/usr/bin/env python3
"""GALI BAI
Script to create configuration file for pyGenomeTracks.
Prints out a tracks_all_vlines.ini
"""

import os
import sys
from collections import defaultdict
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -i [tracks_all.ini] -e [coord extended bed file] -t [tss bed file] -o [tracks_all_vlines.ini]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-i", "--input", help="initial track files for all samples")
    optparser.add_option("-e", "--extend", help="coordinates extended refGene.bed file")
    optparser.add_option("-t", "--tss", help="bed file with only tss coordinates of refGene.bed")
    optparser.add_option("-o", "--output", help="vlines addded track files")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.input:
        optparser.print_help()
        sys.exit(-1)

    with open(options.input) as f:
        with open(options.output, "w") as f1:
            for line in f:
                f1.write(line)
            f1.write("\n")
            f1.write("[spacer]\n")
            f1.write("[gene]\n")
            f1.write("file={ref_file_to_use}\n".format(ref_file_to_use = options.extend))
            f1.write("title=Genes\n")
            f1.write("height=2\n")
            f1.write("fontsize=10\n")
            f1.write("file_type=bed\n")
            f1.write("\n")
            f1.write("[vlines]\n")
            f1.write("file={ref_file_to_use}\n".format(ref_file_to_use = options.tss))
            f1.write("type=vlines\n")
    f.close()
    f1.close()

if __name__ == '__main__':
    main()
