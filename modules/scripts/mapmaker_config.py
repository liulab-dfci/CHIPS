#!/usr/bin/env python
"""Script to configure a mapmaker run based on the chips run
"""

import os
import sys
import ruamel.yaml
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FILE_1] -f [FILE_2] ...-f [FILE_N] -o [output]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-n", "--names", action="append", help="run names")
    optparser.add_option("-b", "--bed_files", action="append", help="list of _sorted_peaks.narrowPeak.bed files")
    optparser.add_option("-w", "--bw_files", action="append", help="list of _treat_pileup.bw files")
    optparser.add_option("-a", "--bam_files", action="append", help="list of _unique.sorted.dedup.bam files")
    optparser.add_option("-i", "--igv_files", action="append", help="list of _treatment.igv.xml files")
    optparser.add_option("-c", "--config", help="mapmaker config file (to use as a template)")
    optparser.add_option("-o", "--output", help="output file")

    (options, args) = optparser.parse_args(sys.argv)

    #if not options.files or not options.output:
    if not options.names or not options.bed_files or not options.bw_files or not options.bam_files or not options.config or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #READ in the mapmaker config.yaml
    config_f = open(options.config)
    config = ruamel.yaml.round_trip_load(config_f.read())
    config_f.close()

    #BED files
    bed = dict(zip(options.names, options.bed_files))
    config['bed'] = bed

    #BW files
    bw = dict(zip(options.names, options.bw_files))
    config['bigwig'] = bw

    #BAM files
    bam = dict(zip(options.names, options.bam_files))
    config['samples'] = bam

    if options.igv_files:
        igv = dict(zip(options.names, options.igv_files))
        config['cnv'] = igv

    #OUTPUT
    out = open(options.output,"w")
    ruamel.yaml.round_trip_dump(config, out)


if __name__=='__main__':
    main()


