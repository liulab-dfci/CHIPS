#!/usr/bin/env python
"""Script to qunatify the reads on chr_epicypher
"""

import os
import sys
from optparse import OptionParser
import subprocess

_SPIKE=['WT-A',
        'WT-B',
        'H3K4me1-A',
        'H3K4me1-B',
        'H3K4me2-A',
        'H3K4me2-B',
        'H3K4me3-A',
        'H3K4me3-B',
        'H3K9me1-A',
        'H3K9me1-B',
        'H3K9me2-A',
        'H3K9me2-B',
        'H3K9me3-A',
        'H3K9me3-B',
        'H3K27me1-A',
        'H3K27me1-B',
        'H3K27me2-A',
        'H3K27me2-B',
        'H3K27me3-A',
        'H3K27me3-B',
        'H3K36me1-A',
        'H3K36me1-B',
        'H3K36me2-A',
        'H3K36me2-B',
        'H3K36me3-A',
        'H3K36me3-B',
        'H4K20me1-A',
        'H4K20me1-B',
        'H4K20me2-A',
        'H4K20me2-B',
        'H4K20me3-A',
        'H4K20me3-B']


def main():
    usage = "USAGE: %prog -b <epicypher bam file>"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-b", "--bam", help="epicypher sorted and indexed aligned reads file--bam")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.bam:
        print(usage)
        sys.exit()

    bam = options.bam
    _START = 100
    _STEP = 1000
    _SIZE = 148 #length of the spike-in sequences
    _STOP = len(_SPIKE)*_STEP+_START

    #cts = []
    total = 0
    for (i, s) in enumerate(range(_START, _STOP, _STEP)):
        ##CALL samtools view -c [bam] REGION for each spikein
        call = "samtools view -c %s chr_epicypher:%s-%s" %(bam,s, s+_SIZE)
        #print(call)
        tmp = subprocess.run(call.split(" "), stdout=subprocess.PIPE)
        #cts.append((_SPIKE[i], int(tmp.stdout)))

        ct = int(tmp.stdout)
        print(_SPIKE[i], ct)
        total += ct
        
    #printout total reads
    print(total)
        
if __name__=='__main__':
    main()
