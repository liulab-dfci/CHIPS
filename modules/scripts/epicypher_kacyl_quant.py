#!/usr/bin/env python
"""Script to qunatify the reads on chr_epicypher
"""

import os
import sys
from optparse import OptionParser
import subprocess

_SPIKE=['WT-A',
        'WT-B',
        'H3K4ac-A',
        'H3K4ac-B',
        'H3K9ac-A',
        'H3K9ac-B',
        'H3K9bu-A',
        'H3K9bu-B',
        'H3K14ac-A',
        'H3K14ac-B',
        'H3K18ac-A',
        'H3K18ac-B',
        'H3.2tetraAc-A',
        'H3.2tetraAc-B',
        'H3K27ac-A',
        'H3K27ac-B',
        'H4K8ac-A',
        'H4K8ac-B',
        'H4K16ac-A',
        'H4K16ac-B',
        'H4tetraAc-A',
        'H4tetraAc-B']

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
    #_STOP = 32*_STEP + _START
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
    #printout
    
    print(total)
        
if __name__=='__main__':
    main()
