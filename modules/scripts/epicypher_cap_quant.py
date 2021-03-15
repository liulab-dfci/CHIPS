#!/usr/bin/env python
"""Script to qunatify the reads on chr_epicypher
"""

import os
import sys
from optparse import OptionParser
import subprocess

_SPIKE=['Unmodified-100',
	'Unmodified-80',
	'Unmodified-60',
	'Unmodified-40',
	'Unmodified-20',
	'Unmodified-10',
	'H3K4me1-100',
	'H3K4me1-80',
	'H3K4me1-60',
	'H3K4me1-40',
	'H3K4me1-20',
	'H3K4me1-10',
	'H3K4me2-100',
	'H3K4me2-80',
	'H3K4me2-60',
	'H3K4me2-40',
	'H3K4me2-20',
	'H3K4me2-10',
	'H3K4me3-100',
	'H3K4me3-80',
	'H3K4me3-60',
	'H3K4me3-40',
	'H3K4me3-20',
	'H3K4me3-10']

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
