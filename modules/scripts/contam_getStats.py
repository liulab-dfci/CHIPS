#!/usr/bin/env python
"""Script to collect the contamination stats
Outputs: 
Sample,Species1,Species2,...,SpeciesN
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FILE_1] -f [FILE_2] ...-f [FILE_N] -o [output]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list files")
    optparser.add_option("-o", "--output", help="output file")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #STORE it as a dictionary of samples {SAMPLEID: {'assembly':val}...}
    contam = {}
    samples = []
    for f in options.files:
        #TRY to infer the SAMPLE NAMES
        sampleID = f.strip().split("/")[-1].split('.')[0]
        if sampleID.endswith('_contamination'):
            sampleID = sampleID.replace("_contamination","")
        samples.append(sampleID)

        f = open(f)
        tmp = [l.strip().split(" ") for l in f]
        #convert the percentage to only 1-decimal number
        tmp = [(species, "%.1f" % float(val)) for (species,val) in tmp]
        f.close()
        #CONVERT to dictionary
        contam[sampleID] = { species : val for (species,val) in tmp}

    #OUTPUT
    samples = sorted(samples)
    assemblies = sorted(list(contam[samples[0]].keys()))

    out = open(options.output,"w")
    out.write("%s\n" % ",".join(["Sample"] + assemblies))
    for s in samples:
        vals = [contam[s][a] for a in assemblies]
        out.write("%s\n" % ",".join([s] + vals))

    out.close()

if __name__=='__main__':
    main()


