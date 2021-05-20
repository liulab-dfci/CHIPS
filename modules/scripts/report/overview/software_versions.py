#!/usr/bin/env python3
"""Len Taing (TGBTG)
Script to check softeare versions for fastp|fastqc|bwa|macs2|bedtools|home
Prints out a tsv with header- Name\tVersion\tBuild\tChannel
"""

import os
import sys
import subprocess
from optparse import OptionParser

    
def main():
    usage = "USAGE: %prog -o [output tsv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-o", "--output", help="output .tsv file")
    (options, args) = optparser.parse_args(sys.argv)
    
    if  not options.output:
        optparser.print_help()
        sys.exit(-1)

    #use conda to get the software versions
    software_ls = ['fastp','fastqc','bwa','macs2','bedtools','homer']
    tmp = """conda list"""
    cmd = tmp.split(" ")
    #print(cmd)
    output, err = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    
    #SPLIT the output by lines
    conda_output = output.decode("utf-8").strip().split("\n")

    #write output
    out = open(options.output, "w")
    #Write header-
    out.write("%s\n" % "\t".join(["Name", "Version", "Build", "Channel"]))

    #Scan for the softwares we are interested in
    for l in conda_output:
        tmp = l.split()
        if tmp and tmp[0] in software_ls:
            out.write("%s\n" % "\t".join(tmp))

    out.close()
 
if __name__ == '__main__':
    main()
