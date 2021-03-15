#!/usr/bin/env python
"""This script parses a set of homerResults.html, takes the TOP hit of each 
run, and generates a motif summary file (csv) which has the following cols:
Run, MotifID, Name, Logo, Z-Score -
Where Motif = motif's name, logo = path of motif logo img
"""
import os
import sys
from optparse import OptionParser

import json
import xml.etree.ElementTree as ET

def parseHomer(motifList):
    """Given a homerResults.html, 
    RETURNS the json represented in the topHit
    NOTE: the version of MDSeqPos (after commit 492b5a3 in the bitbucket repo)
    """
    #TODO: need 
    if os.stat(motifList).st_size == 0:
        #EMPTY file
        return {}
    else:
        #KEY in on Table's second TR
        tree = ET.parse(motifList)
        root = tree.getroot()
        body = root.findall('BODY')[0]
        table = body.findall('TABLE')[0]
        
        #NOTE: TR[0] is HEADER, TR[1] is top hit
        topHit = table[1]
        logo = topHit[1][0].attrib['src']
        pval = topHit[2].text
        logp = topHit[3].text
        #NOTE: the Best Match/Details is separated by /--just take the 1st part
        #NEED to be sure to replace ','
        name = topHit[7].text.split("/")[0].replace(',','-')
        
        return {'name': name, 'logo': logo, 'pval': pval, 'logp': logp}

def main():
    usage = "USAGE: %prog -m [homerResults.html file 1] -m [homerResults.html file 2] ...  -m [homerResults.html file N] -o [output csv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-m", "--files", action="append", help="list of files")
    optparser.add_option("-o", "--output", help="output csv file")
    (options, args) = optparser.parse_args(sys.argv)

    if (not options.files or not options.output):
        optparser.print_help()
        sys.exit(-1)

    #GET run names
    #NOTE: the paths to the motif_list looks like this
    #"analysis/motif/SUM44PE_ER.rep1/results/motif_list.json" 
    #WE want the 3rd to last!
    runs = [m.split("/")[-3] for m in options.files]
    #the path to the homerResults dir--we simply drop motif_list.json
    #used in getting the right LOGO path!
    results_dir = ["/".join(m.split("/")[:-1]) for m in options.files]

    #PARSE data from the runs
    data = { runName: parseHomer(homer_out) \
                 for (runName, homer_out) in zip(runs, options.files)}
    

    #WRITE output 
    f = open(options.output, "w")
    f.write("Runs,MotifName,Logo,Pval,logPval\n")
    for (i, r) in enumerate(runs):
        if data[r]:
            name = data[r]['name']
            #PATH is 
            logo = os.path.join(results_dir[i], data[r]['logo'])
            pval = data[r]['pval']
            logp = data[r]['logp']
        else:
            #NO motifs found
            name = logo = pval = logp = 'NA'

        tmp = ",".join([r,name,logo,pval,logp])
        f.write("%s\n" % tmp)
    f.close()

if __name__=='__main__':
    main()
