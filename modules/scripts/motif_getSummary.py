#!/usr/bin/env python
"""This script parses a set of motif_list.json, takes the TOP hit of each 
run, and generates a motif summary file (csv) which has the following cols:
Run, MotifID, Name, Logo, Z-Score -
Where Motif = motif's name, logo = path of motif logo img
"""
import os
import sys
from optparse import OptionParser

import json

def parseMDseqPos(motifList):
    """Given a motif_list.json, 
    RETURNS the json represented in the file
    NOTE: the version of MDSeqPos (after commit 492b5a3 in the bitbucket repo)
    outputs motif_list.json

    Actually return just the top z-score
    """
    #KEY in on this string 'var motifList_json = [{...}]
    f = open(motifList)
    #rebuild the motif list
    ret = [json.loads(l.strip()) for l in f]
    f.close()

    #check for no motifs
    #NOTE: I know this is weird but for no motifs found, motif_list.json = {}
    #so we have to check that the first elm is {}
    if ret and ret[0]:
        #extract just index and zscore so we can return the top hit
        tmp = sorted([(i,motif["seqpos_results"]["zscore"]) for (i,motif) in enumerate(ret)], key=lambda x: x[1])
        topHit = ret[tmp[0][0]] #use the index found
        #print(topHit)
        return topHit
    else:
        #NO MOTIFS found
        return {}

def main():
    usage = "USAGE: %prog -m [motif_list.json file 1] -m [motif_list.json file 2] ...  -m [motif_list.json file N] -o [output csv file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-m", "--files", action="append", help="list of files")
    optparser.add_option("-o", "--output", help="output csv file")
    (options, args) = optparser.parse_args(sys.argv)

    if (not options.files or not options.output):
        optparser.print_help()
        sys.exit(-1)

    #GET run names
    #NOTE: the paths to the motif_list looks like this
    #"analysis/motif/SUM44PE_ER/results/motif_list.json" 
    #WE want the 3rd to last!
    runs = [m.split("/")[-3] for m in options.files]
    #the path to the mdseqpos results dir--we simply drop motif_list.json
    #used in getting the right LOGO path!
    results_dir = ["/".join(m.split("/")[:-1]) for m in options.files]

    #PARSE data from the runs
    data = { runName: parseMDseqPos(mdseqpos_out) \
                 for (runName, mdseqpos_out) in zip(runs, options.files)}
    
    #WRITE output 
    f = open(options.output, "w")
    f.write("Runs,MotifID,MotifName,Logo,Zscore\n")
    for (i, r) in enumerate(runs):
        
        if data[r]:
            iid = data[r]['id']
            name = data[r]['factors'][0] if data[r]['factors'] else 'NA'
            #HACK the logo path--found in the results/seqLogo/{id}.png
            logo = os.path.join(results_dir[i], "seqLogo", "%s.png" % iid)
            zscore = "%.2f" % float(data[r]['seqpos_results']['zscore'])
        else:
            #NO motifs found
            iid = name = logo = zscore = 'NA'

        tmp = ",".join([r,iid,name,logo,zscore])
        f.write("%s\n" % tmp)
    f.close()

if __name__=='__main__':
    main()
