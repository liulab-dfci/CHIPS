#!/usr/bin/env python3
"""GALI BAI
Script to summary conservation and motif finding results
When homer is on, prints out a csv file with header- 'Sample', 'Conservation', 'Motif', 'Homer Motif Logo', 'Negative Log P-value'
Otherwise, prints out a csv file with header- 'Sample', 'Conservation'
"""

import os
import sys
from collections import defaultdict
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -c [input {runRep}_conserv_thumb.png] -m [input known1.logo.png] -t [input KnownResults.txt] -p [input path to report/Downstream/ ] -o [output conservation_and_top_motifs.csv]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-c", "--conserv", action="append", help="conservation plots from analysis/conserv")
    optparser.add_option("-m", "--homer", action="append", help ="homer motif plots from analysis/motif")
    optparser.add_option("-t", "--motif", action="append", help="motif txt summary file rom analysis/motif")
    optparser.add_option("-p", "--outpath", help="downstream png file output path")
    optparser.add_option("-o", "--output", help="output .csv file")
    (options, args) = optparser.parse_args(sys.argv)

    if  not options.output:
        optparser.print_help()
        sys.exit(-1)

    if not os.path.exists(options.outpath + "/plots_conserv"):
        os.mkdir(options.outpath + "/plots_conserv")
    if not os.path.exists(options.outpath + "/plots_motif"):
        os.mkdir(options.outpath + "/plots_motif")

    conserv_out=defaultdict(list)
    for file in options.conserv:
        os.system("cp {i} {path}/plots_conserv/".format(i = file, path = options.outpath))
        fname = os.path.basename(file)
        sname = fname.split("_")[0]
        path = os.path.join("img:Downstream/plots_conserv/" + fname)
        conserv_out[sname].append(path)

    if not options.homer:
        print("homer is disabled")
        out = open(options.output, "w")
        #write header
        out.write("%s\n" % ",".join(['Sample', 'Conservation']))
        for k in conserv_out.keys():
            s = [str(k)]
            l = conserv_out[k]
            s.extend(l)
            print(s)
            c = ",".join(s)
            out.write("%s\n" % c)

        out.close()

    else:
        print("homer is enabled")
        homer_out=defaultdict(list)
        for file in options.homer:
            fname = file.split("/")[-4:]
            sname = file.split("/")[-4]
            pname = "_".join(fname)
            #print(pname)
            os.system("cp {i} {path}/plots_motif/{pname}".format(i = file, path = options.outpath, pname = pname))
            path = os.path.join("img:Downstream/plots_motif/" + pname)
            homer_out[sname].append(path)

        motif_out=defaultdict(list)
        pvalue_out=defaultdict(list)
        for file in options.motif:
            sname = file.split("/")[-3]
            fhd = open(file, "rt")
            line = fhd.readlines()
            line = line[1].strip().split("\t")
            m = line[0].split('(')[0]
            #print(m)
            p = str(float(line[3])*-1)
            #print(p)
            motif_out[sname].append(m)
            pvalue_out[sname].append(p)
            fhd.close()

        out = open(options.output, "w")
        #write header
        out.write("%s\n" % ",".join(['Sample', 'Conservation', 'Motif', 'Homer Motif Logo', 'Negative Log P-value']))
        for k in conserv_out.keys():
            s = [str(k)]
            l = conserv_out[k] + motif_out[k] + homer_out[k] + pvalue_out[k]
            #print(l)
            s.extend(l)
            #print(s)
            c = ",".join(s)
            out.write("%s\n" % c)

        out.close()
if __name__ == '__main__':
    main()
