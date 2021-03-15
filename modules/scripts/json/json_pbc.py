#!/usr/bin/env python

# porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/library/qc.py
from optparse import OptionParser
import sys
import os
import json

def json_dump(json_dict):   # json
    """
    dump out uniform json files for collecting statistics
    :param json_dict: output python dict to json
    :return: json_file name
    """
    json_file = json_dict["output"]["json"]
    with open(json_file, "w") as f:
        json.dump(json_dict, f, indent=4)
    return json_file

def json_pbc(options):
    """
    collect conservation_plot output Phastcon score
    """
    input_list = []
    for i in options.input:
        input_list.append(str(os.path.abspath(i)))
    input={"pbc": input_list}
    output={"json": str(os.path.abspath(options.output))}

    json_dict = {"stat": {}, "input": input, "output": output}
    
    for i in input_list:
        sampleID = i.strip().split("/")[-1].split('.')[0]
        #remove the _pbc ending
        sampleID = sampleID.replace("_pbc", "")
        f = open(i)
        firstLine = f.readline().strip().split()
        #print(firstLine)
        N1 = int(firstLine[1]) #record number of reads with just one location
        #NOW SUM over the rest of the reads
        Nd = N1
        for l in f:
            tmp = l.strip().split()
            Nd += int(tmp[1])

        json_dict["stat"][sampleID] = {}
        json_dict["stat"][sampleID]["N1"] = N1
        json_dict["stat"][sampleID]["Nd"] = Nd
        json_dict["stat"][sampleID]["PBC"] = round(N1/Nd, 3)
        json_dict["param"]={"samples": []}
        json_dict["param"]['samples'].append(sampleID)
        f.close()
    json_dump(json_dict)

def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", action='append', help="input files")
    optparser.add_option("-o", "--output", help="output files")
    (options, args) = optparser.parse_args(sys.argv)
    json_pbc(options)


if __name__ == '__main__':
    main()
