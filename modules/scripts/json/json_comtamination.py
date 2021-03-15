#!/usr/bin/env python

#porting https://github.com/cfce/chilin/blob/master/chilin2/modules/contamination/qc.py

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

def json_contamination(options):
    # input file :
    """
    hg19 97.8341
    mm9 3.90596
    dm3 1.305
    S_cerevisiae 0.134
    ...
    """
    input={"summaries": str(os.path.abspath(options.input))}
    output={"json": str(os.path.abspath(options.output))}
    param={"samples": options.samples ,"id": options.ID}
    
    with open(os.path.abspath(options.input)) as f:
        data = []
        for line in f.readlines():
            row = line.rstrip("\n").split(" ")
            data.append(row)
        # print(data)
    species = []
    for i in data:
        species.append(i[0])
    library_contamination = {}
    library_contamination["meta"] = {"sample": param["id"], "species": species}
    library_contamination["value"] = {}
    # print(len(species))
    for i in range(len(species)):
        library_contamination["value"][species[i]] = data[i][1]
    # for a_summary, s in zip(input["summaries"], list(map(underline_to_space, param["samples"]))):
    #     ## each bowtie_summary has several species information
    #     library_contamination["value"][s] = {}
    #     for i, j in zip(a_summary, param["species"]):
    #         ## species 1, species2, species3
    #         mapped = int(open(i[0]).readlines()[2].strip().split()[0])
    #         total = int(open(i[1]).read().strip())
    #         library_contamination["value"][s][j] = float(mapped)/total

    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    json_dict["stat"] = library_contamination
    json_dump(json_dict)


def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", help="input files")
    optparser.add_option("-o", "--output", help="output files")
    optparser.add_option("-s", "--samples", help="paramaters: samples")
    optparser.add_option("-I", "--ID", help="ID")
    (options, args) = optparser.parse_args(sys.argv)
    json_contamination(options)


if __name__ == '__main__':
    main()
