#!/usr/bin/env python
#porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/ceas/qc.py
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


def json_meta(options):
    """
    generate json of genomic distribution (given by the bedAnnotate output)
    ***THE key difference between json_meta and this fn is that bedAnnotate
    conveniently outputs the distribution as a dictionary of peak counts
    """
    input={"meta": str(os.path.abspath(options.input))}
    # file "_summary.txt"
    output={"json": str(os.path.abspath(options.output))}
    param={"id": options.ID}
    f = open(input["meta"])
    #f = something like: {'Intron': 68017, 'Exon': 7659, 'Intergenic': 73090, 'Promoter': 11229}
    content = eval(f.read())
    total = 0 
    for k in content.keys():
        total += content[k]
    json_dict = {"input": input, "stat": {}, "output": output, "param": param}
    json_dict["stat"]["exon"] = content['Exon']/float(total) if float(total) else 0.0
    json_dict["stat"]["intron"] = content['Intron']/float(total) if float(total) else 0.0
    json_dict["stat"]["promoter"] = content['Promoter']/float(total) if float(total) else 0.0
    json_dict["stat"]["inter"] = content['Intergenic']/float(total) if float(total) else 0.0
    f.close()
    json_dump(json_dict)


def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", help="input files")
    optparser.add_option("-o", "--output", help="output files")
    optparser.add_option("-I", "--ID", help="ID")
    (options, args) = optparser.parse_args(sys.argv)
    json_meta(options)


if __name__ == '__main__':
    main()
