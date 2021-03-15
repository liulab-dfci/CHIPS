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


def json_dhs(options):
    input={"top_peaks": 5000, "dhs_peaks": str(os.path.abspath(options.input))}
    output={"json": str(os.path.abspath(options.output))}
    param={}
    result_dict = {"stat": {}, "input": input, "output": output, "param": param}
#    result_dict["stat"] = float(open(input["dhs"]).read().strip())
    content = open(input["dhs_peaks"]).read().strip().split(",")
    # first number is 5000, second is dhs, third is 5000+dhs
    result_dict["stat"]["overlap"] = int(content[1])
    result_dict["stat"]["number"] = int(content[0])
    json_dump(result_dict)


def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", help="input files")
    optparser.add_option("-o", "--output", help="output files")
    # optparser.add_option("-I", "--ID", help="ID")
    (options, args) = optparser.parse_args(sys.argv)
    json_dhs(options)


if __name__ == '__main__':
    main()
