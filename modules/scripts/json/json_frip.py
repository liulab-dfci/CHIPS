#!/usr/bin/env python
# porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/frip/qc.py

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

def json_frip(options):
    input_list = []
    for i in options.input:
        input_list.append(str(os.path.abspath(i)))
    input={"frip": input_list}
    output={"json": str(os.path.abspath(options.output))}
    param={"samples":options.samples}
    """
    input is *_frip.txt
    like:
    #reads_under_peaks  501359
    #total_reads    3998454
    output is conf.json_prefix + "_frip.json"
    param for matching samples
    """
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    for i,s in zip(input["frip"], param["samples"]):
        with open(i) as inf:
            data = []
            for item in inf.readlines():
                item = item.rstrip("\n").split("\t")
                data.append(item)
        # print(data)
        json_dict["stat"][s] = {}
        json_dict["stat"][s]["info_tag"] = int(data[0][1])
        json_dict["stat"][s]["total_tag"] = int(data[1][1])
        json_dict["stat"][s]["frip"] = float(int(data[0][1]))/int(data[1][1])
    json_dump(json_dict)


def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", action='append', help="input frip files")
    optparser.add_option("-o", "--output", help="output json files")
    optparser.add_option("-s", "--samples", action='append', help="paramaters: samples")
    (options, args) = optparser.parse_args(sys.argv)
    json_frip(options)


if __name__ == '__main__':
    main()
