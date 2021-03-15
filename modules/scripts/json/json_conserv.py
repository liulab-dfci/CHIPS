#!/usr/bin/env python


#porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/conservation/qc.py

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

def json_conservation(options):
    """
    collect conservation_plot output Phastcon score
    """
    atype = options.TF
    if options.factor:
        atype = options.factor
    if options.basics:
        atype = options.basics
    input={"score": str(os.path.abspath(options.input))}
    output={"json": str(os.path.abspath(options.output))}
    param={"atype": atype, "id": options.run}
    json_dict = {"stat": [], "input": input, "output": output, "param": ""}
    rd = lambda x: str(round(float(x), 3))
    json_dict['stat'] = list(map(rd, open(input['score']).read().strip().split()))
    json_dump(json_dict)


def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", help="input json files")
    optparser.add_option("-o", "--output", help="output json files")
    optparser.add_option("-b", "--basics", help="paramaters: basic")
    optparser.add_option("-f", "--factor", help="paramaters: factor")
    optparser.add_option("-T", "--TF", help="paramaters: TF")
    optparser.add_option("-r", "--run", help="run")
    (options, args) = optparser.parse_args(sys.argv)
    json_conservation(options)


if __name__ == '__main__':
    main()
