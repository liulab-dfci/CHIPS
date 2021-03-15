#!/usr/bin/env python
#porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/replicates/qc.py
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

def json_rep(options):
    """
    input:wigCorrelate of multiple replicates results
          replicates peaks overlap number(percentage: 0.3)
    output: *replicates.json
    """
    overlap_list = [str(os.path.abspath(i)) for i in options.overlap]
    input={"cor": str(os.path.abspath(options.cor)), "overlap": overlap_list}
    output={"json": str(os.path.abspath(options.output))}
    param={"parma": options.run}

    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    json_dict['stat']['cor'] = [ float(i.strip().split()[2]) for i in open(input['cor']).readlines() ]
    json_dict["stat"]['overlap'] = [ float(open(i).read().strip()) for i in input['overlap'] ]
    json_dump(json_dict)


def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-c", "--cor", help="input files")
    optparser.add_option("-o", "--output", help="output files")
    optparser.add_option("-O", "--overlap", aciton='append', help="input overlap")
    optparser.add_option("-r", "--run", help="paramaters: run")
    (options, args) = optparser.parse_args(sys.argv)
    json_rep(options)


if __name__ == '__main__':
    main()
