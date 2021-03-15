#!/usr/bin/env python

#porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/fastqc/qc.py
from optparse import OptionParser
import sys
import os
import json
import re

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

def _fastqc_parse(input, output=None, param=None):
    data = open(input).readlines()
    sequence_length = 0
    quality_dict = {}
    in_seq_quality_section = False
    for line in data:
        if re.search(r"^Sequence length", line):
            assert sequence_length == 0
            sequence_length = int(re.findall(r"^Sequence length\t(\d+)", line)[0])
        elif re.search(r"^>>Per sequence quality", line):
            assert not in_seq_quality_section
            in_seq_quality_section = True
            continue
        if re.search(r"^>>END_MODULE", line) and in_seq_quality_section:
            in_seq_quality_section = False

        if (not line.startswith("#")) and in_seq_quality_section:
            sequence_quality = re.findall("^(\w+)\t(\w+)", line)[0]
            quality_dict[sequence_quality[0]] = float(sequence_quality[1])
    total = sum(quality_dict.values())
    n = 0
    for item in sorted(quality_dict.items(), key=lambda e: e[0], reverse=True):
        n = n + item[1]
        if n / total > 0.5:
            median = int(item[0])
            break
    return {"sequence_length": sequence_length,
            "median": median}

def json_fastqc(options):
    input_files = []
    for i in options.input:
        input_files.append(str(os.path.abspath(i)))
    sample_list = []
    for i in options.samples:
        sample_list.append(i)

    input={"fastqc_summaries": input_files}
    output={"json": str(os.path.abspath(options.output))}
    param={"ids": sample_list,
           "id": options.run}
    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    # print(input["fastqc_summaries"])
    stat = json_dict["stat"]
    for a_summary, a_id in zip(input["fastqc_summaries"], param["ids"]):
        parsed = _fastqc_parse(input=a_summary)
        stat[a_id] = {}
        stat[a_id]["median"] = parsed["median"]
        stat[a_id]["sequence_length"] = parsed["sequence_length"]
    json_dump(json_dict)


def main():
    USAGE="-i [input_txt] -o [output_json] -s [samples] -r [run]"
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", action = 'append', help="input files")
    optparser.add_option("-o", "--output", help="output json files")
    optparser.add_option("-s", "--samples", action = 'append', help="samples")
    optparser.add_option("-r", "--run", help="run")
    (options, args) = optparser.parse_args(sys.argv)
    json_fastqc(options)


if __name__ == '__main__':
    main()
