#!/usr/bin/env python

# porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/macs/qc.py
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

def _peaks_parse(input):
    total = 0
    fc20n = 0
    fc10n = 0
    peaks_info = {}
    with open(input) as peaks_xls:
        for line in peaks_xls:
            if line.startswith('# tags after filtering in treatment'):
                # tags after filtering in treatment: 13438948
                peaks_info["treat_unic"] = int(line.strip().split(":")[1])
            if line.startswith("# Redundant rate in treatment: "):
                peaks_info["treat_unic_ratio"] = 1 - float(line.split(':')[1])

            if line.startswith('# tags after filtering in control'):
                peaks_info["control_unic"] = int(line.strip().split(":")[1])
            if line.startswith("# Redundant rate in control: "):
                peaks_info["control_unic_ratio"] = 1 - float(line.split(':')[1])

            if line.startswith('# d'):
                peaks_info["distance"] = int(line.strip().split("=")[1])
            if line.strip() != "" and not line.startswith("#") and not line.startswith("chr\t"):
                l = line.strip().split("\t")
                total += 1
                ## column 7th denotes fold change value
                fc = float(l[7])
                if fc >= 20:
                    fc20n += 1
                if fc >= 10:
                    fc10n += 1
            if line.startswith("# qvalue cutoff"):
                q_value_cutoff = float(line.split('=')[1])
            if line.startswith("# d"): # parse shift-size, # d =
                shift_size = int(line.strip().split("=")[1])/2

    peaks_info["totalpeak"] = total
    peaks_info["peaksge20"] = fc20n
    peaks_info["peaksge10"] = fc10n
    if peaks_info["totalpeak"] >= 200:
        peaks_info["peaksge20ratio"] = peaks_info["peaksge20"] / peaks_info["totalpeak"]
        peaks_info["peaksge10ratio"] = peaks_info["peaksge10"] / peaks_info["totalpeak"]
    elif 0 < peaks_info["totalpeak"] < 200:
        peaks_info["peaksge20ratio"] = peaks_info["peaksge20"] / peaks_info["totalpeak"]
        peaks_info["peaksge10ratio"] = peaks_info["peaksge10"] / peaks_info["totalpeak"]
        print >> sys.stderr, "Warning: peaks so few for motif scan"
    elif peaks_info["totalpeak"] == 0 :
        peaks_info["peaksge20ratio"] = 0
        peaks_info["peaksge10ratio"] = 0
        print >> sys.stderr, "Warning: 0 peaks"

    peaks_info["qvalue"] = q_value_cutoff
    peaks_info["shiftsize"] = shift_size
    return peaks_info

def json_macs2_rep(options):
    """
    collect replicates macs2 info to json files
    compared to merged one, collect redundant ratio with --keep-dup 1 option
    """
    input_list = []
    for i in options.input:
        input_list.append(str(os.path.abspath(i)))
    sample_list = []
    for i in options.samples:
        sample_list.append(i)
    input={"all_peak_xls": input_list}
    output={"json": str(os.path.abspath(options.output))}
    param={"samples": sample_list}

    json_dict = {"stat": {}, "input": input, "output": output, "param": param}
    parsed = []
    for i in input["all_peak_xls"]:
        if os.path.exists(i): ## in case only broad peaks would break down sometimes, narrowPeak very seldom no peaks
            parsed.append(_peaks_parse(i))

    if all(list(map(os.path.exists, input['all_peak_xls']))):
        for sample, stat in zip(param["samples"], parsed):
            json_dict["stat"][sample] = stat
        json_dump(json_dict)


def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", action="append", help="input xls files")
    optparser.add_option("-o", "--output", help="output json files")
    optparser.add_option("-s", "--samples", action="append", help="paramaters: samples")
    (options, args) = optparser.parse_args(sys.argv)
    json_macs2_rep(options)


if __name__ == '__main__':
    main()
