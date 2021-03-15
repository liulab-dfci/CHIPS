#!/usr/bin/env python

#porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/enrichment/dc.py
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

def json_enrich_meta(options):
    """ enrichment in meta regions
    """
    meta_list = []
    for i in options.input:
        meta_list.append(str(os.path.abspath(i)))
    mapped_list = []
    for i in options.mapped:
        mapped_list.append(str(os.path.abspath(i)))
    dhs_list = []
    for i in options.DHS:
        dhs_list.append(str(os.path.abspath(i)))
    input = {'meta': meta_list, 'mapped': mapped_list}
    # meta "ceas summary txt"; mapped "align mapping"
    output = {"json": str(os.path.abspath(options.output))}
    subsample = False
    if options.down:
        subsample = True
    param = {'dhs': dhs_list, 'down': subsample, 
             'has_dhs':str(os.path.abspath(options.HasDHS)), 'id':options.run, 'samples':options.samples}
    # dhs "ceas dhs"; has_dhs config["DHS"]

    json_dict = {"stat": {}, "input": input, "output": output, "param":param}
    
    for n, s in enumerate(param['samples']):
        ## total mapped reads
        with open(input["mapped"][n]) as input_mapped:
            mapped = float(input_mapped.readlines()[4].split()[0])
        # print(mapped)
        json_dict['stat'][s] = {}
        with open(input['meta'][n]) as input_meta:
            meta_string = input_meta.readline().replace('\'','\"')
            meta_dict = json.loads(meta_string)
        # meta is something like: {'Exon': 68017, 'Promoter': 7659, 'DHS': 73090, 'Total': 11229}
        meta = [meta_dict['Exon'], meta_dict['Promoter'], meta_dict['DHS'], meta_dict['Total']]
        # print(meta)
        if not param["down"]:
            json_dict['stat'][s]['exon'] = meta[0]/mapped
            json_dict['stat'][s]['promoter'] = meta[1]/mapped ## use all mapped reads
        else:
            json_dict['stat'][s]['exon'] = meta[0]/meta[3]
            json_dict['stat'][s]['promoter'] = meta[1]/meta[3] ## use 4M reads

        if param['has_dhs']:
            if not param["down"]:
                json_dict['stat'][s]['dhs'] = meta[2]/mapped
            else:
                json_dict['stat'][s]['dhs'] = meta[2]/meta[3]
    json_dump(json_dict)


def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", action = 'append', help="input meta files")
    optparser.add_option("-m", "--mapped",  action = 'append', help="input mapped files")
    optparser.add_option("-o", "--output", help="output json files")
    optparser.add_option("-D", "--DHS", action = 'append', help="set DHS")
    optparser.add_option("-H", "--HasDHS", help="set has dhs")
    optparser.add_option("-d", "--down", help="set down")
    optparser.add_option("-r", "--run", help="set run")
    optparser.add_option("-s", "--samples", action = 'append', help="set samples")
    (options, args) = optparser.parse_args(sys.argv)
    json_enrich_meta(options)


if __name__ == '__main__':
    main()
