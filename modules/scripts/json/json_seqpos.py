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

def json_seqpos(options):

    input = {"seqpos": str(os.path.abspath(options.input))}
    output = {"json": str(os.path.abspath(options.output))}
    param={"prefix": str(os.path.abspath(options.prefix)), "z_score_cutoff": -1}

    z_score_cutoff = param["z_score_cutoff"]
    seqpos_html_content = open(input['seqpos']).readlines()
    mdseqpos_result = []
    if seqpos_html_content == ['{}\n']:
        top_motifs = 'None'
    else:
        ## parse motif list json file
        for m in seqpos_html_content:
            mdseqpos_result.append(json.loads(m.strip()))
        satisfied_motif_list = []

        for a_motif in mdseqpos_result:
            if a_motif['seqpos_results']['zscore'] == 'None':
                a_motif['seqpos_results']['zscore'] = 65535
            if a_motif['factors'] == None:
                a_motif['factors'] = ['denovo']
            satisfied_motif_list.append(a_motif)

        satisfied_motif_list.sort(key=lambda x:x['seqpos_results']['zscore'])
        satisfied_count = 0
        top_motifs = []
        for a_motif in satisfied_motif_list:

            if a_motif['id'].find('observed')>0:
                continue
            if satisfied_count == 10:
                break

            # z_score is a negative score, the smaller, the better
            if a_motif['seqpos_results']['zscore'] < z_score_cutoff :
                satisfied_count += 1
                top_motifs.append(a_motif)

        ## choose first 5 motifs to fit into latex document

        for n, _ in enumerate(top_motifs):
            top_motifs[n]["logoImg"] = param["prefix"] + top_motifs[n]['id'] + ".png"

    result_dict = {"stat": {}, "input": input, "output": output, "param": param}
    result_dict["stat"]["satisfied_motifs"] = top_motifs
    json_dump(result_dict)


def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-i", "--input", help="input json files")
    optparser.add_option("-o", "--output", help="output json files")
    optparser.add_option("-p", "--prefix", help="seqlogo directory")
    (options, args) = optparser.parse_args(sys.argv)
    json_seqpos(options)


if __name__ == '__main__':
    main()
