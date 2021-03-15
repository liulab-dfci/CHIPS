#!/usr/bin/env python

#porting from https://github.com/cfce/chilin/blob/master/chilin2/modules/macs2_fragment/qc.py

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

def get_size(rscript):
    values = {}
    with open(rscript) as model:
        for line in model:
            if line.startswith("p <- "):
                values["positive"] = line
            if line.startswith("m <- "):
                values["minus"] = line
            if line.startswith("x <- "):
                values["x"] = line
            if line.startswith("xcorr"):
                values['xcorr'] = line
            if line.startswith("ycorr"):
                values['ycorr'] = line
    return values

def json_frag(options):
    """ parse macs2 predictd r file into json file
    """
    input_list = [str(os.path.abspath(i)) for i in options.Rscript]
    output_list = [str(os.path.abspath(i)) for i in options.outputR]
    input={"r": input_list}
    output = {"json": str(os.path.abspath(options.output)), "r": output_list}
    param = {"samples": options.samples, "frag_tool": options.fragTool}

    json_dict = {"input": input, "output": output, "param": param, "stat": {}}
    # for rin, rout, s in zip(input["r"], output["r"], param["samples"]):
    # rin = input["r"]
    # rout = output["r"]
    # s = param["samples"]
    # values = get_size(rin)
    # with open(rout, 'w') as f:
    #     f.write(values['positive'])
    #     f.write(values['minus'])
    #     f.write(values['xcorr'])
    #     f.write(values['ycorr'])
    #     f.write("xcorr.max = xcorr[which(ycorr==max(ycorr))]\n")
    #     f.write(values['x'])
    #     f.write("p.expect = sum(x * p/100) \n")
    #     f.write("m.expect = sum(x * m/100) \n")
    #     f.write("p.sd = sqrt(sum(((x-p.expect)^2)*p/100)) \n")
    #     f.write("m.sd = sqrt(sum(((x-m.expect)^2)*m/100)) \n")
    #     f.write("cat(paste((p.sd + m.sd)/2, '\t', xcorr.max)) \n")
    # f.close()
    # std_frag = os.popen("Rscript %s" % rout).read().strip().split()
    # json_dict["stat"][s] = "%s" % (int(float(std_frag[1])))
    
    # json_dump(json_dict)
    for rin, rout, s in zip(input["r"], output["r"], param["samples"]):
        values = get_size(rin)
        with open(rout, 'w') as f:
            f.write(values['positive'])
            f.write(values['minus'])
            f.write(values['xcorr'])
            f.write(values['ycorr'])
            f.write("xcorr.max = xcorr[which(ycorr==max(ycorr))]\n")
            f.write(values['x'])
            f.write("p.expect = sum(x * p/100) \n")
            f.write("m.expect = sum(x * m/100) \n")
            f.write("p.sd = sqrt(sum(((x-p.expect)^2)*p/100)) \n")
            f.write("m.sd = sqrt(sum(((x-m.expect)^2)*m/100)) \n")
            f.write("cat(paste((p.sd + m.sd)/2, '\t', xcorr.max)) \n")
        f.close()
        std_frag = os.popen("Rscript %s" % rout).read().strip().split()
        json_dict["stat"][s] = "%s" % (int(float(std_frag[1])))
    json_dump(json_dict)


def main():
    USAGE=""
    optparser = OptionParser(usage=USAGE)
    optparser.add_option("-o", "--output", help="output json files")
    optparser.add_option("-r", "--Rscript", action='append', help="input Rscript")
    optparser.add_option("-R", "--outputR", action='append', help="output Rscript")
    optparser.add_option("-s", "--samples", action='append', help="samples")
    optparser.add_option("-f", "--fragTool", help="frag tool")
    (options, args) = optparser.parse_args(sys.argv)
    json_frag(options)


if __name__ == '__main__':
    main()
