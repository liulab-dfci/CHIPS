#!/usr/bin/env python
"""Script to generate an igv session given:
1. paths to all of the .treat.bw files
2. the genome used, e.g. hg19

NOTE: uses cidc_chips/static/chips_igv.session.xml as a TEMPLATE and fills in 
the resource and track info
ref: https://software.broadinstitute.org/software/igv/Sessions

OUTPUT: VALID IGV session that contains these files (to output file)
e.g. analysis/peaks/all_treatments.igv.session
"""

import os
import sys
from optparse import OptionParser

import xml.etree.ElementTree as ET
import copy

_attribs = {'altColor':"0,0,178", 'autoScale':"false",
            'clazz':"org.broad.igv.track.DataSourceTrack",
            'color':"0,0,178", 'displayMode':"COLLAPSED", 
            'featureVisibilityWindow':"-1", 'fontSize':"10",
            'normalize':"false", 'renderer':"BAR_CHART", 
            'sortable':"true", 'visible':"true",'windowFunction':"mean"}


def main():
    usage = "USAGE: %prog -g [genome e.g. hg19] -t [treatment.bw 1] -t [treatment.bw 2] ... -t [treatment.bw N]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-g", "--genome", help="genome")
    optparser.add_option("-t", "--treats", action="append", help="list of paths to treatment bw files")
    optparser.add_option("-x", "--xml", help="xml template to use")
    optparser.add_option("-o", "--out", help="output file")
    optparser.add_option("-l", "--treatIsLocal", help="treatment.bw are in local dir", action="store_true")    
    (options, args) = optparser.parse_args(sys.argv)

    if not options.genome or not options.treats or not options.out:
        optparser.print_help()
        sys.exit(-1)

    if options.xml:
        #USE user-specified template
        _template = options.xml
    else:
        #OTHERWISE: use default template
        _template = "cidc_chips/static/chips_igv.session.xml"

    if not os.path.isfile(_template):
        print("ERROR: Unable to find template XML file OR invalid file")
        sys.exit(-1)

    tree = ET.parse(_template)
    root = tree.getroot()
    
    #SET the genome
    root.attrib['genome'] = options.genome

    #GET the two important elements: Resources and (the first) Panel
    resources = root.findall('Resources')[0]
    panel = root.findall('Panel')[0]
    
    #TEMPLATE track
    track = ET.Element('Track')
    track.attrib = _attribs
    dr = ET.SubElement(track, 'DataRange')
    tmp = {'baseline':"0.0", 'drawBaseline':"true", 'flipAxis':"false",
           'maximum':"0.5331765", 'minimum':"0.0", 'type':"LINEAR"}
    dr.attrib = tmp

    for t in sorted(options.treats):
        #print(t)
        name = t.split("/")[-1] #filename
        #KEY: IGV allows for relative paths!
        relative_path = "/".join(t.split('/')[2:])

        #ADD new resource .bw file
        res = ET.Element('Resource') #singular of Resources
        res.attrib['path'] = name if options.treatIsLocal else relative_path 
        resources.append(res)
        
        #ADD new track to panel
        new_track = copy.deepcopy(track)
        #fixed attribs
    
        new_track.attrib['id'] = name if options.treatIsLocal else relative_path
        new_track.attrib['name'] = name
        
        #create sub-elm: DataRange
        panel.append(new_track)

    #OUTPUT
    tree = ET.ElementTree(root)
    tree.write(options.out)


if __name__=='__main__':
    main()


