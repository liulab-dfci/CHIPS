#!/usr/bin/env python3
"""Len Taing 2020 (TGBTG)
Generalized script to generate analysis/report/report.html
based on the contents of analysis/report sub-dirs: each sub-dir is a 'section'
"""

import os
import sys

import re
import subprocess
import yaml
import json

import jinja2
import markdown
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

#Multiqc stuff
from multiqc.plots import bargraph, linegraph, table
from multiqc.utils import report as mqc_report, config as mqc_config
import multiqc.modules as mqc_modules


#Plotly stuff
import plotly
import plotly.express as px
import plotly.graph_objects as go

from optparse import OptionParser

_mqc_plot_types = {'bar':bargraph, 'line':linegraph, 'table':table}
_resources = {}

def parseYaml(yaml_file):
    """Parses a yaml file and returns a dictionary"""
    ret = {}
    with open(yaml_file, 'r') as stream:
        try:
            ret = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return ret

def prettyprint(s, toUpper=False):
    """Given a string, replaces underscores with spaces and uppercases the
    first letter of each word ONLY if the string is composed of lowercased
    letters.  If the param, toUpper is given then s.upper is returned.

    Examples: "data_quality" -> "Data Quality"
    "copy_number_123" -> "Copy Number 123"

    "My_own_title" -> "My own title"
    "Hla" -> "Hla"
    """
    if toUpper:
        s = s.upper()
        s= s.replace("_", " ")
    else:
        s = s.title()
        s= s.replace("_", " ")
    return s

def is_image_file(s):
    """Checks to see if the string starts with 'img:'"""
    return s.startswith('img:')

def is_color(s):
    """Checks to see if the string starts with 'red:' or hex html color"""
    valid_colors = ['red','green','blue']
    #reference for html hex color reg-ex
    #ref: https://stackoverflow.com/questions/30241375/python-how-to-check-if-string-is-a-hex-color-code
    prog = re.compile('^#(?:[0-9a-fA-F]{3}){1,2}$')

    if ":" in s:
        (prefix, val) = s.split(":")
        return (prefix in valid_colors) or prog.match(prefix)
    else:
        return False

def get_val(s):
    (prefix, val) = s.split(":")
    return val

def get_prefix(s):
    (prefix, val) = s.split(":")
    return prefix

#NOTE: this can easily handle csv files too!
def buildTable(tsv_file, details, jinjaEnv, cssClass=""):
    """Given a tsv file, and a section--converts the tsv file to a table
    assumes the first line is the hdr"""
    #TRY to auto-detect if plot file is , or \t separated
    f = open(tsv_file)
    tmp = f.readline() #take a peek at first line
    separator = "," if "," in tmp else "\t"
    f.close()

    #LOAD jinja2 test and filter for image processing
    jinjaEnv.tests['imagefile'] = is_image_file
    jinjaEnv.tests['color'] = is_color
    jinjaEnv.filters['get_value'] = get_val
    jinjaEnv.filters['get_prefix'] = get_prefix
    template = jinjaEnv.get_template("table.html")

    #Title is the filename prettyprinted
    fname = ".".join(tsv_file.split("/")[-1].split(".")[:-1])
    #REMOVE index from file name, e.g. 01_foo -> foo
    index = fname.split("_")[0] #first save index
    fname = "_".join(fname.split("_")[1:])
    path = "/".join(tsv_file.split("/")[:-1]) #drop the file
    title = prettyprint(fname, toUpper=True)

    vals = {'container': fname+'_container',
            'id': fname, 'title':title, 'class': cssClass}

    #Check for a caption
    caption = details.get('caption', None)
    if caption:
        vals['caption'] = renderMd(caption)
    #check for subcaption
    sub_caption = details.get('subcaption', None)
    if sub_caption:
        vals['sub_caption'] = renderMd(sub_caption)

    table = []
    f = open(tsv_file)
    hdr = f.readline().strip().split(separator)

    for l in f:
        tmp = l.strip().split(separator)
        table.append(tmp)
    f.close()

    vals['header'] = hdr
    vals['table'] = table
    #print(vals)
    return template.render(vals)

def buildPlot(png_file, details, jinjaEnv):
    """Given a png file displays the plot..simple!"""
    template = jinjaEnv.get_template("plot.html")
    fname = ".".join(png_file.split("/")[-1].split(".")[:-1])
    #REMOVE index from file name, e.g. 01_foo -> foo
    index = fname.split("_")[0] #first save index
    fname = "_".join(fname.split("_")[1:])
    path = "/".join(png_file.split("/")[:-1]) #drop the file
    title = prettyprint(fname,toUpper=True)

    #make the png file path REALITIVE to the report.html file!
    png_file_relative = "/".join(png_file.split("/")[2:])
    vals = {'id': fname,
            'title':title,
            'png_file': png_file_relative
    }
    #Check for a caption
    caption = details.get('caption', None)
    if caption:
        vals['caption'] = renderMd(caption)
    #check for subcaption
    sub_caption = details.get('subcaption', None)
    if sub_caption:
        vals['sub_caption'] = renderMd(sub_caption)

    #print(vals)
    return template.render(vals)

def renderMd(s):
    """renders s as markdown text"""
    return  markdown.markdown(s, extensions=['extra'])

def readMqcData01(data_file, separator=","):
    """This file is like a typical csv, where the first line is
    the header and first col = Sample names
    returns a dictionary where the keys are sample names and the
    values are dictionaries that represent the columns
    ref: https://multiqc.info/docs/#bar-graphs
    """
    f = open(data_file)
    hdr = f.readline().strip().split(separator)
    #print(hdr)
    data = {}
    for l in f:
        #Assume col 1 = Samples
        tmp = l.strip().split(separator)
        data[tmp[0]] = dict(zip(hdr[1:],tmp[1:]))
    f.close()
    #print(data)
    return data

def readMqcData02(data_file, separator = ","):
    """This file is like the other csv representation
    the header = X, Sample1, Sample2, ..., SampleN
    Each row represents the x-val, and then y-vals for each sample at
    the x-val
    returns a dictionary where the keys are sample names and the
    values x-y value pairs
    ref: https://multiqc.info/docs/#line-graph
    """
    f = open(data_file)
    hdr = f.readline().strip().split(separator)
    #Init data to keys: [] for each header item
    data = dict([(h,[]) for h in hdr])
    for l in f:
        tmp = zip(hdr, l.strip().split(separator))
        for (k, v) in tmp:
            data[k].append(v)
    f.close()
    #Get x-axis and remove it from data
    #x = map(lambda x: x.zfill(3), data['X'])
    xaxis = list(map(lambda x: int(x), data['X']))
    del data['X']

    #Reprocess each sample
    for s in data:
        #try to infer type
        if '.' in data[s][0]: #float?
            data[s] = list(map(float, data[s])) #DON'T forget to list() map objs!
        else:
            data[s] = list(map(int, data[s]))
        tmp = dict(zip(xaxis, data[s]))
        #print(tmp)
        data[s] = tmp
        #print(data[s])
    #print(data)
    return data

def buildMqcPlot(plot_file, details, jinjaEnv):
    """Given a plotfile which is a csv file with a .plot extension,
    Tries to generate an (interactive) multiqc plot"""
    #TRY to auto-detect if plot file is , or \t separated
    f = open(plot_file)
    tmp = f.readline() #take a peek at first line
    sep = "," if "," in tmp else "\t"
    f.close()
    template = jinjaEnv.get_template("mqc_plot.html")
    #Dropping path and extension to get filename
    fname = ".".join(plot_file.split("/")[-1].split(".")[:-1])
    #REMOVE index from file name, e.g. 01_foo -> foo
    index = fname.split("_")[0] #first save index
    plot_type = fname.split("_")[-1] #and plot type
    fname = "_".join(fname.split("_")[1:-1])
    title = prettyprint(fname,toUpper=True)

    #check to see if we're adding a module or a plot
    if plot_type in _mqc_plot_types:
        #Pickout the proper mqc plot module to use
        mqc_plot = _mqc_plot_types[plot_type]

        #READ in file--for now assume it's of type bar and have other handlers
        #later
        #HERE we should check for plot type
        if plot_type == "bar" or plot_type == "table":
            data = readMqcData01(plot_file, sep)
            #Try to get plot details from details dict
            cats = details.get("cats", None) #Categories
            #pass the rest of details into the plotting fn
            html_plot = mqc_plot.plot(data, cats, details)
        else: #line
            data = readMqcData02(plot_file, sep)
            html_plot = mqc_plot.plot(data, details)



        vals = {'id': fname,
                'title':title,
                'plot': html_plot,
        }

        #Check for a caption
        caption = details.get('caption', None)
        if caption:
            vals['caption'] = renderMd(caption)
        #check for subcaption
        sub_caption = details.get('subcaption', None)
        if sub_caption:
            vals['sub_caption'] = renderMd(sub_caption)

        #print(vals)
        ret = template.render(vals)

    else: #TRY incoporate MQC module-- plot_type should define which one
        #NOTE: the plot_file itself is empty
        #print(dir(mqc_modules))
        #print(mqc_modules.__file__)
        mod = __import__("multiqc.modules.%s" % plot_type)
        klass = getattr(getattr(mqc_modules, plot_type), "MultiqcModule")
        results_dir = details.get("results_dir", "./analysis")
        #initialize where to search for MQC files
        if results_dir not in mqc_config.analysis_dir:
            mqc_config.analysis_dir.append(results_dir)
        #initialize report.get_filelist
        mqc_report.get_filelist([plot_type])
        #TRY to render out each section
        mqcKlass = klass()
        ret = ""
        for s in mqcKlass.sections:
            vals = {'id': s['name'],
                    'title': s['name'],
                    'plot': s['plot'],
                    'caption': s['description']}
            ret += template.render(vals)

    return ret

def formatter(x, pos):
    'The two args are the value and tick position'
    return '%1.1fM' % (float(x)*1e-6)

#SOON to be deprecated!
def plot(plot_f):
    """Given a path to a csv file, generate a plot with the same name
    NOTE: the last _ of the name gives the type of plot"""
    #PARSE the filename
    tmp = plot_f.split(".")[0].split("_")
    name = "_".join(tmp[1:-1])

    df = pd.read_csv(plot_f)
    plot = getattr(df.plot, tmp[-1])
    plot(x="Sample")
    #df.plot.barh(x="Sample")
    #plt.xlabel("Reads") # Text for Y-Axis
    #plt.ylabel("Sample") # Text for Y-Axis
    #plt.title("Mapping statistics plot")
    ax = plt.gca()
    ax.xaxis.set_major_formatter(FuncFormatter(formatter))
    plt.savefig('%s.png' % "_".join(tmp))

def loadJson(json_file):
    """"Simply loads the json data into _resources which will be added
    into the html directly for use"""
    #get the filename, which will be the key in _resources
    fname = ".".join(json_file.split("/")[-1].split(".")[:-1])
    ffile = open(json_file)
    tmp = json.load(ffile)
    ffile.close()
    _resources[fname] = tmp

def buildPlotly(plotly_file, details, jinjaEnv):
    """Given a plotfile which is a csv file with a .plotly extension,
    Tries to generate an (interactive) plotly chart"""
    #TRY to auto-detect if plot file is , or \t separated
    f = open(plotly_file)
    tmp = f.readline() #take a peek at first line
    sep = "," if "," in tmp else "\t"
    f.close()

    #NOTE: it's simple enough to just use the same form
    template = jinjaEnv.get_template("mqc_plot.html")
    #Dropping path and extension to get filename
    fname = ".".join(plotly_file.split("/")[-1].split(".")[:-1])
    #REMOVE index from file name, e.g. 01_foo -> foo
    index = fname.split("_")[0] #first save index
    plot_type = fname.split("_")[-1] #and plot type
    fname = "_".join(fname.split("_")[1:-1])
    title = prettyprint(fname, toUpper=True)

    #Pickout the proper mqc plot module to use
    plot = getattr(px, plot_type) if plot_type != "oncoplot" else oncoplot
    df = pd.read_csv(plotly_file, index_col=0)
    #Colors: red, green, blue, purple, gray, gold
    colors = ['#e84118','#44bd32','#0097e6', '#8c7ae6', '#7f8fa6', '#e1b12c']
    #print(details)
    if 'plotly' in details:
        if not 'color_discrete_sequence' in details['plotly']:
            details['plotly']['color_discrete_sequence'] = colors
    else:
        details['plotly'] = {'color_discrete_sequence':colors}
    fig = plot(df, **details['plotly'])
    fig.update_layout(plot_bgcolor='#f6f6f6')
    html_plot = plotly.offline.plot(fig, include_plotlyjs=False, output_type='div')

    vals = {'id': fname,
            'title':title,
            'plot': html_plot,
    }

    #Check for a caption
    caption = details.get('caption', None)
    if caption:
        vals['caption'] = renderMd(caption)
    #check for subcaption
    sub_caption = details.get('subcaption', None)
    if sub_caption:
        vals['sub_caption'] = renderMd(sub_caption)

    #print(vals)
    return template.render(vals)

def oncoplot(df, **kwargs):
    """Given a dataframe returns a plotly oncoplot figure"""
    #print(kwargs)

    #Sort samples by hits
    #ref: https://stackoverflow.com/questions/20480238/getting-top-3-rows-that-have-biggest-sum-of-columns-in-pandas-dataframe
    top_ngenes = min(kwargs['top_ngenes'], len(df.columns))
    sample_idx = df.sum(axis=1).sort_values(ascending=False).index
    df = df.reindex(sample_idx)

    #Sort by genes/annotations by hits
    #top_genes = df.sum(axis=0).sort_values(ascending=False).tail(top_ngenes)
    top_genes = df.sum(axis=0).sort_values(ascending=False).head(top_ngenes)
    #print(top_genes)
    df = df.reindex(top_genes.index, axis=1)

    #moving subplots closer-
    #ref: https://stackoverflow.com/questions/31526045/remove-space-between-subplots-in-plotly
    fig = plotly.subplots.make_subplots(rows=1, cols=2, column_widths=[0.9, 0.1], specs = [[{}, {}]], horizontal_spacing=0.005)
    #Colors- 0 = lightgray, 1 = green
    colors = kwargs.get('colors', ['#b5b5b5', '#44bd32'])
    fig.add_trace(go.Heatmap(x=df.index,
                             y=list(reversed(top_genes.index)),
                             z = df.transpose().iloc[::-1], #KEY: NEED to transpose matrix and reverse rows
                             type = 'heatmap',
                             colorscale = colors,
                             xgap = 3, ygap=3, #Add gridlines
                             showscale=False, #turn off scale
                             hoverinfo=["x",'y'], #turn off hover data
    ), row=1, col=1)
    fig.add_trace(go.Bar(x=list(reversed(top_genes)),
                         orientation='h',
                         texttemplate="%{x}",
                         textposition='inside',
                         textfont_color="white",
                         marker_color= "#0097e6", #blue
                         hoverinfo="skip"), row=1, col=2)
    fig.update_xaxes(showticklabels=False, row=1, col=1) #turn off samplenames
    fig.update_yaxes(showticklabels=False, row=1, col=2) #turn off percents
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)')
    #fig.write_html("wes_oncoplot.html")
    return fig

def main():
    usage = "USAGE: %prog -o [output html file]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-d", "--dir", help="report directory path")
    optparser.add_option("-s", "--sections", help="sections list in order of appearance")
    optparser.add_option("-j", "--template", help="path to jinja2 template")
    optparser.add_option("-t", "--title", help="report title", default="WES Summary Report")
    optparser.add_option("-o", "--output", help="output html file")
    (options, args) = optparser.parse_args(sys.argv)

    if not options.dir or not options.sections or not options.template or not options.output:
        optparser.print_help()
        sys.exit(-1)

    #SET report title-default to "WES Summary Report"
    title = options.title

    #ASSUMING it's being run as WES project level
    template_dir = os.path.dirname(options.template)
    template_fname = os.path.basename(options.template)
    templateLoader = jinja2.FileSystemLoader(searchpath=template_dir)
    templateEnv = jinja2.Environment(loader=templateLoader)
    template = templateEnv.get_template(template_fname)

    #Build up this dictionary
    wes_report = {'title':title}

    #infer sections from the analysis/report dir structure
    #sections = os.listdir(options.dir)
    sections = options.sections.split(",")
    wes_panels = {}
    first_section = ""
    #ONLY sections so far--no further recursion

    #for sect in os.listdir(options.dir):
    for (i, sect) in enumerate(sections):
        if sect == "static": #SKIP static content if it's there
            continue
        #Check for {section}.yaml file for overrides
        path = os.path.join(options.dir, sect)
        ordering = sorted(os.listdir(path))

        #Check for meta files- runs_meta.json and samples_meta.json
        if os.path.exists(os.path.join(options.dir, 'runs_meta.json')):
            loadJson(os.path.join(options.dir, 'runs_meta.json'))
        if os.path.exists(os.path.join(options.dir, 'samples_meta.json')):
            loadJson(os.path.join(options.dir, 'samples_meta.json'))

        #Build container
        if i == 0: #first element is shown
            first_section = sect
            tmp = """<div id="%s" class="container wes_container">\n""" % sect
        else:
            tmp = """<div id="%s" class="container wes_container" style="display:none">\n""" % sect

        for ffile in ordering:
            filepath = os.path.join(path, ffile)
            #I will soon deprecate these 'plots' sub-dir--is not currenlty used by either report
            if ffile == 'plots' and os.path.isdir(filepath): #handle plots
                #PLOT all csv files
                csv_files = [a for a in os.listdir(filepath) if a.endswith('.csv')]
                for f in csv_files:
                    plot_f = os.path.join(filepath, f)
                    plot(plot_f)

            elif os.path.isfile(filepath): #SKIP directories
                #CHECK for file existance
                if not os.path.exists(filepath):
                    print("WARNING: report.py- file %s is not found. SKIPPED for rendering" % filepath)
                    continue
                #check for details
                index = ffile.split("_")[0] #Get part index
                if os.path.exists(os.path.join(path, "%s_details.yaml" % index)):
                    details = parseYaml(os.path.join(path, "%s_details.yaml" % index))
                else:
                    details = {}

                if ffile.endswith(".tsv") or ffile.endswith(".csv"): #MAKE a table
                    tmp += buildTable(filepath, details, templateEnv)
                elif ffile.endswith(".dt"):
                    tmp += buildTable(filepath, details, templateEnv, "wes_datatable")
                elif ffile.endswith(".png"): #Make a plot
                    tmp += buildPlot(filepath, details, templateEnv)
                elif ffile.endswith(".mqc"): #Make a Multiqc plot
                    tmp += buildMqcPlot(filepath, details, templateEnv)
                elif ffile.endswith(".json"): #Load json data
                    loadJson(filepath)
                elif ffile.endswith(".plotly"): #Make a plotly plot
                    tmp += buildPlotly(filepath, details, templateEnv)

        #END container
        tmp += "\n</div>"
        wes_panels[sect] = tmp

    wes_sections = [(s, prettyprint(s)) for s in sections]
    wes_report['sections'] = wes_sections
    wes_report['panels'] = wes_panels
    wes_report['first_section'] = first_section
    wes_report['plot_compressed_json'] = mqc_report.compress_json(mqc_report.plot_data)
    wes_report['wes_resources'] = json.dumps(_resources)
    wes_report['config'] = mqc_config
    template.stream(wes_report).dump(options.output)

if __name__ == '__main__':
    main()
