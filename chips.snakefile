#!/usr/bin/env python

import os
import sys
import subprocess
import pandas as pd
import yaml
import re
import errno

from string import Template

def getRuns(config):
    """parse metasheet for Run groupings"""
    ret = {}

    #LEN: Weird, but using pandas to handle the comments in the file
    #KEY: need skipinitialspace to make it fault tolerant to spaces!
    metadata = pd.read_csv(config['metasheet'], index_col=0, sep=',', comment='#', skipinitialspace=True)
    f = metadata.to_csv().split() #make it resemble an actual file with lines
    #SKIP the hdr
    for l in f[1:]:
        tmp = l.strip().split(",")
        ret[tmp[0]] = tmp[1:]
    config['runs'] = ret
    return config

# def addPy2Paths_Config(config):
#     """ADDS the python2 paths to config"""
#     conda_root = subprocess.check_output('conda info --root',shell=True).decode('utf-8').strip()
#     conda_path = os.path.join(conda_root, 'pkgs')
#     config["python2_pythonpath"] = os.path.join(conda_root, 'envs', 'chips_py2', 'lib', 'python2.7', 'site-packages')
#     if not "python2" in config or not config["python2"]:
#         config["python2"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'python2.7')
    # if not "mdseqpos_path" in config or not config["mdseqpos_path"]:
    #     config["mdseqpos_path"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'MDSeqPos.py')
    # if not "macs2_path" in config or not config["macs2_path"]:
    #     config["macs2_path"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'macs2')

def check_bwa_index_exist(path):
    """check if bwa index files exist or not
    given a path to the fasta.
    e.g., given ./ref_files/hg38/bwa_indices/hg38/hg38.fa
    check if
    ./ref_files/hg38/bwa_indices/hg38/hg38.fa.amb
    ./ref_files/hg38/bwa_indices/hg38/hg38.fa.ann
    ./ref_files/hg38/bwa_indices/hg38/hg38.fa.bwt
    ./ref_files/hg38/bwa_indices/hg38/hg38.fa.pac
    ./ref_files/hg38/bwa_indices/hg38/hg38.fa.sa
    exist or not

    """
    bwa_suffix = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    bwa_index_files = [ path + suffix for suffix in bwa_suffix ]
    missing_index_files = []
    for file in bwa_index_files:
        if not os.path.isfile(file):
            missing_index_files.append(file)
    return missing_index_files


def loadRef(config):
    """Adds the static reference paths found in config['ref']
    NOTE: if the elm is already defined, then we DO NOT clobber the value

    Also check if reference files exist
    e.g. The config['assembly'] is hg38
    In the ref.yaml file under hg38:
    bwa_index, geneTable, geneBed, conservation, DHS, exons, promoters, velcro_regions
    and chrom_lens files should exist
    """
    f = open(config['ref'])
    ref_info = yaml.safe_load(f)
    f.close()

    missing_ref = []
    #print(ref_info[config['assembly']])
    if ref_info.get(config['assembly']):
        for (k,v) in ref_info[config['assembly']].items():
        #NO CLOBBERING what is user-defined!
            if k not in config:
                config[k] = v
            if k in ['geneTable', 'geneBed', 'conservation', 'DHS', 'exons', 'promoters', 'chrom_lens']:
                if not os.path.isfile(v):
                    missing_ref.append(v)
            elif k == "bwa_index":
                missing_ref.extend(check_bwa_index_exist(v))
    else:
        print("assembly {} specified in config.yaml file does not exist in ref.yaml file".format(config['assembly']))
        sys.exit(1)

    if config.get('contamination_panel_qc'):
        # check if contamination reference files exist, The bwa index files should exist
        for contamination in ref_info['contamination_panel']:
            missing_ref.extend(check_bwa_index_exist(contamination))
        config['contamination_panel'] = ref_info['contamination_panel']

    return missing_ref

def check_fastq_exist(config):
    """check if the fastq files listed in the config[samples]
    exist or not
    """
    missing_fqs = []
    samples = config['samples']
    for sample in samples.keys():
        for fq in samples[sample]:
            if not os.path.isfile(fq):
                missing_fqs.append(fq)
    return missing_fqs


#---------  CONFIG set up  ---------------
configfile: "config.yaml"   # This makes snakemake load up yaml into config
config = getRuns(config)
# addPy2Paths_Config(config)

#NOW load ref.yaml - SIDE-EFFECT: loadRef CHANGES config. Also returns a list of missing reference files
missing_refs = loadRef(config)
if missing_refs:
    for reference in missing_refs:
        print( "\n" + "ERROR!! file {} specified in the ref.yaml does not exist!".format(reference) + "\n")
    sys.exit(1)



# preflight check for fastqs exist or not
missing_fqs = check_fastq_exist(config)
if missing_fqs:
    for fq in missing_fqs:
        print( "\n" + "ERROR!! fastq file {} does not exist! make sure you have the right path.".format(fq) + "\n")
    sys.exit(1)


#-----------------------------------------

#------------------------------------------------------------------------------
# Handle replicates
#------------------------------------------------------------------------------
#used to define the replicate structure
_reps = {}
for run in config['runs'].keys():
    r = config['runs'][run]
    tmp = []
    for (rep, i) in enumerate(range(0, len(r), 2)):
        if r[i]: tmp.append("rep%s" % str(rep+1))
    _reps[run] = tmp
#print(_reps)

# Set output path
if ("output_path" not in config) or config["output_path"] == "":
    output_path = "analysis"
else:
    output_path = re.sub("^\./","", config["output_path"].rstrip("/"))

#NOTE: Template class allows for _ in the variable names, we want to DISALLOW
#that for replicates
#ref: http://stackoverflow.com/questions/2326757/string-templates-in-python-what-are-legal-characters
class RepTemplate(Template):
    idpattern = r'[a-z][a-z0-9]*'

#THIS helper fn is used in several of the modules peaks, ceas, frips
#Instead of an expand, we need this fn to create the CORRECT input-list
def _getRepInput(temp, suffix=""):
    """generalized input fn to get the replicate files
    CALLER passes in temp: a python string template that has the var runRep
    e.g. analysis/ceas/$runRep/$runRep_DHS_stats.txt
    Return: list of the string template filled with the correct runRep names
    """
    #print(temp)
    s = RepTemplate(temp)
    ls = []
    for run in config['runs'].keys():
        for rep in _reps[run]:
            #GENERATE Run name: concat the run and rep name
            runRep = "%s.%s" % (run, rep)
            ls.append(s.substitute(runRep=runRep,))
    #print(ls)
    return ls

#------------------------------------------------------------------------------
# TARGETS
#------------------------------------------------------------------------------

def all_targets(wildcards):
    ls = []
    if config.get('trim_adapter'):
        ls.extend(trim_targets(wildcards))
    #IMPORT all of the module targets
    ls.extend(align_targets(wildcards))
    ls.extend(peaks_targets(wildcards))
    ls.extend(fastqc_targets(wildcards))
    ls.extend(conservation_targets(wildcards))
    ls.extend(ceas_targets(wildcards))
    ls.extend(frips_targets(wildcards))
    ls.extend(targets_targets(wildcards))
    #Check to see if motif is enabled
    if ("macs2_broadpeaks" not in config) or config["macs2_broadpeaks"] == False:
        if 'motif' in config:
            ls.extend(motif_targets(wildcards))

    #HANDLE CNV/qdnaseq analysis
    _qdnaseq = config["cnv_analysis"]
    if _qdnaseq:
        #ls.extend(qdnaseq_targets(wildcards))
        #check for some inputs
        hasInput = False
        #HACK: for some reason, using the following line causes errors
        #for (run, ls) in config['runs'].items():
        #SO we call getRuns (from above) using a simplified config
        tmp_config = {'metasheet': config['metasheet']}
        runs = getRuns(tmp_config)['runs'].copy()
        for (run) in runs.keys():
            #NOTE: if i do this, this is an error!
            #ls = runs[run]
            if runs[run][1] or runs[run][3]:
                #these are the control sample indices
                hasInput = True
                break

        if hasInput:
            ls.extend(qdnaseq_targets(wildcards))

    if config.get('contamination_panel_qc'):
        ls.extend(contamination_targets(wildcards))
    # skip running modules that useless in cistrome db
    if config.get('CistromeApi'):
        ls.extend(json_targets(wildcards))
        ls.extend(cistrome_targets(wildcards))
    else:
        # ls.extend(mapmaker_targets(wildcards))
        # ls.extend(bam_snapshots_targets(wildcards))
        #ls.extend(report_targets(wildcards))
        if config.get('epicypher_analysis'):
            ls.extend(epicypher_targets(wildcards))
    ls.extend(checking_targets(wildcards))
    ls.append(output_path + "/report/report.html")
    #ls.extend(report_targets(wildcards))
    return ls


rule target:
    input:
        all_targets,

    message: "Compiling all output"
# if config['aligner'] == 'bwt2':
#     include: "./modules/align_bwt2.snakefile"     # rules specific to Bowtie2
# else:
if config.get('trim_adapter'):
    include: "./modules/trim_adapter.snakefile"
    include: "./modules/align_common.snakefile"
    include: "./modules/align_bwa_trim.snakefile"
else:
    include: "./modules/align_bwa.snakefile"      # rules specific to BWA
    include: "./modules/align_common.snakefile"  # common align rules

include: "./modules/peaks.snakefile"         # peak calling rules
include: "./modules/fastqc.snakefile"        # fastqc (sequence qual) rules
include: "./modules/conservation.snakefile"  # generate conservation plot
include: "./modules/ceas.snakefile"          # annotate peak regions
include: "./modules/frips.snakefile"         # fraction of reads in peaks

if ("macs2_broadpeaks" not in config) or config["macs2_broadpeaks"] != True:
    if 'motif' in config and config['motif'] == 'mdseqpos':
        include: "./modules/motif_mdseqpos.snakefile"     # mdseqpos motif module
    else:
        include: "./modules/motif_homer.snakefile"        # homer motif module

if config.get('contamination_panel_qc'):
    include: "./modules/contamination.snakefile" # contamination panel module
include: "./modules/qdnaseq.snakefile"       # qdnaseq (CNV) module
include: "./modules/mapmaker.snakefile"      # chips-mapmaker interface module
include: "./modules/epicypher.snakefile"     # epicypher spike-in module
include: "./modules/bam_snapshots.snakefile" # generate bam snapshots module
include: "./modules/targets.snakefile"       # targets module
#include: "./modules/report.snakefile"        # report module
include: "./modules/json.snakefile"          # json module
include: "./modules/cistrome.snakefile"      # cistrome adapter module
include: "./modules/emptychecking.snakefile" # checking empty file module
include: "./modules/new_report.snakefile"
