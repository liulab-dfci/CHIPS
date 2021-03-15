import os
import subprocess
from snakemake.report import data_uri_from_file
import json
import pandas as pd

def parse_mapped_rate(sample):
    file = output_path + "/align/%s/%s_mapping.txt" %(sample,sample)
    f = open(file)
    total = int(f.readline().strip().split()[0])
    #skip 3 lines
    l = f.readline()
    l = f.readline()
    l = f.readline()
    mapped = int(f.readline().strip().split()[0])
    #skip 8 lines
    l = f.readline()
    l = f.readline()
    l = f.readline()
    l = f.readline()
    l = f.readline()
    l = f.readline()
    l = f.readline()
    l = f.readline()
    uniq_mapped = int(f.readline().strip())
    f.close()
    return total,mapped,uniq_mapped

def parse_GC_rate(sample):
    file = output_path + "/fastqc/%s/%s_stats.csv" %(sample,sample)
    f = open(file)
    MQ = int(f.readline().strip().split(",")[1])
    GC = int(f.readline().strip().split(",")[1])
    f.close()
    return MQ,GC

def parse_contamination(sample):
    file = output_path + "/contam/%s/%s_contamination.txt"%(sample,sample)
    f = open(file)
    percentage = [item.strip().split() for item in f.readlines()]
    f.close()
    return percentage

def parse_PBC(sample):
    file = output_path + "/frips/%s/%s_pbc.txt"%(sample,sample)
    f = open(file)
    firstLine = f.readline().strip().split()
    N1 = int(firstLine[1]) 
    Nd = N1
    for l in f:
        tmp = l.strip().split()
        Nd += int(tmp[1])
    f.close()
    return N1,Nd

def parse_FRiP(runRep):
    file = output_path + "/frips/%s/%s_frip.txt"%(runRep,runRep)
    f = open(file)
    ReadsInPeaks = int(f.readline().strip().split()[1])
    Total = int(f.readline().strip().split()[1])
    f.close()
    return ReadsInPeaks,Total

def parse_peaks(runRep):
    try:
        f = open(output_path+"/peaks/%s/%s_peaks.narrowPeak"%(runRep,runRep))
    except FileNotFoundError:
        f = open(output_path+"/peaks/%s/%s_peaks.broadPeak"%(runRep,runRep))
    #start the counts
    tot = fc_10 = fc_20 = 0
    for l in f:
        tmp = l.strip().split("\t") 
        #note FC is 7th col
        fc = float(tmp[6])
        if fc >= 20.0:
            fc_20 += 1
        if fc >= 10.0:
            fc_10 += 1
        tot += 1
    f.close()
    with open(output_path+"/ceas/%s/%s_summary.txt"%(runRep,runRep),"r") as ceas_meta:
        ceas_meta=(ceas_meta.read().replace("\'","\""))
        ceas = json.loads(ceas_meta)
        prom = ceas['Promoter']
        exon = ceas['Exon']
        intr = ceas['Intron']
        inte = ceas['Intergenic']
    with open(output_path+"/ceas/%s/%s_DHS_summary.dhs"%(runRep,runRep),"r") as dhs_meta:
        dhs_list = dhs_meta.readline().strip().split(",")
        if dhs_list[0] and int(dhs_list[0]) > 0:
            dhs="%.2f%%" % (int(dhs_list[1])*100/int(dhs_list[0]))
        else:
            dhs="%.2f%%" % 0.0
    return tot,fc_10,fc_20,dhs,prom,exon,intr,inte

def parse_targets(runRep):
    file = output_path +"/targets/%s/%s_gene_score.txt" % (runRep,runRep)
    table = pd.read_csv(file,comment = "#",sep = "\t",header = None)
    table = table.iloc[:,[0,1,2,4,6]]
    table.columns = ["chr","start","end","score","gene"]
    table = table.drop_duplicates().head(2000)
    table.index = range(0,2000)
    coordinate = []
    for i in range(len(table)):
        tmp = table.loc[i,["chr","start","end"]].values.tolist()
        coordinate.append("%s:%s-%s" % (tmp[0],tmp[1],tmp[2]))
    table["coordinate"] = coordinate
    table = table.loc[:,["gene","score","coordinate"]]
    html = "<tr><td>{Gene}</td><td>{Score:.2f}</td><td>{Coordinate}</td></tr>"
    txt = ""
    for i in range(table.shape[0]):
        gene = table.loc[i,"gene"]
        score = table.loc[i,"score"]
        coor = table.loc[i,"coordinate"]
        txt += html.format(Gene=gene, Score=score, Coordinate=coor)
    return txt

def parse_version(software):
    if software == "bwa":
        # out = subprocess.Popen("bwa", shell=True, stderr=subprocess.PIPE).stderr.read().decode('utf-8')
        # version = out.strip().split("\n")[1].split()[1]
        with open(output_path+"/align/run_info.txt") as out:
            version = out.read().strip()
    # elif software == "bowtie2":
    #     out = subprocess.check_output("bowtie2 --version", shell=True).decode('utf-8')
    #     version = out.strip().split("\n")[0].split("version")[1].strip()
    elif software == "macs2":
        with open(output_path+"/peaks/run_info.txt") as out:
            version = out.readlines()[0].strip().split()[1]
    else:
        version = ""
    return version


def result_dict(wildcards):
    report_dict={}
    report_dict["Config"]={}
    # ChipsVersion
    git_commit_string = "XXXXXX"
    if os.path.exists("cidc_chips/.git"):
        git_commit_string = subprocess.check_output('git --git-dir="cidc_chips/.git" rev-parse --short HEAD',shell=True).decode('utf-8').strip()
    report_dict["Config"]["ChipsVersion"]=git_commit_string
    # result path
    report_dict["Config"]["ResultsPath"]=os.path.abspath(output_path)
    # assembly
    report_dict["Config"]["AssemblyVersion"]=config["assembly"]
    # sentieon
    UsingSentieon = "No"
    if "sentieon" in config and config["sentieon"]:
        UsingSentieon = "Yes"
    report_dict["Config"]["UsingSentieon"]=UsingSentieon
    # aligner
    report_dict["Config"]["Aligner"]="bwa " + parse_version("bwa")
    # Peaks caller
    report_dict["Config"]["PeaksCaller"]="MACS2" + " " + parse_version("macs2")
    # cutoff of filtering 
    report_dict["Config"]["Cutoff"]=config["cutoff"]
    # motif
    MotifFinder = "No"
    report_dict["Motifs"]={}
    if "motif" in config and config["motif"]:
        MotifFinder = config["motif"]
        for run in config["runs"].keys():
            report_dict["Motifs"][run]={}
            for rep in _reps[run]:
                runRep = "%s.%s" % (run, rep)
                if config["motif"] == "mdseqpos":
                    report_dict["Motifs"][run][runRep]={"MotifHtml": "motif/%s/results/table.html"% runRep}
                if config["motif"] == "homer":
                    report_dict["Motifs"][run][runRep]={"MotifHtml": "motif/%s/results/homerResults.html"% runRep}
    report_dict["Config"]["MotifFinder"]=MotifFinder
    # contamination
    ContaminationPanel = "None"
    report_dict["Contam"]={}
    if "contamination_panel" in config and len(config["contamination_panel"]) > 0:
        ContaminationPanel = "; ".join([c.split("/")[-1] for c in config["contamination_panel"]])
        for run in config["runs"].keys():
            report_dict["Contam"][run] = {}
            for sample in config["runs"][run]:
                if sample:
                    report_dict["Contam"][run][sample]={}
                    for i in parse_contamination(sample):
                        report_dict["Contam"][run][sample][i[0]]="%.2f%%" % float(i[1])
    report_dict["Config"]["ContaminationPanel"]=ContaminationPanel
    # CNV
    CNVAnalysis = "No"
    if "cnv_analysis" in config and config["cnv_analysis"]:
        CNVAnalysis = "Yes"
    report_dict["Config"]["CNVAnalysis"]=CNVAnalysis
    # Run information
    RunsInformation=[]
    for run in config["runs"].keys():
        RunsInformation.append("%s: %s" % (run, ", ".join([s for s in config["runs"][run] if s])))
    report_dict["Config"]["Runlist"]=list(config["runs"].keys())
    report_dict["Config"]["RunsInformation"]="; ".join(RunsInformation)
    # Mapped GC PBC FRiP Fragments Conserv
    report_dict["Mapped"]={}
    report_dict["GC"]={}
    report_dict["PBC"]={}
    report_dict["FRiP"]={}
    report_dict["Fragments"]={}
    report_dict["Conserv"]={}
    report_dict["Peaks"]={}
    report_dict["Targets"]={}
    for run in config["runs"].keys():
        report_dict["Mapped"][run]={}
        report_dict["GC"][run]={}
        report_dict["PBC"][run]={}
        report_dict["Fragments"][run]={}
        report_dict["Conserv"][run]={}
        report_dict["FRiP"][run]={}
        report_dict["Peaks"][run]={}
        report_dict["Targets"][run]={}
        for rep in _reps[run]:
            runRep = "%s.%s" % (run, rep)
            ReadsInPeaks,Total=parse_FRiP(runRep)
            report_dict["FRiP"][run][runRep]={"ReadsInPeaks":ReadsInPeaks, "FRiP":"%.2f%%" % (ReadsInPeaks*100/Total), 
                                              "FRiPFigure": data_uri_from_file(output_path+"/frips/%s/%s_frip.png"%(runRep,runRep))[0]}
            report_dict["Conserv"][run][runRep]={"ConservationFigure": data_uri_from_file(output_path+"/conserv/%s/%s_conserv.png"%(runRep,runRep))[0]}
            tot,fc_10,fc_20,dhs,prom,exon,intr,inte=parse_peaks(runRep)
            report_dict["Peaks"][run][runRep]={"TotalPeaks":tot, "10FoldChangePeaks":fc_10, "20FoldChangePeaks":fc_20, "DHSPeaks":dhs, 
                                               "PromoterPeaks":prom, "ExonPeaks":exon, "IntronPeaks":intr, "IntergenicPeaks":inte, 
                                               "PeaksFigure":data_uri_from_file(output_path+"/peaks/%s/%s_peaks.png"%(runRep,runRep))[0],
                                               "PeaksFile":"files/%s_sorted_peaks.bed" % runRep}
            html_table = parse_targets(runRep)
            report_dict["Targets"][run][runRep]={"TargetsTable": html_table, "TargetsFile": "files/%s_gene_score.txt" % runRep}
        for sample in config["runs"][run]:
            if sample:
                total,mapped,uniq_mapped = parse_mapped_rate(sample)
                report_dict["Mapped"][run][sample]={"TotalReads":total, "MappedReads":mapped, "UniquelyMappedReads":uniq_mapped, 
                                                    "MappedRate":"%.2f%%" % (mapped*100/total), "UniquelyMappedRate":"%.2f%%" % (uniq_mapped*100/total),
                                                    "MappedFigure":data_uri_from_file(output_path + "/align/%s/%s_mapping.png"% (sample,sample))[0]}
                MQ,GC=parse_GC_rate(sample)
                report_dict["GC"][run][sample]={"MedianQuality":MQ, "GCMedian":GC, "GCFigure": data_uri_from_file(output_path + "/fastqc/%s/%s_perSeqGC.png"% (sample,sample))[0]}
                N1,Nd=parse_PBC(sample)
                report_dict["PBC"][run][sample]={"N1Score":N1,"NdScore":Nd,"PBCScore":"%.2f%%" % (N1*100/Nd)}
                if len(config["samples"][sample]) > 1:
                    report_dict["Fragments"][run][sample]={"FragmentFigure": data_uri_from_file(output_path + "/frag/%s/%s_fragDist.png"% (sample,sample))[0]}
    report_dict["Group"]={}
    report_dict["Group"]["Mapping"]=data_uri_from_file(output_path + "/report/image/mapping.png")[0]
    report_dict["Group"]["PBC"]=data_uri_from_file(output_path + "/report/image/pbc.png")[0]
    report_dict["Group"]["Peaks"]=data_uri_from_file(output_path + "/report/image/peakFoldChange.png")[0]
    return report_dict


def getReportInputs(wildcards):
    ret = []
    ret.append(output_path + '/report/image/mapping.png')
    ret.append(output_path + '/report/image/pbc.png')
    ret.append(output_path + '/report/image/peakFoldChange.png')
    ret.append(output_path + '/peaks/run_info.txt')
    for run in config["runs"].keys():
        for rep in _reps[run]:
            runRep = "%s.%s" % (run, rep)
            if "motif" in config and config["motif"]:
                if config["motif"] == "mdseqpos":
                    ret.append(output_path + "/motif/%s/results/table.html"% runRep)
                if config["motif"] == "homer":
                    ret.append(output_path + "/motif/%s/results/homerResults.html"% runRep)
            ret.append(output_path + "/frips/%s/%s_frip.txt"%(runRep,runRep))
            ret.append(output_path+"/frips/%s/%s_frip.png"%(runRep,runRep))
            if "macs2_broadpeaks" in config and config["macs2_broadpeaks"]:
                ret.append(output_path+"/peaks/%s/%s_peaks.broadPeak"%(runRep,runRep))
            else:
                ret.append(output_path+"/peaks/%s/%s_peaks.narrowPeak"%(runRep,runRep))
            ret.append(output_path+"/ceas/%s/%s_summary.txt"%(runRep,runRep))
            ret.append(output_path+"/ceas/%s/%s_DHS_summary.dhs"%(runRep,runRep))
            ret.append(output_path +"/targets/%s/%s_gene_score.txt"%(runRep,runRep))
            ret.append(output_path+"/conserv/%s/%s_conserv.png"%(runRep,runRep))
            ret.append(output_path + "/targets/%s/%s_gene_score.txt" % (runRep,runRep))
            ret.append(output_path + "/peaks/%s/%s_sorted_peaks.bed" % (runRep,runRep))
        for sample in config["runs"][run]:
            if sample:
                ret.append(output_path + "/align/%s/%s_mapping.txt" %(sample,sample))
                ret.append(output_path + "/frips/%s/%s_pbc.txt"%(sample,sample))
                ret.append(output_path + "/align/%s/%s_mapping.png"% (sample,sample))
                ret.append(output_path + "/fastqc/%s/%s_perSeqGC.png"% (sample,sample))
                ret.append(output_path + "/frag/%s/%s_fragDist.png"% (sample,sample))
                if "contamination_panel" in config and len(config["contamination_panel"]) > 0:
                    ret.append(output_path + "/contam/%s/%s_contamination.txt"% (sample,sample))
    return ret

def getZipReportInput(wildcards):
    ret = []
    for run in config["runs"].keys():
        for rep in _reps[run]:
            runRep = "%s.%s" % (run, rep)
            ret.append(output_path + "/report/files/%s_gene_score.txt" % runRep)
            ret.append(output_path + "/report/files/%s_sorted_peaks.bed" % runRep)
            if "motif" in config and config["motif"]:
                if config["motif"] == "mdseqpos":
                    ret.append(output_path + "/motif/%s/results/table.html"% runRep)
                if config["motif"] == "homer":
                    ret.append(output_path + "/motif/%s/results/homerResults.html"% runRep)
    return ret


def report_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    ls.append(output_path + '/report/report.html')
    ls.append(output_path + '/report/image/mapping.png')
    ls.append(output_path + '/report/image/pbc.png')
    ls.append(output_path + '/report/image/peakFoldChange.png')
    for run in config["runs"].keys():
        for rep in _reps[run]:
            runRep = "%s.%s" % (run, rep)
            ls.append(output_path + "/report/files/%s_gene_score.txt" % runRep)
            ls.append(output_path + "/report/files/%s_sorted_peaks.bed" % runRep)
    if ("macs2_broadpeaks" not in config) or config["macs2_broadpeaks"] != True:
        if 'motif' in config and config['motif']:
            ls.append(output_path + '/report/motif')
    ls.append(output_path + '/report/report.zip')
    return ls


rule report_all:
    input:
        report_targets

rule report_generate:
    input:
        getReportInputs
    output:
        output_path + '/report/report.html'
    message: "REPORT: Generate report for whole runs"
    run:
        report_dict=result_dict(wildcards)
        template = "cidc_chips/static/chipsTemplate.html"
        report = open(template)
        with open(str(output),"w") as o:
            o.write(report.read().replace("{ RESULT_DICT }",json.dumps(report_dict)))
        report.close()

if ("macs2_broadpeaks" not in config) or config["macs2_broadpeaks"] != True:
    if 'motif' in config and config['motif']:
        rule report_copyMotif:
            input:
                motif_targets
            output:
                directory(output_path + '/report/motif')
            message: "REPORT: Copy motif results to report"
            params:
                motif_path = output_path + "/motif/"
            shell:
                "cp -r {params.motif_path} {output}"

rule report_copyTargets:
    input:
        output_path + "/targets/{run}.{rep}/{run}.{rep}_gene_score.txt"
    output:
        score=output_path + "/report/files/{run}.{rep}_gene_score.txt",
    message: "REPORT: Copy targets results to report"
    shell:
        "cp {input} {output}"

rule report_copyPeaks:
    input:
        output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.bed"
    output:
        output_path + "/report/files/{run}.{rep}_sorted_peaks.bed",
    message: "REPORT: Copy peaks results to report"
    shell:
        "cp {input} {output}"

rule report_plotMapStat:
    input:
        output_path + "/align/mapping.csv"
    output:
        output_path + "/report/image/mapping.png"
    log: _logfile
    conda: "../envs/report/report.yaml"
    shell:
        "Rscript cidc_chips/modules/scripts/report_mapStats.R {input} {output}"

rule report_plotPBCStat:
    input:
        #output_path + "/align/pbc.csv"
        output_path + "/frips/pbc.csv"
    output:
        output_path + "/report/image/pbc.png"
    log: _logfile
    conda: "../envs/report/report.yaml"
    shell:
        "Rscript cidc_chips/modules/scripts/report_plotPBC.R {input} {output}"

rule report_plotPeakFoldChange:
    input: 
        output_path + "/peaks/peakStats.csv"
    output:
        output_path + "/report/image/peakFoldChange.png"
    log: _logfile
    conda: "../envs/report/report.yaml"
    shell:
        "Rscript cidc_chips/modules/scripts/report_plotFoldChange.R {input} {output}"

rule report_zipReport:
    input:
        report=output_path + '/report/report.html',
        files=getZipReportInput,
    output:
        output_path + '/report/report.zip'
    message: "REPORT: Compress report"
    params: 
        motif = output_path + '/report/motif' if ("motif" in config and config["motif"]) else "" ,
        files = output_path + '/report/files'
    shell:
        "zip -q -r {output} {input.report} {params.files} {params.motif}"















