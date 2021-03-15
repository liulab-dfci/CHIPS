#MODULE: conservation- module to create conservation plots
import math

# _logfile=output_path + "/logs/conservation.log"
#_numPngs is used in conservation_plot rule to see how many pngs to expect
#note: the rule plots 3 runs per png, so for example, 12 runs results in 4 pngs
_nPerPlot = 3
_numPngs = math.ceil(len(config['runs'].keys())/float(_nPerPlot))
_nPngs = [n+1 for n in range(_numPngs)]

#NOTE: using the _refs from chips.snakefile
def conservation_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config["runs"].keys():
        for rep in _reps[run]:
            #GENERATE Run name: concat the run and rep name
            runRep = "%s.%s" % (run, rep)
            # ls.append(output_path + "/peaks/%s/%s_sorted_5k_summits.bed" % (runRep,runRep))
            if ("macs2_broadpeaks" in config) and config["macs2_broadpeaks"]:
                ls.append(output_path + "/peaks/%s/%s_sorted_5k_peaks.bed" % (runRep,runRep))
            else:
                ls.append(output_path + "/peaks/%s/%s_sorted_5k_summits.bed" % (runRep,runRep))
            ls.append(output_path + "/conserv/%s/%s_conserv.R" % (runRep,runRep))
            ls.append(output_path + "/conserv/%s/%s_conserv.png" % (runRep,runRep))
            ls.append(output_path + "/conserv/%s/%s_conserv_thumb.png" % (runRep,runRep))
    return ls

def conservationInput(wildcards):
    run = wildcards.run
    rep = wildcards.rep
    runRep = "%s.%s" % (run,rep)
    if ("macs2_broadpeaks" in config) and config["macs2_broadpeaks"]:
        temp = output_path + "/peaks/%s/%s_sorted_peaks.bed" % (runRep,runRep)
    else:
        temp = output_path + "/peaks/%s/%s_sorted_summits.bed" % (runRep,runRep)
    return temp

rule conservation_all:
    input:
        conservation_targets

rule conservation_top5kBroadPeaks:
    """take the top 5000 peaks, sorted by score"""
    input:
        output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.bed"
        # output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_summits.bed"
    output:
        output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_5k_peaks.bed"
    params:
        lines = 5000
    message: "CONSERVATION: top5k_peaks"
    log: output_path + "/logs/conservation/{run}.{rep}.log"
    benchmark: output_path + "/Benchmark/{run}.{rep}_conservation_top5kBroadPeaks.benchmark"
    conda: "../envs/conservation/conservation.yaml"
    shell:
        "head -n {params.lines} {input} > {output}"

rule conservation_top5kPeaks:
    """take the top 5000 peaks, sorted by score"""
    input:
        output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_summits.bed"
    output:
        output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_5k_summits.bed"
    params:
        lines = 5000
    message: "CONSERVATION: top5k_peaks"
    log: output_path + "/logs/conservation/{run}.{rep}.log"
    benchmark: output_path + "/Benchmark/{run}.{rep}_conservation_top5kPeaks.benchmark"
    conda: "../envs/conservation/conservation.yaml"
    shell:
        "head -n {params.lines} {input} > {output}"

rule conservation_plotConservation:
    """generate conservation plots"""
    input:
        # output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_summits.bed"
        # output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.bed"
        conservationInput
    output:
        png=output_path + "/conserv/{run}.{rep}/{run}.{rep}_conserv.png",
        thumb=output_path + "/conserv/{run}.{rep}/{run}.{rep}_conserv_thumb.png",
        r=output_path + "/conserv/{run}.{rep}/{run}.{rep}_conserv.R",
        score=output_path + "/conserv/{run}.{rep}/{run}.{rep}_conserv.txt",
    params:
        db=config['conservation'],
        # script="conservation_plot.py" if os.path.isdir(config['conservation']) else "conservation_onebw_plot.py",
        width=4000,
        #run = lambda wildcards: wildcards.run,
        run="{run}.{rep}" ,
        main_output_path=output_path
    message: "CONSERVATION: calling conservation script"
    log: output_path + "/logs/conservation/{run}.{rep}.log"
    benchmark: output_path + "/Benchmark/{run}.{rep}_conservation_plotConservation.benchmark"
    conda: "../envs/conservation/conservation.yaml"
    shell:
        "cidc_chips/modules/scripts/conservation_onebw_plot.py -t Conservation_at_summits -d {params.db} -o  {params.main_output_path}/conserv/{params.run}/{params.run}_conserv -l Peak_summits {input} -w {params.width} > {output.score} 2>>{log}"
