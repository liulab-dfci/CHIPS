#MODULE: targets module--use BETA to calculate the Regular Potential

#PARAMETERS
# _logfile=output_path + "/logs/targets.log"

target_decay_rate = 10000

def targets_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config["runs"].keys():
        #NOTE: using the fact that _reps has this info parsed already!
        for rep in _reps[run]:
            #GENERATE Run name: concat the run and rep name
            runRep = "%s.%s" % (run, rep)
            ls.append(output_path + "/targets/%s/%s_gene_score_5fold.txt" % (runRep,runRep))
            ls.append(output_path + "/targets/%s/%s_gene_score.txt" % (runRep,runRep))
            ls.append(output_path + "/targets/%s/%s_gene_score_1k.txt" % (runRep,runRep))
            ls.append(output_path + "/targets/%s/%s_gene_score_100k.txt" % (runRep,runRep))
    return ls

rule targets_all:
    input:
        targets_targets

rule targets_get5FoldPeaksRPScore:
    input:
        output_path + "/peaks/{run}.{rep}/{run}.{rep}_5fold_peaks.bed"
    output:
        output_path + "/targets/{run}.{rep}/{run}.{rep}_gene_score_5fold.txt"
    params:
        genome=config['geneBed'],
        decay=target_decay_rate,
        scripts="targets_getRP.py"
    message: "REGULATORY: get RP score of 5 fold peaks"
    log:output_path + "/logs/targets/{run}.{rep}.log"
    benchmark: output_path + "/Benchmark/{run}.{rep}_targets_get5FoldPeaksRPScore.benchmark"
    shell:
        "python cidc_chips/modules/scripts/targets_RegPotential_Version2.py -p {input} -a {params.genome} -n {output} -d {params.decay}"

rule targets_getTopPeaks:
    input:
        output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.bed"
    output:
        output_path + "/targets/{run}.{rep}/{run}.{rep}_peaks_top_reg.bed"
    params:
        peaks = 10000
    log:output_path + "/logs/targets/{run}.{rep}.log"
    benchmark: output_path + "/Benchmark/{run}.{rep}_targets_getTopPeaks.benchmark"
    message: "REGULATORY: get top summits for regpotential"
    shell:
        "head -n {params.peaks} {input} > {output}"

rule targets_getTopPeaksRPScore:
    input:
        output_path + "/targets/{run}.{rep}/{run}.{rep}_peaks_top_reg.bed"
    output:
        output_path + "/targets/{run}.{rep}/{run}.{rep}_top10k_gene_score.txt"
    params:
        genome=config['geneBed'],
        decay=target_decay_rate
    message: "REGULATORY: get RP score of top peaks"
    log:output_path + "/logs/targets/{run}.{rep}.log"
    benchmark: output_path + "/Benchmark/{run}.{rep}_targets_getTopPeaksRPScore.benchmark"
    shell:
        "python cidc_chips/modules/scripts/targets_RegPotential_Version2.py -p {input} -a {params.genome} -n {output} -d {params.decay}"


rule targets_getAllPeaksRPScore:
    input:
        output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.bed"
    output:
        output_path + "/targets/{run}.{rep}/{run}.{rep}_gene_score.txt"
    params:
        genome=config['geneBed'],
        decay=target_decay_rate
    message: "REGULATORY: get RP score of all peaks with 10k decay rate"
    log:output_path + "/logs/targets/{run}.{rep}.log"
    benchmark: output_path + "/Benchmark/{run}.{rep}_targets_getAllPeaksRPScore.benchmark"
    shell:
        "python cidc_chips/modules/scripts/targets_RegPotential_Version2.py -p {input} -a {params.genome} -n {output} -d {params.decay}"

rule targets_getAllPeaksRPScore1k:
    input:
        output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.bed"
    output:
        output_path + "/targets/{run}.{rep}/{run}.{rep}_gene_score_1k.txt"
    params:
        genome=config['geneBed'],
        decay=1000
    message: "REGULATORY: get RP score of all peaks with 1k decay rate"
    log:output_path + "/logs/targets/{run}.{rep}.log"
    benchmark: output_path + "/Benchmark/{run}.{rep}_targets_getAllPeaksRPScore1k.benchmark"
    shell:
        "python cidc_chips/modules/scripts/targets_RegPotential_Version2.py -p {input} -a {params.genome} -n {output} -d {params.decay}"

rule targets_getAllPeaksRPScore100k:
    input:
        output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.bed"
    output:
        output_path + "/targets/{run}.{rep}/{run}.{rep}_gene_score_100k.txt"
    params:
        genome=config['geneBed'],
        decay=100000
    message: "REGULATORY: get RP score of all peaks with 100k decay rate"
    log:output_path + "/logs/targets/{run}.{rep}.log"
    benchmark: output_path + "/Benchmark/{run}.{rep}_targets_getAllPeaksRPScore100k.benchmark"
    shell:
        "python cidc_chips/modules/scripts/targets_RegPotential_Version2.py -p {input} -a {params.genome} -n {output} -d {params.decay}"
