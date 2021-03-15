#MODULE: contamination module to check for sample contamination
# _logfile=output_path + "/logs/contamination.log"
_bwa_q="5"
_bwa_l="32"
_bwa_k="2"
_bwa_threads=8

###############################################################################
# HELPERS
###############################################################################
_contaminationPanel= config.get('contamination_panel', [])

def extractIndexName(path):
    """Given a contamination panel path, e.g. /some/path/to/BWA/hg19.fa,
    this fn returns the index filename, e.g. hg19"""
    tmp = path.split("/")[-1]
    #lop off the suffix
    return ".".join(tmp.split(".")[:-1])
_contaminationNames = list(map(extractIndexName, _contaminationPanel))
_contaminationDict = dict(zip(_contaminationNames, _contaminationPanel))
###############################################################################

def contamination_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        for panel in _contaminationNames:
            ls.append(output_path + "/contam/%s/%s.%s.txt" % (sample, sample, panel))
            ls.append(output_path + "/contam/%s/%s_contamination.txt" % (sample, sample))
    ls.append(output_path + "/contam/contamination.csv")
    return ls

def get100k_sample(wildcards):
    """NOTE: there are two different 100k files, either _100k.fastq or
    _100k.bam.fastq.  The former comes from runs that have fastqs as input,
    while latter started w/ bams as inputs.

    This fn distinguishes which. (based on fastqc.snakefile getFastqcInput
    """
    s = wildcards.sample
    first_file = config["samples"][wildcards.sample][0]
    if first_file.endswith('.bam'):
        ret = [output_path + "/align/%s/%s_100k.bam.fastq" % (s,s)]
    else:
        ret = [output_path + "/align/%s/%s_100k.fastq" % (s,s)]
    return ret

rule contamination_all:
    input:
        contamination_targets

rule contamination_alignContamination:
    """For each sample, run an alignment for each entry in the
    contaminationPanel and get the mapping rate"""
    input:
        get100k_sample
    output:
        temp(output_path + "/contam/{sample}/{sample}.{panel}.bam")
    params:
        index=lambda wildcards: _contaminationDict[wildcards.panel],
        bwa_q = _bwa_q,
        bwa_l = _bwa_l,
        bwa_k = _bwa_k,
        sample = lambda wildcards: wildcards.sample,
        panel = lambda wildcards: wildcards.panel
    threads: _bwa_threads
    message: "CONTAMINATION: checking {params.sample} against {params.panel}"
    log: output_path + "/logs/contamination/{sample}.{panel}.log"
    benchmark: output_path + "/Benchmark/{sample}.{panel}_contam.benchmark"
    conda: "../envs/contamination/contamination.yaml"
    shell:
        "bwa mem -t {threads} {params.index} {input} | samtools view -Sb - > {output} 2>>{log}"

rule contamination_Stats:
    """Extract the mapping stats for each species"""
    input:
        output_path + "/contam/{sample}/{sample}.{panel}.bam"
    output:
        #TEMP
        output_path + "/contam/{sample}/{sample}.{panel}.txt"
    params:
        sample = lambda wildcards: wildcards.sample,
        panel = lambda wildcards: wildcards.panel,
        #READ out the 5th row, the first element (and divide by 100/100000)
        awk_cmd = "\'BEGIN {RS=\'\\t\'}{print $23 / $1 * 100}\'"
    message: "CONTAMINATION: get mapping stats {params.sample}:{params.panel}"
    # log: output_path + "/logs/contamination/{sample}.log"
    conda: "../envs/contamination/contamination.yaml"
    shell:
        "samtools flagstat {input} | awk {params.awk_cmd} > {output} "#2>>{log}"

rule contamination_collectStats:
    """Collect the mapping stats across the entire panel"""
    input:
        expand(output_path + "/contam/{{sample}}/{{sample}}.{panel}.txt", panel=_contaminationNames)
    output:
        output_path + "/contam/{sample}/{sample}_contamination.txt"
    #conda: "../envs/contamination/contamination.yaml"
    run:
        for (n, f) in zip(_contaminationNames, input):
            shell("per=$(cat {f}) && echo {n} $per >> {output}")

rule contamination_collectAllContamination:
    """Aggregate all of the contamination stats into one large table/panel"""
    input:
        expand(output_path + "/contam/{sample}/{sample}_contamination.txt", sample=config['samples'])
    message: "Contamination: collecting contamination panel"
    # log: output_path + "/logs/contamination/{sample}.log"
    output:
        output_path + "/contam/contamination.csv"
    params:
        files = lambda wildcards, input: [" -f %s" % i for i in input]
    conda: "../envs/contamination/contamination.yaml"
    # run:
    #     files = " -f ".join(input)
    #     shell("cidc_chips/modules/scripts/contam_getStats.py -f {files} -o {output} 2>>{log}")
    shell:
        "cidc_chips/modules/scripts/contam_getStats.py {params.files} -o {output} "#2>>{log}"
