#MODULE: epicypher- quantify epicypher spike-ins
# _logfile=output_path + "/logs/epicypher.log"
_threads=8

def epicypher_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    if 'epicypher_analysis' in config and config['epicypher_analysis'] == 'SNAP':
        for sample in config["samples"]:
            ls.append(output_path + "/epicypher.snap/%s/%s.epicypher.quant.txt" % (sample,sample))
    elif 'epicypher_analysis' in config and config['epicypher_analysis'] == 'CAP':
        for sample in config["samples"]:
            ls.append(output_path + "/epicypher.cap/%s/%s.epicypher.quant.txt" % (sample,sample))
    elif 'epicypher_analysis' in config and config['epicypher_analysis'] == 'KACYL':
        for sample in config["samples"]:
            ls.append(output_path + "/epicypher.kacyl/%s/%s.epicypher.quant.txt" % (sample,sample))
    else:
        print("ERROR: epicypher_analysis not specified to either 'SNAP', 'CAP', or 'KACYL'")
        sys.exit()
    return ls

#OBSOLETE
# def getUnmappedReads(wildcards):
#     sample = wildcards.sample
#     ls = [output_path + "/align/%s/%s.unmapped.fq.gz" % (sample, sample)]
#     if len(config["samples"][wildcards.sample]) == 2:
#         ls.append(output_path + "/align/%s/%s.unmapped.fq2.gz" % (sample, sample))
#     return ls

rule epicypher_all:
    input:
        epicypher_targets

rule epicypher_alignToEpicypher:
    """Align unmapped reads to epicypher assembly"""
    input:
        output_path + "/align/{sample}/{sample}.unmapped.fq.gz"
    params:
        epicypher_index=lambda wildcards: "cidc_chips/static/epicypher/%s/epicypher.fa" % wildcards.ttype,
        #check for PE mate
        mate2 =  lambda wildcards: output_path + "/align/{sample}/{sample}.unmapped.fq2.gz" if len(config["samples"][wildcards.sample]) == 2 else ""
    output:
        output_path + "/epicypher.{ttype}/{sample}/{sample}.epicypher.bam"
    threads: _threads
    message: "EPICYPHER: aligning unmapped reads to epicypher assembly"
    # log: output_path + "/logs/epicypher/{sample}.log"
    conda: "../envs/epicypher/epicypherl.yaml"
    shell:
        "bwa mem -t {threads} {params.epicypher_index} {input} {params.mate2} | samtools view -Sb - > {output}"

rule epicypher_sortEpicypher:
    input:
        output_path + "/epicypher.{ttype}/{sample}/{sample}.epicypher.bam"
    output:
        output_path + "/epicypher.{ttype}/{sample}/{sample}.epicypher.sorted.bam"
    message: "EPICYPHER: sorting epicypher bam file"
    # log: output_path + "/logs/epicypher/{sample}.log"
    threads: _threads
    conda: "../envs/epicypher/epicypherl.yaml"
    shell:
        "sambamba sort {input} -o {output} -t {threads}"# 2>>{log}"

rule epicypher_uniquelyMapped:
    """Get uniquely mapped reads from epicypher.bam"""
    input:
        output_path + "/epicypher.{ttype}/{sample}/{sample}.epicypher.sorted.bam"
    output:
        output_path + "/epicypher.{ttype}/{sample}/{sample}.epicypher.sorted.unique.bam"
    message: "EPICYPHER: get uniquely mapped reads"
    threads: _threads
    # log: output_path + "/logs/epicypher/{sample}.log"
    conda: "../envs/epicypher/epicypherl.yaml"
    shell:
        "samtools view -bq 1 -@ {threads} {input} > {output}"

rule epicypher_indexEpicypher:
    input:
        output_path + "/epicypher.{ttype}/{sample}/{sample}.epicypher.sorted.unique.bam"
    output:
        output_path + "/epicypher.{ttype}/{sample}/{sample}.epicypher.sorted.unique.bam.bai"
    message: "EPICYPHER: indexing epicypher bam file"
    # log:output_path + "/logs/epicypher/{sample}.log"
    conda: "../envs/epicypher/epicypherl.yaml"
    shell:
        "sambamba index {input} {output}"

rule epicypher_quantify:
    """Quantify the reads for each spike-in mark"""
    input:
        bam=output_path + "/epicypher.{ttype}/{sample}/{sample}.epicypher.sorted.unique.bam",
        bai=output_path + "/epicypher.{ttype}/{sample}/{sample}.epicypher.sorted.unique.bam.bai"
    params:
        script = lambda wildcards: "cidc_chips/modules/scripts/epicypher_%s_quant.py" % wildcards.ttype
    output:
        output_path + "/epicypher.{ttype}/{sample}/{sample}.epicypher.quant.txt"
    message: "EPICYPHER: quantifying epicypher marks"
    threads: _threads
    # log: output_path + "/logs/epicypher/{sample}.log"
    conda: "../envs/epicypher/epicypherl.yaml"
    shell:
        "{params.script} -b {input.bam} > {output}"
