#MODULE: Align fastq files to genome - BWA specific calls
#PARAMETERS:
# _logfile=output_path + "/logs/align.log"
import subprocess

_bwa_threads=8
_bwa_q="5"
_bwa_l="32"
_bwa_k="2"


def getFastq(wildcards):
    return config["samples"][wildcards.sample]


def getAlnFastq(wildcards):
    return config["samples"][wildcards.sample][int(wildcards.mate)]


def getMates(wildcards):
    s = wildcards.sample
    files = config["samples"][s]
    return ["%s/align/%s/%s_%s_aln.sai" % (output_path,s,s,m) for m in range(len(files))]


def getRunType(wildcards):
    s = wildcards.sample
    return "sampe" if len(config["samples"][s]) > 1 else "samse"


checkpoint align_readsLength:
    input:
        output_path + "/fastqc/{sample}/{sample}_100k_fastqc/fastqc_data.txt"
    output:
        # direc=directory(output_path + "/fastqc/{sample}/"),
        length=temp(output_path + "/fastqc/{sample}/{sample}_read_length.txt")
    shell:
        # "mkdiroutput_path +  /fastqc/{wildcards.sample}/;"
        "python cidc_chips/modules/scripts/align_getReadLength.py -f {input} -o {output.length}"


rule align_bwaMem:
    input:
        getFastq
    output:
        temp(output_path + "/align/{sample}/{sample}_mem.bam")
    params:
        index=config['bwa_index'],
        sentieon=config.get("sentieon", ""),
        read_group=lambda wildcards: "@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA" % (wildcards.sample, wildcards.sample)
    threads: _bwa_threads
    message: "ALIGN: Running BWA mem for alignment for {input}"
    log: output_path + "/logs/align/{sample}.log"
    benchmark: output_path + "/Benchmark/{sample}_align_bwaMem.benchmark"
    conda: "../envs/align/align_bwa.yaml"
    shell:
        "{params.sentieon} bwa mem -t {threads} -R \"{params.read_group}\" {params.index} {input} | samtools view -Sb - > {output} 2>>{log}"


rule align_bwaAln:
    input:
        getAlnFastq
    output:
        sai=temp(output_path + "/align/{sample}/{sample}_{mate}_aln.sai")
    params:
        index=config['bwa_index'],
        bwa_q = _bwa_q,
        bwa_l = _bwa_l,
        bwa_k = _bwa_k,
        sentieon=config["sentieon"] if ("sentieon" in config) and config["sentieon"] else ""
    threads: _bwa_threads
    message: "ALIGN: Running BWA aln for alignment for {input}"
    # log: output_path + "/logs/align/{sample}.log"
    benchmark: output_path + "/Benchmark/{sample}_{mate}_align_bwaAln.benchmark"
    shell:
        "{params.sentieon} bwa aln -q {params.bwa_q} -l {params.bwa_l} -k {params.bwa_k} -t {threads} {params.index} {input} > {output.sai}"

rule align_bwaConvert:
    input:
        sai=getMates,
        fastq=getFastq
    output:
        temp(output_path + "/align/{sample}/{sample}_aln.bam")
    params:
        run_type= getRunType,
        index=config['bwa_index'],
        #NOTE: this is a hack b/c snakemake didn't like the - in the shell cmd
        hack="view -bS -",
        sentieon=config["sentieon"] if ("sentieon" in config) and config["sentieon"] else "",
        read_group=lambda wildcards: "@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA" % (wildcards.sample, wildcards.sample)
    threads: _bwa_threads
    message: "ALIGN: Converting BWA alignment to BAM"
    # log: output_path + "/logs/align/{sample}.log"
    benchmark: output_path + "/Benchmark/{sample}_align_bwaConvert.benchmark"
    shell:
        """{params.sentieon} bwa {params.run_type} -r \"{params.read_group}\" {params.index} {input.sai} {input.fastq} | samtools {params.hack} > {output}"""


def aggregate_align_input(wildcards):
    # decision based on content of output file
    with open(checkpoints.align_readsLength.get(sample=wildcards.sample).output[0]) as f:
        if int(f.read().strip()) >= 40:
            return output_path + "/align/{sample}/{sample}_mem.bam"
        else:
            return output_path + "/align/{sample}/{sample}_aln.bam"


checkpoint align_aggregate:
    input:
        aggregate_align_input
    output:
        temp(output_path + "/align/{sample}/{sample}.bam")
    shell:
        "mv {input} {output}"


rule align_macsRunInfo:
    """Dump the current version of bwa into a text file for the report"""
    output:
        output_path + "/align/run_info.txt"
    message: "ALIGN/REPORT - collection bwa version info"
    conda: "../envs/align/align_bwa.yaml"
    shell:
        "cidc_chips/modules/scripts/align_parseBwaVersion.py -o {output}"
