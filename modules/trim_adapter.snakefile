_fastp_threads=2

# def trim_targets(wildcards):
#     """Generates the targets for this module"""
#     ls = []
#     for sample in config["samples"]:
#       # fastp output json and html summarization
#       ls.append(output_path + "/trim_adaptor/%s/%s_fastp.json" % (sample,sample))
#       ls.append(output_path + "/trim_adaptor/%s/%s_fastp.html" % (sample,sample))
#       fastq_num = len(config["samples"][sample]) 
#       if fastq_num > 1: #paired-end
#           ls.append(expand(output_path + "/trim_adaptor/%s/%s_{mate}.trimmed.fq" % (sample,sample), mate = range(fastq_num)))
#       else: # single-end
#           ls.append(output_path + "/trim_adaptor/%s/%s_1.trimmed.fq" % (sample,sample))
#     return ls


# only make the {sample}_0.trimmed.fq as required target
# in the trim_fastp rule, we generate {sample}_0.trimmed.fq for single end data
# and {sample}_0.trimmed.fq {sample}_1.trimmed.fq for paired end data

def trim_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        # fastp output json and html summarization
        ls.append(output_path + "/trim_adaptor/%s/%s_fastp.json" % (sample,sample))
        ls.append(output_path + "/trim_adaptor/%s/%s_fastp.html" % (sample,sample))
        ls.append(output_path + "/trim_adaptor/%s/%s_0.trimmed.fq" % (sample,sample))
        ls.append(output_path + "/trim_adaptor/%s/%s_1.trimmed.fq" % (sample,sample))
    return ls

def getFastq(wildcards):
    return config["samples"][wildcards.sample]

## it is not possible to write an output function to determine the output file number
## based on the wildcards. use touch to trick snakemake instead.
## see https://bitbucket.org/snakemake/snakemake/issues/865/pre-determined-dynamic-output
## and https://stackoverflow.com/questions/56861913/split-bam-by-clusters-and-then-merge-bam-by-cluster-using-checkpoint

FASTP_VERSION = subprocess.check_output("fastp -v", shell=True)

# rule trim_fastp:
#   input: getFastq
#   output: 
#       fastqs = expand(output_path + "/trim_adaptor/{{sample}}/{{sample}}_{mate}.trimmed.fq", mate = range(2)),
#       json = output_path + "/trim_adaptor/{sample}/{sample}_fastp.json",
#       html = output_path + "/trim_adaptor/{sample}/{sample}_fastp.html" 
#   threads: _fastp_threads
#   message: "trimming adaptors for {input} using fastp"
#   log: output_path + "/logs/trim_adaptor/{sample}.log"
#   conda: "../envs/trimming/trim_fastp.yaml"
#   run:
#       if len(input) > 1:
#           shell("fastp --thread {threads} --detect_adapter_for_pe --in1 {input[0]} --in2 {input[1]} --out1 {output.fastqs[0]} --out2 {output.fastqs[1]} -h {output.html} -j {output.json} > {log} 2>&1 ")
#       else:
#           shell("fastp --thread {threads} --in1 {input[0]} --out1 {output.fastqs[0]} -h {output.html} -j {output.json} > {log} 2>&1")
#           shell("touch {output.fastqs[1]}")
    

rule trim_fastp:
    input: getFastq
    output: 
        fq1 = temp(output_path + "/trim_adaptor/{sample}/{sample}_0.trimmed.fq"),
        fq2 = temp(output_path + "/trim_adaptor/{sample}/{sample}_1.trimmed.fq"),
        json = output_path + "/trim_adaptor/{sample}/{sample}_fastp.json",
        html = output_path + "/trim_adaptor/{sample}/{sample}_fastp.html"
    threads: _fastp_threads
    message: "trimming adaptors for {input} using fastp"
    version: FASTP_VERSION
    log: output_path + "/logs/trim_adaptor/{sample}.log"
    conda: "../envs/trimming/trim_fastp.yaml"
    benchmark: output_path + "/Benchmark/{sample}_trim_adapter.benchmark"
    run:
        if len(input) > 1:
            shell("fastp --thread {threads} --detect_adapter_for_pe --in1 {input[0]} --in2 {input[1]} --out1 {output.fq1} --out2 {output.fq2} -h {output.html} -j {output.json} > {log} 2>&1")
        else:
            shell("fastp --thread {threads} --in1 {input[0]} --out1 {output.fq1} -h {output.html} -j {output.json} > {log} 2>&1")
            shell("touch{output.fq2}")




