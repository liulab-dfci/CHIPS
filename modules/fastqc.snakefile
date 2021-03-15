#MODULE: fastqc- sequence quality scores
#RULES:
#    sample_fastq: subsample 100k reads from each sample to perform fastqc analysis
# _logfile=output_path + "/logs/fastqc.log"

def fastqc_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        ls.append(output_path + "/fastqc/%s/%s_perSeqGC.txt" % (sample,sample))
        ls.append(output_path + "/fastqc/%s/%s_perSeqQual.txt" % (sample,sample))
        ls.append(output_path + "/fastqc/%s/%s_stats.csv" % (sample,sample))
        ls.append(output_path + "/fastqc/%s/%s_perSeqGC.png" % (sample,sample))
        ls.append(output_path + "/fastqc/%s/%s_perSeqGC_thumb.png" % (sample,sample))
        ls.append(output_path + "/fastqc/%s/%s_100k_fastqc/" % (sample, sample))
        ls.append(output_path + "/fastqc/%s/%s_100k_fastqc/fastqc_data.txt" % (sample, sample))
    ls.append(output_path + "/fastqc/fastqc.csv")
    return ls

def getFastqcInput(wildcards):
    """Get the input file for fastqc. It's either the 100k.fastq sample or the
     original bam"""
    s = wildcards.sample
    first_file = config["samples"][wildcards.sample][0]
    if first_file.endswith('.bam'):
        #CLEANER to check for .bam vs (.fastq, fq, fq.gz, etc)
        #ret = [first_file]
        #HACK to get fastqc to name things correctly.  IF we returned
        #the bam file, fastqc will name the output files accordingly
        ret = [output_path + "/align/%s/%s_100k.bam" % (s,s)]
    else:
        #HACK: need to give the EVALUATED cannonical path
        ret = [output_path + "/align/%s/%s_100k.fastq" % (s,s)]
    #print(ret)
    return ret

def getFastqcBam(wildcards):
    """Used in conjunction w/ getFastqcInput to just return the sample bam
    """
    s = wildcards.sample
    ret = [config["samples"][wildcards.sample][0]]
    #print(ret)
    return ret

def getFastq3(wildcards):
    """Get associated fastqs for each sample.
    NOTE: if PE, take the first pair"""
    s = config["samples"][wildcards.sample]
    return s[0]

rule fastqc_all:
    input:
        fastqc_targets

rule fastqc_downsampleFastq:
    """Subsample fastq"""
    input:
        getFastq3
    output:
        temp(output_path + "/align/{sample}/{sample}_100k.fastq")
    params:
        seed=11,
        #how many to sample
        size=100000
    message: "FASTQC: downsample {input} to 100k reads"
    log: output_path + "/logs/fastqc/{sample}.log"
    benchmark: output_path + "/Benchmark/{sample}_fastqc_downsample.benchmark"
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "seqtk sample -s {params.seed} {input} {params.size} > {output} 2>>{log}"

rule fastqc_sampleBam:
    """USED only when the sample-input is a bam file.
    Subsample bam to 100k reads"""
    input:
        getFastqcBam
    params:
        n=100000
    output:
        temp(output_path + "/align/{sample}/{sample}_100k.bam")
    message: "FASTQC: sampling 100k reads from bam for {input}"
    log: output_path + "/logs/fastqc/{sample}.log"
    benchmark: output_path + "/Benchmark/{sample}_fastqc_samplebam.benchmark"
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "cidc_chips/modules/scripts/fastqc_sampleBam.sh -i {input} -n {params.n} -o {output}"

rule fastqc_convertBamToFastq:
    """USED only when the sample-input is a bam file."""
    input:
        output_path + "/align/{sample}/{sample}_100k.bam"
    output:
        temp(output_path + "/align/{sample}/{sample}_100k.bam.fastq")
    message: "FASTQC: converting 100k.bam to 100k.fastq for {input}"
    log: output_path + "/logs/fastqc/{sample}.log"
    benchmark: output_path + "/Benchmark/{sample}_fastqc_bamtofastq.benchmark"
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "bamToFastq -i {input} -fq {output} 2>> {log}"

rule fastqc_callFastqc:
    """CALL FASTQC on each sub-sample"""
    input:
        getFastqcInput
    output:
        #MAKE temp
        directory(output_path + "/fastqc/{sample}/{sample}_100k_fastqc/"),
        output_path + "/fastqc/{sample}/{sample}_100k_fastqc/fastqc_data.txt",
        output_path + "/fastqc/{sample}/{sample}_100k_fastqc.html",
        output_path + "/fastqc/{sample}/{sample}_100k_fastqc.zip"
    #threads:
    params:
        sample = lambda wildcards: wildcards.sample,
        main_output_path = output_path
    message: "FASTQC: call fastqc for {input}"
    log: output_path + "/logs/fastqc/{sample}.log"
    benchmark: output_path + "/Benchmark/{sample}_fastqc_call.benchmark"
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "fastqc {input} --extract -o {params.main_output_path}/fastqc/{params.sample} 2>>{log}"

rule fastqc_getPerSequenceQual:
    """extract per sequence quality from fastqc_data.txt"""
    input:
        output_path + "/fastqc/{sample}/{sample}_100k_fastqc/fastqc_data.txt"
    output:
        output_path + "/fastqc/{sample}/{sample}_perSeqQual.txt"
    params:
        #DON'T forget quotes
        section="'Per sequence quality'"
    message: "FASTQC: get_PerSequenceQual for {input}"
    log: output_path + "/logs/fastqc/{sample}.log"
    benchmark: output_path + "/Benchmark/{sample}_fastqc_SeqQual.benchmark"
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "cidc_chips/modules/scripts/fastqc_dataExtract.py -f {input} -s {params.section} > {output} 2>>{log}"

rule fastqc_getPerSequenceGC:
    """extract per sequence GC contentfrom fastqc_data.txt"""
    input:
        output_path + "/fastqc/{sample}/{sample}_100k_fastqc/fastqc_data.txt"
    output:
        output_path + "/fastqc/{sample}/{sample}_perSeqGC.txt"
    params:
        #DON'T forget quotes
        section="'Per sequence GC content'"
    message: "FASTQC: get_PerSequenceGC for {input}"
    log: output_path + "/logs/fastqc/{sample}.log"
    benchmark: output_path + "/Benchmark/{sample}_fastqc_GC.benchmark"
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "cidc_chips/modules/scripts/fastqc_dataExtract.py -f {input} -s {params.section} > {output} 2>>{log}"

rule fastqc_extractFastQCStats:
    """extract per sequence GC content, and seq qual stats from fastqc run"""
    input:
        gc = output_path + "/fastqc/{sample}/{sample}_perSeqGC.txt",
        qual = output_path + "/fastqc/{sample}/{sample}_perSeqQual.txt"
    output:
        output_path + "/fastqc/{sample}/{sample}_stats.csv"
    message: "FASTQC: extract_FastQCStats from {input}"
    log:output_path + "/logs/fastqc/{sample}.log"
    benchmark: output_path + "/Benchmark/{sample}_fastqc_extractFastQCStats.benchmark"
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "cidc_chips/modules/scripts/fastqc_stats.py -a {input.qual} -b {input.gc} > {output} 2>>{log}"

rule fastqc_collectFastQCStats:
    """Collect and parse out the fastqc stats for the ALL of the samples"""
    input:
        expand(output_path + "/fastqc/{sample}/{sample}_stats.csv", sample=sorted(config["samples"]))
    output:
        output_path + "/fastqc/fastqc.csv"
    message: "FASTQC: collect and parse ALL mapping stats for {input}"
    # log: output_path + "/logs/fastqc/{sample}.log"
    #conda: "../envs/fastqc/fastqc.yaml"
    run:
        files = " -f ".join(input)
        shell("cidc_chips/modules/scripts/fastqc_getFastQCStats.py -f {files} > {output}")

rule fastqc_plotFastQCGC:
    """Plots the GC distribution of the sample according to data in
    perSeqGC.txt.
    Generates a full-size image and *thumbnail image* (embedded into report)
    """
    input:
        gc = output_path + "/fastqc/{sample}/{sample}_perSeqGC.txt",
    output:
        png=output_path + "/fastqc/{sample}/{sample}_perSeqGC.png",
        thumb=output_path + "/fastqc/{sample}/{sample}_perSeqGC_thumb.png",
    message:
        "FASTQC: generating GC content distrib. plots for {input}"
    log: output_path + "/logs/fastqc/{sample}.log"
    conda: "../envs/fastqc/fastqc.yaml"
    shell:
        "Rscript cidc_chips/modules/scripts/fastqc_plotGC.R {input.gc} {output.png} {output.thumb}"
