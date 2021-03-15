#MODULE: Align fastq files to genome - common rules
#import os
_align_threads=8
_sambamba_sort_mem=4

def align_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        ls.append(output_path + "/align/%s/%s.sorted.bam" % (sample,sample))
        ls.append(output_path + "/align/%s/%s.sorted.bam.bai" % (sample,sample))
        ls.append(output_path + "/align/%s/%s_unique.sorted.bam" % (sample,sample))
        ls.append(output_path + "/align/%s/%s_unique.sorted.bam.bai"%(sample,sample))
        # ls.append(output_path + "/align/%s/%s_unique.sorted.dedup.bam" % (sample,sample))
        # ls.append(output_path + "/align/%s/%s_unique.sorted.dedup.bam.bai" % (sample,sample))
        ls.append(output_path + "/align/%s/%s_mapping.png" % (sample,sample))
        if len(config["samples"][sample]) > 1 and ('cutoff' in config) and config['cutoff']:
            ls.append(output_path + "/align/%s/%s_unique.sorted.sub%s.bam"%(sample,sample,config['cutoff']))
        ls.append(output_path + "/align/%s/%s.unmapped.fq.gz" % (sample,sample))
        ls.append(output_path + "/align/%s/%s_readsPerChrom.txt" % (sample,sample))
    ls.append(output_path + "/align/mapping.csv")
    ls.append(output_path + "/align/run_info.txt")
    return ls


def getBam(wildcards):
    """This input fn will check to see if the user specified a .fastq or a .bam
    for the sample.  IF the former (.fastq), will simply return the canonical
    path, otherwise (.bam) will return the user-specified (bam) path"""
    #CHECK first entry's file suffix
    s = wildcards.sample
    first_file = config["samples"][wildcards.sample][0]
    ret = output_path + "/align/%s/%s.bam" % (s,s)
    if first_file.endswith('.bam'):
        #CLEANER to check for .bam vs (.fastq, fq, fq.gz, etc)
        ret = first_file
    return [ret]

## this will require sambamba to use 2x size of the bam file for sorting, which is too big and may cause memory issues.
## define the memory _sambamba_sort_mem on top of this snakefile instead
def getSortMemory(wildcards):
    bam_size = math.ceil(os.path.getsize(checkpoints.align_aggregate.get(sample=wildcards.sample).output[0])/1024/1024/1024)
    memory = bam_size*2
    return str(memory)

def getUniqueSortMemory(wildcards):
    bam_size = math.ceil(os.path.getsize(checkpoints.align_uniquelyMappedReads.get(sample=wildcards.sample).output[0])/1024/1024/1024)
    memory = bam_size*2
    return str(memory)

rule align_all:
    input:
        align_targets

checkpoint align_uniquelyMappedReads:
    """Get the uniquely mapped reads"""
    input:
        #output_path + "/align/{sample}/{sample}.bam"
        # getBam
        output_path + "/align/{sample}/{sample}.sorted.bam"
    output:
        temp(output_path + "/align/{sample}/{sample}_unique.bam")
    message: "ALIGN: Filtering for uniquely mapped reads for {input}"
    log: output_path + "/logs/align/{sample}.log"
    threads: _align_threads
    conda: "../envs/align/align_common.yaml"
    shell:
        #NOTE: this is the generally accepted way of doing this as multiply
        #mapped reads have a Quality score of 0
        #NOTE: -@ = --threads
        "samtools view -bq 1 -@ {threads} {input} > {output}"

rule align_mapStats:
    """Get the mapping stats for each aligment run"""
    input:
        #bam=output_path + "/align/{sample}/{sample}.bam",
        # bam=getBam,
        bam=output_path + "/align/{sample}/{sample}.sorted.bam",
        uniq_bam=output_path + "/align/{sample}/{sample}_unique.bam"
    output:
        #temp(output_path + "/align/{sample}/{sample}_mapping.txt")
        output_path + "/align/{sample}/{sample}_mapping.txt"
    threads: _align_threads
    message: "ALIGN: get mapping stats for {input}"
    log: output_path + "/logs/align/{sample}.log"
    conda: "../envs/align/align_common.yaml"
    #CAN/should samtools view be multi-threaded--
    #UPDATE: tricky on how to do this right w/ compounded commands
    shell:
        #FLAGSTATS is the top of the file, and we append the uniquely mapped
        #reads to the end of the file
        "samtools flagstat {input.bam} > {output} 2>>{log}"
        " && samtools view -c {input.uniq_bam} >> {output} 2>> {log}"

rule align_mapFigure:
    input:
        output_path + "/align/{sample}/{sample}_mapping.txt"
    output:
        output_path + "/align/{sample}/{sample}_mapping.png"
    message: "ALIGN: plot mapping rate for {input}"
    conda: "../envs/align/align_common.yaml"
    shell:
        "cidc_chips/modules/scripts/align_mappedFigure.py -f {input} -o {output}"

rule align_collectMapStats:
    """Collect and parse out the mapping stats for the ALL of the samples"""
    input:
        #samples sorted to match order of rest of report
        expand(output_path + "/align/{sample}/{sample}_mapping.txt", sample=sorted(config["samples"]))
    output:
        output_path + "/align/mapping.csv"
    params:
        files = lambda wildcards, input: [" -f %s" % i for i in input]
    message: "ALIGN: collect and parse ALL mapping stats for {input}"
    # log: output_path + "/logs/align/{sample}.log"
    conda: "../envs/align/align_common.yaml"
    #NOTE: can't do conda envs with run
    #run:
    #    files = " -f ".join(input)
    #    shell("cidc_chips/modules/scripts/align_getMapStats.py -f {files} > {output} 2>>{log}")
    shell:
        "cidc_chips/modules/scripts/align_getMapStats.py {params.files} > {output}"


if config.get("sentieon"):
    rule align_sortBamsBySentieon:
        """General sort rule--take a bam {filename}.bam and
        output {filename}.sorted.bam"""
        input:
            #output_path + "/align/{sample}/{filename}.bam"
            getBam
        output:
            output_path + "/align/{sample}/{sample}.sorted.bam",
            #output_path + "/align/{sample}/{sample}.sorted.bam.bai"
        message: "ALIGN: sort bam file for {input}"
        log: output_path + "/logs/align/{sample}.log"
        conda: "../envs/align/align_common.yaml"
        params:
            sentieon=config['sentieon'],
        threads: _align_threads
        shell:
            "{params.sentieon} util sort -o {output} -i {input} 2>>{log}"

    rule align_sortUniqueBamsBySentieon:
        """General sort rule--take a bam {filename}.bam and
        output {filename}.sorted.bam"""
        input:
            output_path + "/align/{sample}/{sample}_unique.bam"
        output:
            #CANNOT temp this b/c it's used by qdnaseq!
            output_path + "/align/{sample}/{sample}_unique.sorted.bam",
            #output_path + "/align/{sample}/{sample}_unique.sorted.bam.bai"
        message: "ALIGN: sort unique bam file for {input}"
        log: output_path + "/logs/align/{sample}.log"
        conda: "../envs/align/align_common.yaml"
        params:
            sentieon=config['sentieon'],
        threads: _align_threads
        shell:
            "{params.sentieon} util sort -o {output} -i {input} 2>>{log}"

    rule align_scoreSampleBySentieon:
        "Calls sentieon driver  --fun score_info on the sample"
        input:
            bam=output_path + "/align/{sample}/{sample}_unique.sorted.bam",
            bai=output_path + "/align/{sample}/{sample}_unique.sorted.bam.bai",
        output:
            score=temp(output_path + "/align/{sample}/{sample}_unique.sorted.score.txt"),
            idx=temp(output_path + "/align/{sample}/{sample}_unique.sorted.score.txt.idx"),
        message: "ALIGN: score sample for {input.bam}"
        log: output_path + "/logs/align/{sample}/align.scoreSample.{sample}.log"
        threads: _align_threads
        params:
            sentieon=config['sentieon'],
        # group: "align"
        # benchmark:
        #     "benchmarks/align/{sample}/{sample}.scoreSample.txt"
        shell:
            """{params.sentieon} driver -t {threads} -i {input.bam} --algo LocusCollector --fun score_info {output.score} """

    rule align_dedupSortedUniqueBamBySentieon:
        """Dedup sorted unique bams using sentieon
         output {sample}_unique.sorted.dedup.bam"""
        input:
            bam=output_path + "/align/{sample}/{sample}_unique.sorted.bam",
            bai=output_path + "/align/{sample}/{sample}_unique.sorted.bam.bai",
            score=output_path + "/align/{sample}/{sample}_unique.sorted.score.txt",
            idx=output_path + "/align/{sample}/{sample}_unique.sorted.score.txt.idx"
        output:
            bamm=output_path + "/align/{sample}/{sample}_unique.sorted.dedup.bam",
            met=temp(output_path + "/align/{sample}/{sample}_unique.sorted.dedup.metric.txt"),
        message: "ALIGN: dedup sorted unique bam file for {input.bam}"
        log: output_path + "/logs/align/{sample}/align.dedupSortedUniqueBam.{sample}.log"
        threads: _align_threads
        params:
            sentieon=config['sentieon'],
        # group: "align"
        # benchmark:
        #     "benchmarks/align/{sample}/{sample}.dedupSortedUniqueBam.txt"
        shell:
            "{params.sentieon} driver -t {threads} -i {input.bam} --algo Dedup --rmdup --score_info {input.score} --metrics {output.met} {output.bamm}"

else:
    rule align_dedupSortedUniqueBams:
        """Dedup sorted unique bams using PICARD
        output {sample}_unique.sorted.dedup.bam"""
        input:
            output_path + "/align/{sample}/{sample}_unique.sorted.bam"
        output:
            bam = output_path + "/align/{sample}/{sample}_unique.sorted.dedup.bam",
            # bai = output_path + "/align/{sample}/{sample}_unique.sorted.dedup.bam.bai"
        message: "ALIGN: dedup sorted unique bam file for {input}"
        log: output_path + "/logs/align/{sample}.log"
        benchmark: output_path + "/Benchmark/{sample}_align_dedup.benchmark"
        conda: "../envs/align/align_common.yaml"
        threads: _align_threads
        shell:
            "picard MarkDuplicates I={input} O={output} REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE={log} 2>> {log}"
            # "sambamba markdup -t {threads} -r --overflow-list-size 600000 --hash-table-size 800000 --sort-buffer-size 4096 {input} {output.bam}"

    rule align_sortBams:
        """General sort rule--take a bam {filename}.bam and
        output {filename}.sorted.bam"""
        input:
            #output_path + "/align/{sample}/{filename}.bam"
            getBam
        output:
            output_path + "/align/{sample}/{sample}.sorted.bam",
            #output_path + "/align/{sample}/{sample}.sorted.bam.bai"
        message: "ALIGN: sort bam file for {input}"
        log: output_path + "/logs/align/{sample}.log"
        benchmark: output_path + "/Benchmark/{sample}_align_sortBam.benchmark"
        conda: "../envs/align/align_common.yaml"
        params:
            memory = _sambamba_sort_mem
        threads: _align_threads
        shell:
            "sambamba sort {input} -o {output} -t {threads} -m {params.memory}G 2>>{log}"

    rule align_sortUniqueBams:
        """General sort rule--take a bam {filename}.bam and
        output {filename}.sorted.bam"""
        input:
            output_path + "/align/{sample}/{sample}_unique.bam"
        output:
            #CANNOT temp this b/c it's used by qdnaseq!
            output_path + "/align/{sample}/{sample}_unique.sorted.bam",
            #output_path + "/align/{sample}/{sample}_unique.sorted.bam.bai"
        message: "ALIGN: sort bam file for {input}"
        log: output_path + "/logs/align/{sample}.log"
        benchmark: output_path + "/Benchmark/{sample}_align_sortUnique.benchmark"
        conda: "../envs/align/align_common.yaml"
        params:
            memory = _sambamba_sort_mem
        threads: _align_threads
        shell:
            "sambamba sort {input} -o {output} -t {threads} -m {params.memory}G 2>>{log}"

rule align_filterBams:
    """Filter out the long reads to get more accurate results in peaks calling"""
    input:
        output_path + "/align/{sample}/{sample}_unique.sorted.bam"
    output:
        output_path + "/align/{sample}/{sample}_unique.sorted.sub%s.bam" % str(config['cutoff'])
    message: "ALIGN: filter bam files for {input}"
    log: output_path + "/logs/align/{sample}.log"
    params:
        cutoff = config['cutoff']
    conda: "../envs/align/align_common.yaml"
    shell:
        "samtools view -h {input} | awk '($9 <= {params.cutoff} && $9 >= (-1)*{params.cutoff}) || $1 ~ /^@/' | samtools view -bS - > {output}"

rule align_indexBam:
    """Index bam file"""
    input:
        output_path + "/align/{sample}/{prefix}.bam"
    output:
        output_path + "/align/{sample}/{prefix}.bam.bai"
    message: "ALIGN: indexing bam file {input}"
    # log: output_path + "/logs/align/{sample}.log"
    conda: "../envs/align/align_common.yaml"
    shell:
        "sambamba index {input} {output}"

rule align_extractUnmapped:
    """Extract the unmapped reads and save as {sample}.unmapped.bam"""
    input:
        #output_path + "/align/{sample}/{sample}.bam"
        # getBam
        output_path + "/align/{sample}/{sample}.sorted.bam"
    output:
        temp(output_path + "/align/{sample}/{sample}.unmapped.bam")
    message: "ALIGN: extract unmapped reads for {input}"
    log: output_path + "/logs/align/{sample}.log"
    conda: "../envs/align/align_common.yaml"
    threads: _align_threads
    shell:
        #THIS extracts all unmapped reads
        #"samtools view -b -f 4 --threads {threads} {input} >{output} 2>>{log}"
        #THIS extracts all READ (pairs) where at least one in unmapped
        #ref: https://www.biostars.org/p/56246/ search "rgiannico"
        #NOTE: -@ = --threads
        "samtools view -b -F 2 -@ {threads} {input} > {output} 2>>{log}"

rule align_bamToFastq:
    """Convert the unmapped.bam to fastq"""
    input:
        output_path + "/align/{sample}/{sample}.unmapped.bam"
    output:
        output_path + "/align/{sample}/{sample}.unmapped.fq"
    params:
        #handle PE alignments!
        mate2 = lambda wildcards: "-fq2 " + output_path + "/align/%s/%s.unmapped.fq2" % (wildcards.sample,wildcards.sample) if len(config["samples"][wildcards.sample]) == 2 else ""
    message: "ALIGN: convert unmapped bam to fastq for {input}"
    log: output_path + "/logs/align/{sample}.log"
    conda: "../envs/align/align_common.yaml"
    shell:
        "bamToFastq -i {input} -fq {output} {params.mate2}"

rule align_gzipUnmappedFq:
    """gzip unmapped fq(s)"""
    input:
        output_path + "/align/{sample}/{sample}.unmapped.fq"
    output:
        output_path + "/align/{sample}/{sample}.unmapped.fq.gz"
    params:
        #handle PE alignments!
        mate2 = lambda wildcards: output_path + "/align/%s/%s.unmapped.fq2" % (wildcards.sample,wildcards.sample) if len(config["samples"][wildcards.sample]) == 2 else ""
    message: "ALIGN: gzip unmapped fq files for {input}"
    log: output_path + "/logs/align/{sample}.log"
    conda: "../envs/align/align_common.yaml"
    shell:
        "gzip -f {input} {params} 2>>{log}"

rule align_readsPerChromStat:
    """For each sample, generates a _readsPerChrom.txt file, which is:
    chr1   #readsOnChr1
    ...
    chrX   #readsOnChrX
    """
    input:
        bam = output_path + "/align/{sample}/{sample}.sorted.bam",
        #NOTE: even though we don't use the bai, we need to ensure bam sorted
        bai = output_path + "/align/{sample}/{sample}.sorted.bam.bai"
    params:
        awk_call = """awk '{print $1 \"\\t\" $3; s+=$3} END {print \"total reads = \" s}'"""
    output:
        output_path + "/align/{sample}/{sample}_readsPerChrom.txt"
    message: "ALIGN: collecting the number of reads per chrom for {input.bam}"
    log: output_path + "/logs/align/{sample}.log"
    conda: "../envs/align/align_common.yaml"
    shell:
        "cidc_chips/modules/scripts/align_readsPerChrom.sh -a {input.bam} > {output} 2>> {log}"
