#MODULE: CEAS- annotating where the peaks fall (promoter, exon, intron, interg)
# _logfile=output_path + "/logs/ceas.log"

#NOTE: using the _refs from chips.snakefile
def ceas_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config["runs"].keys():
        for rep in _reps[run]:
            #GENERATE Run name: concat the run and rep name
            runRep = "%s.%s" % (run, rep)
            ls.append(output_path + "/ceas/%s/%s_summary.txt" % (runRep,runRep))
            ls.append(output_path + "/ceas/%s/%s_DHS_peaks.bed" % (runRep,runRep))
            ls.append(output_path + "/ceas/%s/%s_DHS_stats.txt" % (runRep,runRep))
            ls.append(output_path + "/ceas/%s/%s_DHS_summary.dhs" % (runRep,runRep))
            if config["velcro_regions"]:
                ls.append(output_path + "/ceas/%s/%s_velcro_peaks.bed" % (runRep,runRep))
                ls.append(output_path + "/ceas/%s/%s_velcro_stats.txt" % (runRep,runRep))
    #ADD bam_regionStats
    for sample in config["samples"]:
        if config['exons']:
            ls.append(output_path + "/ceas/samples/%s/%s.exons" % (sample,sample))
        if config['promoters']:
            ls.append(output_path + "/ceas/samples/%s/%s.promoters" % (sample,sample))
        if config['DHS']:
            ls.append(output_path + "/ceas/samples/%s/%s.DHS" % (sample,sample))
        if config['exons'] and config['promoters'] and config['DHS']:
            ls.append(output_path + "/ceas/samples/%s/%s_meta.json" % (sample,sample))
    ls.append(output_path + "/ceas/samples/bamRegionStats.csv")
    ls.append(output_path + "/ceas/dhs.csv")
    ls.append(output_path + "/ceas/meta.csv")
    return ls

def collect_BamRegionStats_dirs(file_paths):
    """Given a list of file paths, returns the dirname of the filepaths"""
    #NOTE: relying on os to be imported in chips.snakefile
    dirs = [os.path.dirname(f) for f in file_paths]
    return [" -d %s" % d for d in dirs]

def ceasInput(wildcards):
    run = wildcards.run
    rep = wildcards.rep
    runRep = "%s.%s" % (run,rep)
    if ("macs2_broadpeaks" in config) and config["macs2_broadpeaks"]:
        temp = output_path + "/peaks/%s/%s_peaks.bed" % (runRep,runRep)
    else:
        temp = output_path + "/peaks/%s/%s_summits.bed" % (runRep,runRep)
    return temp

rule ceas_all:
    input:
        ceas_targets

rule ceas_annotatePeaksRegions:
    """Annotate peak regions"""
    input:
        # output_path + "/peaks/{run}.{rep}/{run}.{rep}_summits.bed"
        # output_path + "/peaks/{run}.{rep}/{run}.{rep}_peaks.bed"
        ceasInput
    output:
        promoter=output_path + "/ceas/{run}.{rep}/{run}.{rep}_summits_promoter.bed",
        exon=output_path + "/ceas/{run}.{rep}/{run}.{rep}_summits_exon.bed",
        intron=output_path + "/ceas/{run}.{rep}/{run}.{rep}_summits_intron.bed",
        intergenic=output_path + "/ceas/{run}.{rep}/{run}.{rep}_summits_intergenic.bed",
        summary=output_path + "/ceas/{run}.{rep}/{run}.{rep}_summary.txt",
    message: "CEAS: annotating peak regions"
    log: output_path + "/logs/ceas/{run}.{rep}.log"
    benchmark: output_path + "/Benchmark/{run}.{rep}_ceas_annotatePeaksRegions.benchmark"
    conda: "../envs/ceas/ceas.yaml"
    params:
        db=config['geneTable'],
        path=output_path + "/ceas/{run}.{rep}/",
        name= "{run}.{rep}_summits",
    shell:
        #TWO ways to run bedAnnotate: w/ basename param (-n) or w/o
        #For now we keep the -n explictly defined
        "cidc_chips/modules/scripts/ceas_bedAnnotate.v2.py -g {params.db} -b {input} -o {params.path} -n {params.name} > {output.summary} 2>>{log}"
        #"cidc_chips/modules/scripts/bedAnnotate.v2.py -g {params.db} -b {input} -o {params.path} > {output.summary} 2>>{log}"


rule ceas_takeTop5k:
    """Take the top 5000 sites"""
    input:
        output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.bed"
    params:
        n=5000
    message: "DHS: Take top sites"
    log: output_path + "/logs/ceas/{run}.{rep}.log"
    conda: "../envs/ceas/ceas.yaml"
    output:
        temp(output_path + '/ceas/{run}.{rep}/{run}.{rep}_sorted_5k_peaks.bed')
    shell:
        "head -n {params.n} {input} > {output} " #2>>{log}"

#------------------------------------------------------------------------------
rule ceas_DHSIntersect:
    """Intersect PEAKS with DHS regions"""
    input:
        output_path + '/ceas/{run}.{rep}/{run}.{rep}_sorted_5k_peaks.bed'
    params:
        #check for config['DHS'] defined, otherwise, use null
        dhs=config['DHS'] if config['DHS'] else "/dev/null"
    message: "DHS: intersect PEAKS with DHS regions"
    log: output_path + "/logs/ceas/{run}.{rep}.log"
    conda: "../envs/ceas/ceas.yaml"
    output:
        output_path + '/ceas/{run}.{rep}/{run}.{rep}_DHS_peaks.bed'
    shell:
        "intersectBed -wa -u -a {input} -b {params.dhs} > {output} " #2>>{log}"
    #run:
        #HACK! IN the CASE where no DHS is defined/available for species/assemb
        #USE AN EMPTY FILE.  CONSEQUENCE- everything downstream will count 0
        #DHS regions, NOT N/A
    #    if config['DHS']:
    #        shell("intersectBed -wa -u -a {input} -b {params.dhs} > {output} 2>>{log}")
    #    else:
    #        #make empty file
    #        shell("touch .snakemake/null.dhs.txt")
    #        shell("intersectBed -wa -u -a {input} -b .snakemake/null.dhs.txt > {output} 2>>{log}")

rule ceas_DHSStat:
    """collect DHS stats"""
    input:
        n=output_path + '/ceas/{run}.{rep}/{run}.{rep}_sorted_5k_peaks.bed',
        dhs=output_path + '/ceas/{run}.{rep}/{run}.{rep}_DHS_peaks.bed'
    message: "DHS: collecting stats"
    log: output_path + "/logs/ceas/{run}.{rep}.log"
    conda: "../envs/ceas/ceas.yaml"
    output:
        output_path + '/ceas/{run}.{rep}/{run}.{rep}_DHS_stats.txt'
    shell:
        "wc -l {input.n} {input.dhs} > {output} " #2>>{log}"

rule ceas_DHSSummary:
    """get DHS summary"""
    input:
        output_path + '/ceas/{run}.{rep}/{run}.{rep}_DHS_stats.txt'
    output:
        output_path + '/ceas/{run}.{rep}/{run}.{rep}_DHS_summary.dhs'
    message: "DHS: collecting stats"
    log: output_path + "/logs/ceas/{run}.{rep}.log"
    conda: "../envs/ceas/ceas.yaml"
    shell:
        """cat {input} | awk '{{printf"%s,", $1}}' > {output}"""

#------------------------------------------------------------------------------
rule ceas_VELCROIntersect:
    """Intersect PEAKS with velcro regions"""
    input:
        output_path + '/ceas/{run}.{rep}/{run}.{rep}_sorted_5k_peaks.bed'
    params:
        #CHECK for if config is set, otherwise use /dev/null
        velcro=config['velcro_regions'] if config['velcro_regions'] else "/dev/null"
    message: "VELCRO: intersect PEAKS with velcro regions"
    log: output_path + "/logs/ceas/{run}.{rep}.log"
    conda: "../envs/ceas/ceas.yaml"
    output:
        output_path + '/ceas/{run}.{rep}/{run}.{rep}_velcro_peaks.bed'
    shell:
        #NOTE: if no velcro defined, then maybe we should print out warning
        "intersectBed -wa -u -a {input} -b {params.velcro} > {output} " # 2>>{log}"
    #run:
    #    #CHECK for the existence of this file!
    #    if params.velcro:
    #        shell("intersectBed -wa -u -a {input} -b {params.velcro} > {output} 2>>{log}")
    #    else:
    #        #No velcro file defined --> empty output
    #        shell("touch {output} && echo 'WARNING: no velcro region defined' >>{log}")

rule ceas_VELCROStat:
    """collect VELCRO stats"""
    input:
        n=output_path + '/ceas/{run}.{rep}/{run}.{rep}_sorted_5k_peaks.bed',
        velcro=output_path + '/ceas/{run}.{rep}/{run}.{rep}_velcro_peaks.bed'
    message: "VELCRO: collecting stats"
    log: output_path + "/logs/ceas/{run}.{rep}.log"
    conda: "../envs/ceas/ceas.yaml"
    output:
        output_path + '/ceas/{run}.{rep}/{run}.{rep}_velcro_stats.txt'
    shell:
        "wc -l {input.n} {input.velcro} > {output}" # 2>>{log}"

rule ceas_bamRegionStat:
    """count the number of reads in promoter, exon, dhs--these regions
    are defined in the config.yaml"""
    input:
        output_path + "/align/{sample}/{sample}_4M_unique_nonChrM.bam"
    params:
        bed = lambda wildcards: config[wildcards.region],
        #for use in message only
        msg = lambda wildcards: "%s:%s" % (wildcards.sample, wildcards.region)
    message: "CEAS: bam stat region {params.msg}"
    log: output_path + "/logs/ceas/{sample}.{region}.log"
    benchmark: output_path + "/Benchmark/{sample}.{region}_ceas_bamRegionStat.benchmark"
    conda: "../envs/ceas/ceas.yaml"
    output:
        #make temp
        output_path + '/ceas/samples/{sample}/{sample}.{region}'
    shell:
        "cidc_chips/modules/scripts/ceas_meta_bamRegionCount.sh -i {input} -b {params.bed} -o {output}" # 2>> {log}"

rule ceas_collectBamRegionStatsToJson:
    """collect the BAM region stats into a single file"""
    input:
        #INPUT the stats directories--
        #hack just add all of the file so we dont get a missing input exception
        #and collect the directories down below
        dhs = output_path + "/ceas/samples/{sample}/{sample}.DHS",
        prom = output_path + "/ceas/samples/{sample}/{sample}.promoters",
        exon = output_path + "/ceas/samples/{sample}/{sample}.exons",
    message: "CEAS: collect bam region stats into json file"
    log: output_path + "/logs/ceas/{sample}.log"
    benchmark: output_path + "/Benchmark/{sample}_ceas_collectBamRegionStatsToJson.benchmark"
    conda: "../envs/ceas/ceas.yaml"
    output:
        output_path + "/ceas/samples/{sample}/{sample}_meta.json"
    shell:
        "cidc_chips/modules/scripts/ceas_collectBamRegStatsIntoJson.py -d {input.dhs} -p {input.prom} -e {input.exon} -o {output}"

rule ceas_collectBamRegionStats:
    """collect the BAM region stats into a single file"""
    input:
        #INPUT the stats directories--
        #hack just add all of the file so we dont get a missing input exception
        #and collect the directories down below
        dhs = expand(output_path + "/ceas/samples/{sample}/{sample}.DHS", sample=config['samples']),
        prom = expand(output_path + "/ceas/samples/{sample}/{sample}.promoters", sample=config['samples']),
        exon = expand(output_path + "/ceas/samples/{sample}/{sample}.exons", sample=config['samples'])
    params:
        dirs = lambda wildcards,input: collect_BamRegionStats_dirs(input.exon)
    message: "CEAS: collect bam region stats"
    # log: output_path + "/logs/ceas/{run}.{rep}.log"
    conda: "../envs/ceas/ceas.yaml"
    output:
        output_path + '/ceas/samples/bamRegionStats.csv'
    shell:
        "cidc_chips/modules/scripts/ceas_collectBamRegStats.py {params.dirs} > {output}" # 2>>{log}"

rule ceas_collectDHSstats:
    """collect the DHS stats into a single file"""
    input:
        #Generalized INPUT fn defined in chips.snakefile
        _getRepInput(output_path + "/ceas/$runRep/$runRep_DHS_stats.txt")
    params:
        files = lambda wildcards, input: [" -f %s" % i for i in input]
    message: "CEAS: collect DHS stats"
    # log: output_path + "/logs/ceas/{run}.{rep}.log"
    conda: "../envs/ceas/ceas.yaml"
    output:
        output_path + '/ceas/dhs.csv'
    shell:
        "cidc_chips/modules/scripts/ceas_peaks_getDHSstats.py {params.files} -o {output} " # 2>>{log}"

rule ceas_collectCEASstats:
    """collect the CEAS stats into a single file"""
    input:
        #Generalized INPUT fn defined in chips.snakefile
        _getRepInput(output_path + "/ceas/$runRep/$runRep_summary.txt")
    params:
        files = lambda wildcards, input: [" -f %s" % i for i in input]
    message: "CEAS: collect CEAS stats"
    #log: output_path + "/logs/ceas/{run}.{rep}.log"
    conda: "../envs/ceas/ceas.yaml"
    output:
        output_path + '/ceas/meta.csv'
    shell:
        "cidc_chips/modules/scripts/ceas_getMetaStats.py {params.files} -o {output} " #2>>{log}"
