#MODULE: mapmaker- generating the analysis/mapmaker folder
#The folder should be a fully runnable version of the latest mapmaker USING
#the chips run information
_logfile=output_path + "/logs/mapmaker.log"

def mapmaker_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    # ls.append(output_path + "/mapmaker")
    ls.append(output_path + "/mapmaker/config.yaml")
    ls.append(output_path + "/mapmaker/metasheet.csv")
    return ls

def make_relative(path):
    """GIVEN a path that is relative to the PROJECT directory, e.g.
    analysis/peaks/{run}.rep1/.., we convert the path to be relative to the 
    analysis/mapmaker directory.  In other words, we replace analysis with ..
    """
    return path.replace("analysis","..")

def getTreat(run):
    """GIVEN a run, will return the sample name of the first treatment"""
    treat1 = config['runs'][run][0]
    return treat1

rule mapmaker_all:
    input:
        mapmaker_targets

rule mapmaker_clone:
    """Clones the mapmaker project"""
    message: "MAPMAKER: cloning the mapmaker project from github"
    log: _logfile
    output:
        # output_path + "/mapmaker",
        output_path + "/mapmaker/.config.yaml",
        output_path + "/mapmaker/.metasheet.csv"
    conda: "../envs/mapmaker/mapmaker.yaml"
    shell:
        # "git clone git@bitbucket.org:cfce/mapmaker analysis/mapmaker > {log} 2>&1 && "
        "git clone https://bitbucket.org/cfce/mapmaker.git analysis/mapmaker > {log} 2>&1 && "
        "sleep 20 && "
        #NOTE: the two files below are OTHER targets of this snakefile, so
        #we trick snakemake into thinking they don't exist by renaming them
        "mv analysis/mapmaker/config.yaml analysis/mapmaker/.config.yaml && "
        "mv analysis/mapmaker/metasheet.csv analysis/mapmaker/.metasheet.csv"

rule mapmaker_config:
    """Tries to configure the mapmaker run based on the chips run info"""
    input:
        bed_files= expand(output_path + "/peaks/{run}.rep1/{run}.rep1_sorted_peaks.bed", run=config['runs']),
        bw_files= expand(output_path + "/peaks/{run}.rep1/{run}.rep1_treat_pileup.bw", run=config['runs']),
        bam_files = expand(output_path + "/align/{sample}/{sample}_unique.sorted.dedup.bam", sample=[getTreat(r) for r in config['runs']]),
        igv_files= expand(output_path + "/peaks/{run}.rep1/{run}.rep1_treatment.igv.xml", run=config['runs']),
        config=output_path + "/mapmaker/.config.yaml"
    params:
        # run_names=config['runs']
        bed_files = lambda wildcards, input: [" -b %s" % make_relative(b) for b in input.bed_files],
        bw_files  = lambda wildcards, input: [" -w %s" % make_relative(w) for w in input.bw_files],
        bam_files = lambda wildcards, input: [" -a %s" % make_relative(bam) for bam in input.bam_files],
        igv_files = lambda wildcards, input: [" -i %s" % make_relative(igv) for igv in input.igv_files],
        names = lambda wildcards : [" -n %s" % i for i in config['runs']]
    message: "MAPMAKER: configuring mapmaker"
    log: _logfile
    conda: "../envs/mapmaker/mapmaker.yaml"
    output:
        output_path + "/mapmaker/config.yaml"
    # run:
        #TRANSFORM to analysis/mapmaker relative
        #input.bed_files = [make_relative(b) for b in bed_files]
        #bed_files = " -b ".join([make_relative(b) for b in input.bed_files])
        # bw_files = " -w ".join([make_relative(bw) for bw in input.bw_files])
        # bam_files= " -a ".join([make_relative(bam) for bam in input.bam_files])
        # igv_files= " -i ".join([make_relative(igv) for igv in input.igv_files])
        # names = " -n ".join(params.run_names)
        # shell("cidc_chips/modules/scripts/mapmaker_config.py -n {names} -b {bed_files} -w {bw_files} -a {bam_files} -i {igv_files} -c {input.config} -o {output}")
    shell:
      "cidc_chips/modules/scripts/mapmaker_config.py {params.names} {params.bed_files} {params.bw_files} {params.bam_files} {params.igv_files} -c {input.config} -o {output}"  

rule mapmaker_meta:
    """Tries to configure the mapmaker meta based on the chips run info"""
    input:
        meta=output_path + "/mapmaker/.metasheet.csv"
    params:
        #run_names=config['runs']
        names = lambda wildcards : [" -n %s" % i for i in config['runs']]
    message: "MAPMAKER: configuring mapmaker"
    log: _logfile
    conda: "../envs/mapmaker/mapmaker.yaml"
    output:
        output_path + "/mapmaker/metasheet.csv"
    # run:
    #     names = " -n ".join(params.run_names)
    #     shell("cidc_chips/modules/scripts/mapmaker_meta.py -n {names} -o {output}")
    shell:
        "cidc_chips/modules/scripts/mapmaker_meta.py {params.names} -o {output}"
