#MODULE: bam_snapshots- takes bam read pile up snapshots for a set of genes
# _logfile=output_path + "/logs/bam_snapshot.log"

#LIST of target genes and zoom factor
_gene_list=[("ACTB",3),
            ("KLK3",3),
            ("TMPRSS2",3),
            ("GAPDH", 10)]

        
def bam_snapshots_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for sample in config["samples"]:
        for (gene, zoom_f) in _gene_list:
            ls.append(output_path + "/bam_snapshots/%s/%s_%s_%sXzoom.png" % (sample,sample,gene,zoom_f))
    return ls

rule bam_snapshots_all:
    input:
        bam_snapshots_targets

rule bam_snapshot:
    """Generates a snapshot of a sample's read pileup in given gene"""
    input:
        output_path + "/align/{sample}/{sample}.sorted.bam"
    output:
        output_path + "/bam_snapshots/{sample}/{sample}_{gene}_{zoom_f}Xzoom.png"
    message: "BAM_SNAPSHOT: generate snapshot for sample {wildcards.sample} in gene {wildcards.gene}"
    # log: output_path + "/logs/bam_snapshot/{sample}.log"
    conda: "../envs/bam_snapshots/bam_snapshots.yaml"
    params:
        isPaired= lambda wildcards: "TRUE" if len(config["samples"][wildcards.sample]) == 2 else "FALSE",
        name= lambda wildcards: wildcards.sample,
        gene= lambda wildcards: wildcards.gene,
        species= config['motif_path'],
        zoom= lambda wildcards: wildcards.zoom_f
    shell:
        "Rscript cidc_chips/modules/scripts/bam_snapshot.R {input} {params.isPaired} {params.name} {params.gene} {params.species} {params.zoom} {output}"
