import yaml
configfile: "chips_assembler.yaml"

def get_filename(path, withoutExt=False):
    """GIVEN a path, returns the filename, i.e. last element of path"""
    f = path.split("/")[-1]
    if withoutExt:
        #REMOVE the .ext
        return ".".join(f.split(".")[:-1])
    else:
        return f

def calcGenomeSize(chrom_lens_file):
    """Given a chrom_lens_file, returns in scientific notation the 
    genome size"""
    ssum = 0
    f = open(chrom_lens_file)
    for l in f:
        tmp = l.strip().split()
        ssum += int(tmp[1])
    return("%.2e" %ssum)

def writeRef(name, bwa, bwt, ln, refGene, prom, exon, conserv, out_file):
    """Given paths to the various elms, writes out the ref.yaml file"""
    #PROCESS the input down
    #bwa = ref_files/mm10/bwa_indices/mm10.fa.amb --CHOP off .amb
    tmp = get_filename(bwa, withoutExt=True)
    bwa = os.path.join("/".join(bwa.split("/")[:-1]), tmp)    
    #bwt = ref_files/mm10/bowtie2/mm10.1.bt2 --chomp off .1.bt2
    bwt = os.path.join("/".join(bwt.split("/")[:-1]), name)
    #conservation = ref_files/mm10/conservation/convert_done.txt --drop file

    conserv = "/".join(conserv.split("/")[:-1])
    out_f = open(out_file, "w")
    #BUILD dict
    tmp = {'bwa_index': bwa}
    tmp['bwt2_index'] = bwt
    tmp['geneTable'] = refGene
    tmp['conservation'] = conserv
    tmp['DHS'] = ""
    tmp['exons'] = exon
    tmp['promoters'] = prom
    tmp['velcro_regions'] = ""
    tmp['chrom_lens'] = ln
    #NOTE: IF not defined, calculate it automatically from chrom lens file
    if 'genome_size' in config and config['genome_size']:
        tmp['genome_size'] = config['genome_size']
    else:
        tmp['genome_size'] = calcGenomeSize(ln)
    tmp['motif_path'] = name

    #write to file
    out_f.write(yaml.dump({name:tmp}, default_flow_style=False))
    out_f.close()

assembly = config['assembly_name']
assembly_fname = get_filename(config['fasta'])

conservation_fname = get_filename(config['conservation_bw'], withoutExt=True)

rule all:
    input:
        f="%s.ref.yaml" % (assembly)

rule bwa_index:
    "GENERATE the bwa_index"
    input:
        fa = config['fasta']
    params:
        assembly=config['assembly_name'],
        fname = assembly_fname
    output:
        "ref_files/%s/bwa_indices/%s.amb" % (assembly, assembly_fname)
    shell:
        #NOTE: BWA generates the files in the local directory--so we to do this
        #HACK to make sure everthing is made in the right place!  GRRR!
        "cd ref_files/{params.assembly}/bwa_indices/ && "
        "ln -s {input} && """ #to bring the file locally; avoid permission errs
        "bwa index {params.fname} 2> bwa_index.log.out"

rule bowtie2_index:
    "Generate the bowtie2 index"
    input:
        fa = config['fasta']
    params:
        outpath = "ref_files/%s/bowtie2/%s" % (assembly, assembly)
    output:
        "ref_files/%s/bowtie2/%s.1.bt2" % (assembly, assembly)
    shell:
        """bowtie2-build {input} {params.outpath}"""

rule getChromLen:
    """Fetches the genome chromosome size info.  Uses UCSC's fetchChromSizes,
    removes the chromosomes that are not standardized, e.g. chr4_random and
    sorts into KARYOTYPIC order"""
    params:
        assembly = assembly,
        grep_pat = "\"\\_\""
    output:
        "ref_files/%s/regions/%s.len" % (assembly, assembly)
    shell:
        "chips/modules/scripts/assembler/fetchChromSizes {params.assembly} | "
        "grep -v {params.grep_pat} | sort -k1,1 > {output}"

rule refGene:
    "Given a refGene file, copy it to the right place"
    input:
        config['refGene']
    output:
        "ref_files/%s/%s.refGene" % (assembly, get_filename(assembly))
    shell:
        "cp {input} {output}"


rule generateExonPromoter:
    """GIVEN a refGene, uses a script to parse out the exons and promoters"""
    input:
        refGene=config['refGene'],
        lengths="ref_files/%s/regions/%s.len" % (assembly, assembly)
    params:
        outpath="ref_files/%s/regions/" % (assembly),
        name = assembly
    output:
        "ref_files/%s/regions/%s.promoter.bed" % (assembly, assembly),
        "ref_files/%s/regions/%s.exon.bed" % (assembly, assembly)
    shell:
        "chips/modules/scripts/assembler/parseRefGene.py -g {input.refGene} -l {input.lengths} -o {params.outpath} -n {params.name}"

#CONSERVATION:
rule conservation_bwToBedGraph:
    """Convert the bw file to a bedGraph"""
    input:
        config['conservation_bw']
    output:
        "ref_files/%s/conservation/%s.bedGraph" % (assembly,conservation_fname)
    shell:
        "chips/modules/scripts/assembler/bigWigToBedGraph {input} {output}"

#NOTE: non-canonical snakemake--where we just look for one file at the end of
#the script
rule conservation_splitBedGraph:
    """Split the large bedgraph into individual chroms"""
    input:
        "ref_files/%s/conservation/%s.bedGraph" % (assembly,conservation_fname)
    params:
        outdir="ref_files/%s/conservation/tmp/" % (assembly)
    output:
        #FILE to test the script is done
        "ref_files/%s/conservation/tmp/split_done.txt" % (assembly),
        "ref_files/%s/conservation/tmp/chroms.txt" % (assembly)
    shell:
        "chips/modules/scripts/assembler/conservation_splitBedgraph.sh {input} {params.outdir}"

#NOTE: non-canonical snakemake--where we just look for one file at the end of
#the script
rule conservation_convertToBw:
    """Convert the individual .bedGraph files back to .bw"""
    input:
        splitBedGraph_flag="ref_files/%s/conservation/tmp/split_done.txt" % (assembly),
        length="ref_files/%s/regions/%s.len" % (assembly, assembly)
    params:
        indir="ref_files/%s/conservation/tmp/" % (assembly),
        outdir="ref_files/%s/conservation/" % (assembly)
    output:
        "ref_files/%s/conservation/chr1.bw" % (assembly),
        "ref_files/%s/conservation/convert_done.txt" % (assembly)
    shell:
        "chips/modules/scripts/assembler/conservation_convertToBw.sh {params.indir} {params.outdir} {input.length}"

rule write_ref_yaml:
    """When we finish generating all of the files, generate a ref.yaml for the
    species to integrate into ref.yaml"""
    input:
        bwa="ref_files/%s/bwa_indices/%s.amb" % (assembly, assembly_fname),
        bwt="ref_files/%s/bowtie2/%s.1.bt2" % (assembly, assembly),
        ln="ref_files/%s/regions/%s.len" % (assembly, assembly),
        refGene="ref_files/%s/%s.refGene" % (assembly, get_filename(assembly)),
        prom="ref_files/%s/regions/%s.promoter.bed" % (assembly, assembly),
        exon="ref_files/%s/regions/%s.exon.bed" % (assembly, assembly),
        conserv="ref_files/%s/conservation/convert_done.txt" % (assembly)
    output:
        f="%s.ref.yaml" % (assembly)
    run:
        writeRef(assembly, input.bwa, input.bwt, input.ln, input.refGene, input.prom, input.exon, input.conserv, output.f)


#TODO: 
#1. write a cleanup script to delete conservation/tmp and the large bedGraph 
    
