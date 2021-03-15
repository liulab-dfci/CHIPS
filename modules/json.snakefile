#MODULE: json- generating the output_path/json folder for cistrome database
# _logfile=output_path + "/logs/json.log"

def json_targets(wildcards):
    ls = []
    # ls.append(output_path + "/json")
    for run in config["runs"]:
        # ls.append(output_path + "/json/%s/" % run)
        ls.append(output_path + "/json/%s/%s_conserv.json" % (run, run))
        if config["DHS"]:
            ls.append(output_path + "/json/%s/%s_dhs.json" % (run, run))        
        # if config["velcro_regions"]:
        #     ls.append(output_path + "/json/%s/%s_velcro.json" % (run, run))
        ls.append(output_path + "/json/%s/%s_frip.json" % (run, run))
        ls.append(output_path + "/json/%s/%s_macs2.json" % (run, run))
        if config['runs'][run][2]:
            ls.append(output_path + "/json/%s/%s_macs2_rep.json" % (run, run))
            # ls.append(output_path + "/json/%s/%s_rep.json" % (run, run))
        ls.append(output_path + "/json/%s/%s_meta.json" % (run, run))
        ls.append(output_path + "/json/%s/%s_fastqc.json" % (run, run))
        ls.append(output_path + "/json/%s/%s_map.json" % (run, run))
        ls.append(output_path + "/json/%s/%s_enrich_meta.json" % (run, run))
        ls.append(output_path + "/json/%s/%s_pbc.json"% (run, run))
        ls.append(output_path + "/json/%s/%s_frag.json" % (run, run))
        if ("macs2_broadpeaks" not in config) or config["macs2_broadpeaks"] == False:
            if ("motif" in config) and config["motif"] == "mdseqpos":
                ls.append(output_path + "/json/%s/%s_seqpos.json" % (run, run))
        # Cistrome DB do not need contamination part for now
        # run_samples = config['runs'][run]
        # for sample in run_samples:
        #     if sample:
        #         #contam is sample-- see cistrome.snakefile- fastqc as example
        #         ls.append(output_path + "/json/%s/%s_contam.json" % (run, sample))

    return ls

def json_getRunAndRep(wildcards):
    # item = 0
    # replicate = []
    # for i in range(2,len(r)+2,2):
    #     a = int(i/2)
    #     replicate.append("rep%s" % str(a))
    #     replicate.append("rep%s" % str(a))
    #     item += 2
    # # print(replicate)
    # sample = wildcards.sample
    # # print(sample)
    # key_list=[]
    # value_list=[]
    # for key,value in config["runs"].items():
    #     key_list.append(key)
    #     value_list.append(value)
    # for i in value_list:
    #     if sample in i:
    #         run = key_list[value_list.index(i)]
    #         rep = replicate[i.index(sample)]
    #         break
    #     else:
    #         continue
    # return "%s.%s" % (run, rep)

    # assume chips's {run}.{rep1} = cistrome {run}
    run_rep = "%s.rep1" % wildcards.run
    return run_rep

def json_getSample(wildcards):
    json_sample = config['runs'][wildcards.run]
    return json_sample

def json_checkBAMPE(wildcards):
    """Fn returns '-f BAMPE' IF the run's FIRST treatment replicate (sample) is
    Paired-END.
    NOTE: this flag is required for macs2 callpeak, AUTO detect format does not
    work with PE bams!
    """
    r = config['runs'][wildcards.run]
    #GET first treatement sample
    first = config['samples'][r[0]]
    ret = "-f BAMPE" if len(first) == 2 else "-f BAMSE"
    return ret

# def velcro_target(wildcards):
#     input = output_path + "/ceas/%s/%s_velcro_summary.txt" % (json_getRunAndRep(wildcards),json_getRunAndRep(wildcards))
#     return input

rule json_all:
    input:
        json_targets

rule json_conservation:
    input:
        lambda wildcards:output_path + "/conserv/%s/%s_conserv.txt" % (json_getRunAndRep(wildcards),json_getRunAndRep(wildcards))
    output:
        json=output_path + "/json/{run}/{run}_conserv.json",
        # dir=directory(output_path + "/json/{run}/")
    params:
        basics = "-b %s " % config["basics"] if "basics" in config else "",
        factor = "-f %s " % config["factor"] if "factor" in config else "",
        TF = "-T %s " % config["TF"] if "TF" in config else ""
    message: "JSON: generate conservation json"
    # log: output_path + "/logs/json/{sample}.log"
    shell:
        "cidc_chips/modules/scripts/json/json_conserv.py -i {input} -o {output.json} {params.basics} {params.factor} {params.TF} -r {wildcards.run}"


# rule json_comtamination:
#     input:
#         output_path + "/contam/{sample}/{sample}_contamination.txt"
#     output:
#         output_path + "/json/{run}/{sample}_contam.json"
#     params:

#     message: "JSON: generate comtamination json"
#     log: _logfile
#     shell:
#         "cidc_chips/modules/scripts/json/json_comtamination.py -i {input} -o {output} -I {wildcards.sample} -s {wildcards.sample}"


rule json_dhs:
    input:
        lambda wildcards:output_path + "/ceas/%s/%s_DHS_summary.dhs" % (json_getRunAndRep(wildcards),json_getRunAndRep(wildcards))
    output:
        output_path + "/json/{run}/{run}_dhs.json"
    message: "JSON: generate DHS json"
    # log: output_path + "/logs/json/{sample}.log"
    shell:
        "cidc_chips/modules/scripts/json/json_dhs.py -i {input} -o {output}"

# rule json_velcro:
#     input:
#         velcro_target
#     output:
#         output_path + "/json/{sample}/{sample}_velcro.json"
#     params:

#     message: "JSON: generate velcro json"
    # log: output_path + "/logs/json/{sample}.log"
#     shell:
#         "cidc_chips/modules/scripts/json/json_velcro.py  "

rule json_enrichMeta:
    input:
        meta = lambda wildcards: output_path + "/ceas/samples/%s/%s_meta.json" % (json_getSample(wildcards)[0], json_getSample(wildcards)[0]),
        mapped = lambda wildcards: output_path + "/align/%s/%s_mapping.txt" % (json_getSample(wildcards)[0], json_getSample(wildcards)[0]),
        dhs = lambda wildcards: output_path + "/ceas/%s/%s_DHS_summary.dhs" % (json_getRunAndRep(wildcards),json_getRunAndRep(wildcards))
    output:
        output_path + "/json/{run}/{run}_enrich_meta.json"
    params:
        mapped_files = lambda wildcards, input: [" -m %s" % i for i in {input.mapped}],
        has_dhs = "-H %s " % config["DHS"],
        down = True,
        sample_name = lambda wildcards: [" -s %s" % i for i in json_getSample(wildcards) if i] if json_getSample(wildcards) else "",
    message: "JSON: generate meta enrichment json"
    # log: output_path + "/logs/json/{sample}.log"
    shell:
        "cidc_chips/modules/scripts/json/json_enrich_meta.py -i {input.meta} -o {output} {params.mapped_files} -D {input.dhs} {params.has_dhs} -d {params.down} -r {run} {params.sample_name}"


rule json_fastqc:
    input:
        lambda wildcards: output_path + "/fastqc/%s/%s_100k_fastqc/fastqc_data.txt" % (json_getSample(wildcards)[0], json_getSample(wildcards)[0])
    output:
        output_path + "/json/{run}/{run}_fastqc.json"
    params:
        fastqc_data = lambda wildcards: [" -i " + output_path + "/fastqc/%s/%s_100k_fastqc/fastqc_data.txt" % (i,i) for i in json_getSample(wildcards) if i] if json_getSample(wildcards) else "",
        sample_name = lambda wildcards: [" -s %s" % i for i in json_getSample(wildcards) if i] if json_getSample(wildcards) else "",
    message: "JSON: generate fastqc json"
    # log: output_path + "/logs/json/{sample}.log"
    shell:
        "cidc_chips/modules/scripts/json/json_fastqc.py {params.fastqc_data} -o {output}{params.sample_name} -r {run}"


rule json_frag:
    input:
        lambda wildcards: output_path + "/peaks/%s/%s_model.R" % (json_getRunAndRep(wildcards),json_getRunAndRep(wildcards))
    output:
        json = output_path + "/json/{run}/{run}_frag.json",
    params:
        sd_R = lambda wildcards: output_path + "/peaks/%s/%s_model_sd.R" % (json_getRunAndRep(wildcards),json_getRunAndRep(wildcards)),
        frag_tool = lambda wildcards: "%s" % json_checkBAMPE(wildcards),
        sample_name = lambda wildcards: [" -s %s" % i for i in json_getSample(wildcards) if i] if json_getSample(wildcards) else "",
    message: "JSON: generate frag json"
    # log: output_path + "/logs/json/{sample}.log"
    shell:
        "cidc_chips/modules/scripts/json/json_frag.py -r {input} -o {output.json} -R {params.sd_R} {params.frag_tool} {params.sample_name}"


rule json_frip:
    input:
        lambda wildcards: output_path + "/frips/%s/%s_frip.txt" % (json_getRunAndRep(wildcards),json_getRunAndRep(wildcards))
    output:
        output_path + "/json/{run}/{run}_frip.json"
    message: "JSON: generate frip json"
    params:
        input_file = lambda wildcards, input: [" -i %s" % i for i in input],
        sample_name = lambda wildcards: [" -s %s" % i for i in json_getSample(wildcards) if i] if json_getSample(wildcards) else "",
    # log: output_path + "/logs/json/{sample}.log"
    shell:
        "cidc_chips/modules/scripts/json/json_frip.py {params.input_file} -o {output}{params.sample_name}"


rule json_macs2:
    input:
        lambda wildcards: output_path + "/peaks/%s/%s_peaks.xls" % (json_getRunAndRep(wildcards),json_getRunAndRep(wildcards))
    output:
        output_path + "/json/{run}/{run}_macs2.json"
    message: "JSON: generate macs2 json"
    # log: output_path + "/logs/json/{sample}.log"
    shell:
        "cidc_chips/modules/scripts/json/json_macs2.py -i {input} -o {output} -I {wildcards.run}"


rule json_macs2Rep:
    input:
        lambda wildcards: output_path + "/peaks/%s/%s_peaks.xls" % (json_getRunAndRep(wildcards),json_getRunAndRep(wildcards))
    output:
        output_path + "/json/{run}/{run}_macs2_rep.json"
    params:
        sample_name = lambda wildcards: [" -s %s" % i for i in json_getSample(wildcards) if i] if json_getSample(wildcards) else "",
    message: "JSON: generate macs2 rep json"
    # log: output_path + "/logs/json/{sample}.log"
    shell:
        "cidc_chips/modules/scripts/json/json_macs2_rep.py -i {input} -o {output} {params.sample_name} "


rule json_map:
    input:
        lambda wildcards: output_path + "/align/%s/%s_mapping.txt" % (json_getSample(wildcards)[0], json_getSample(wildcards)[0])
    output:
        output_path + "/json/{run}/{run}_map.json"
    message: "JSON: generate map json"
    params:
        input_file = lambda wildcards: [" -i " + output_path + "/align/%s/%s_mapping.txt" % (i,i) for i in json_getSample(wildcards) if i] if json_getSample(wildcards) else "",
    # log: output_path + "/logs/json/{sample}.log"
    shell:
        "cidc_chips/modules/scripts/json/json_map.py{params.input_file} -o {output} "


rule json_meta:
    input:
        lambda wildcards: output_path + "/ceas/%s/%s_summary.txt" % (json_getRunAndRep(wildcards),json_getRunAndRep(wildcards))
    output:
        output_path + "/json/{run}/{run}_meta.json"
    message: "JSON: generate meta json"
    # log: output_path + "/logs/json/{sample}.log"
    shell:
        "cidc_chips/modules/scripts/json/json_meta.py -i {input} -o {output} -I {wildcards.run} "


rule json_pbc:
    input:
        lambda wildcards: output_path + "/frips/%s/%s_pbc.txt" % (json_getSample(wildcards)[0], json_getSample(wildcards)[0])
    output:
        output_path + "/json/{run}/{run}_pbc.json"
    message: "JSON: generate pbc json"
    params:
        input_file = lambda wildcards: [" -i " + output_path + "/frips/%s/%s_pbc.txt" % (i,i) for i in json_getSample(wildcards) if i] if json_getSample(wildcards) else "",
    # log: output_path + "/logs/json/{sample}.log"
    shell:
        "cidc_chips/modules/scripts/json/json_pbc.py{params.input_file} -o {output}"

rule json_seqpos:
    input:
        lambda wildcards: output_path + "/motif/%s/results/motif_list.json" % json_getRunAndRep(wildcards)
    output:
        output_path + "/json/{run}/{run}_seqpos.json"
    message: "JSON: generate seqpos json"
    params:
        prefix = lambda wildcards: output_path + "/motif/%s/results/seqLogo/" % json_getRunAndRep(wildcards)
    # log: output_path + "/logs/json/{sample}.log"
    shell:
        "cidc_chips/modules/scripts/json/json_seqpos.py -i {input} -o {output} -p {params.prefix}"
# rule json_rep:
#     input:
        
#     output:
#         output_path + "/json/{run}/{run}_rep.json"
#     message: "JSON: generate rep json"
#     log: output_path + "/logs/json/{sample}.log"
#     shell:
#         "cidc_chips/modules/scripts/json/json_rep.py -c {input.cor} -o {output} -O {input.overlap} -i {run}" 


