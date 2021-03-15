# Author: Len Taing, Clara Cousins, Gali Bai
# Last modified: 01/11/2021
#MODULE: Chips report module

# Import packages
import yaml
from yaml import dump as yaml_dump
import pandas as pd
import os
import csv
import matplotlib.pyplot as plt
import numpy as np
import re
import glob



#define ref_bed to be used
ref_bed = "ref_files/%s/%s_refGene.bed" % (config['assembly'],config['assembly'])
#define prams from config.yaml file
output_path = config['output_path']
upstream = config["upstream"]
downstream = config["downstream"]
extend_bed = output_path + "/report/extend.bed"
tss_bed = output_path + "/report/tss.bed"

#function to generate bed files to be used for pygenometrack view
def pygenomeTacks_bed(ref_bed, extend_bed, tss_bed, upstream, downstream):
    #deduplicate based on gene symbol
    df_bed= pd.read_csv(ref_bed, sep = '\t',header=0, index_col=False).drop_duplicates(subset=['symbol'], keep='first')
    #extend the region of GenomeTrack view based on users' input
    df_bed["left_extend"] = df_bed["start"] - upstream
    df_bed["right_extend"] = df_bed["end"] + downstream
    extend_list = []
    tss_list = []
    for i in df_bed.index.tolist():
        chromosome = df_bed.loc[i, "chromosome"]
        start = df_bed.loc[i, "start"]
        end = df_bed.loc[i, "end"]
        if df_bed.loc[i, "left_extend"] < 0:
            istart = 0
        else:
            istart = df_bed.loc[i, "left_extend"]
        iend = df_bed.loc[i, "right_extend"]
        coordinate = '%s:%s-%s' % (chromosome, istart, iend)
        product_accession = df_bed.loc[i, "product_accession"]
        strand = df_bed.loc[i,"strand"]
        symbol = df_bed.loc[i, "symbol"]
        TSS = df_bed.loc[i, "TSS"]
        try:
            TSS_end = df_bed.loc[(i+1),"TSS"]
        except:
            TSS_end = df_bed.loc[i, "TSS"] + 1
        extend_list.append([chromosome, str(start), str(end), str(symbol), str(coordinate), str(product_accession), str(strand), str(TSS)])
        if TSS < TSS_end:
            tss_list.append([chromosome, str(TSS), str(TSS_end), str(coordinate), str(product_accession), str(strand), str(symbol), str(TSS)])
    print(extend_list[0:10])
    print(tss_list[0:10])
    f1 = open(extend_bed, "w")
    f1.write("\n".join(list(map(lambda x: "\t".join(x), extend_list))))
    f1.close()
    f2 = open(tss_bed, "w")
    f2.write("\n".join(list(map(lambda x: "\t".join(x), tss_list))))
    f2.close()


# Create input file lists
HOMERCONSERV_FILES_CONSERV_PROJECT = []
HOMERCONSERV_FILES_LOGOS_PROJECT = []
LOGOS_TXT_FILES = []
GENOMETRACK_FILES = []
for run in config["runs"].keys():
    for rep in _reps[run]:
        runRep = "%s.%s" % (run, rep)
        HOMERCONSERV_FILES_CONSERV_PROJECT.append(output_path + "/conserv/{runRep}/{runRep}_conserv_thumb.png".format(runRep = runRep))
        HOMERCONSERV_FILES_LOGOS_PROJECT.append(output_path + "/motif/{runRep}/results/knownResults/known1.logo.png".format(runRep = runRep))
        LOGOS_TXT_FILES.append(output_path + "/motif/{runRep}/results/knownResults.txt".format(runRep= runRep))
        GENOMETRACK_FILES.append(output_path + "/peaks/{runRep}/{runRep}_treat_pileup.bw".format(runRep = runRep))

# Write the name of the assembly into a csv to display in the report
File_object=open("cidc_chips/report/ref_assembly.csv",'w')
File_object.write(config["assembly"])
File_object.close()

#Define TrackView output files based on the genes to plot
TRACK_PNG_LIST = []
GENE_INCLUDED_LIST = []

for list_num, gene in enumerate(config["genes_to_plot"].strip().split()):
    if gene in pd.read_csv(ref_bed, sep = '\t',header=None, index_col=None).iloc[:,-2].values:
        GENE_INCLUDED_LIST.append(gene)
        TRACK_PNG_LIST.append((output_path + "/report/Genome_Track_View/01_genome_track_for_{track}.png").format(track = gene))
        print(gene + ' found in bed')
    else:
        print('Gene not in the lookup')


# Define wildcards for report
def report_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    ls.append(output_path + "/report/Overview/03_assembly.csv")
    ls.append(output_path + "/report/Overview/02_select_software_versions.dt")
    ls.append(output_path + "/report/Reads_Level_Quality/01_read_level_summary_table.dt")
    ls.append(output_path + "/report/Reads_Level_Quality/02_mapped_reads_bar.plotly")
    ls.append(output_path + "/report/Reads_Level_Quality/04_contamination_table.dt")
    ls.append(output_path + "/report/Reads_Level_Quality/05_contamination_bar.plotly")
    ls.append(output_path + "/report/Reads_Level_Quality/03_pcr_bottleneck_coefficient_bar.plotly")
    ls.append(output_path + "/report/Reads_Level_Quality/06_fragment_length_line.plotly")
    ls.append(output_path + "/report/Peaks_Level_Quality/01_peak_level_summary_table.dt")
    ls.append(output_path + "/report/Peaks_Level_Quality/02_number_of_peaks_bar.plotly")
    ls.append(output_path + "/report/Peaks_Level_Quality/03_fraction_of_reads_in_peaks_bar.plotly")
    ls.append(output_path + "/report/Peaks_Level_Quality/03_details.yaml")
    ls.append(output_path + "/report/Peaks_Level_Quality/04_peak_annotations_bar.plotly")
    ls.append(output_path + "/report/Peaks_Level_Quality/05_DNAse_I_hypersensitivity_bar.plotly")
    ls.append(output_path + "/report/Downstream/01_details.yaml")

    try:
        if config["motif"] == "homer":
            ls.append(output_path + "/report/Downstream/01_conservation_and_top_motifs.csv")
        else:
            print('not using homer')
    except KeyError:
        print('Homer not used.')
        ls.append(output_path + "/report/Downstream/01_conservation.csv")
        ls.extend(expand(output_path + "/report/Reads_Level_Quality/frag_files/{sample}_frags_hist.csv", sample = list(config["samples"].keys())))
    ls.extend(TRACK_PNG_LIST)
    ls.append(output_path + "/data/pbc_parsed_samplenames.csv")
    ls.append(output_path + "/data/frip_score_parsed_samplenames.csv")
    ls.append(output_path + "/data/putative_targets_parsed_samplenames.csv")
    ls.append(output_path + "/data/dhs_parsed_samplenames.csv")
    ls.append(output_path + "/data/pbc_parsed.csv")
    ls.append(output_path + "/data/dhs_parsed.csv")
    ls.append(output_path + "/report/tracks_all.ini")
    ls.append(output_path + "/report/tracks_all_vlines.ini")
    return ls

###############################################################################

# Copy correct details files to report directory, depending on whether homer module is being used
try:
    if config["motif"] == "homer":
        print("Homer module is being used and details file was copied.")
        rule copy_for_homerconservation_details:
            """Copy the input csv files for homer conservation table"""
            input:
                "cidc_chips/report/not_using_homer_details.yaml" #typo with file naming
            output:
                output_path + "/report/Downstream/01_details.yaml"
            run:
                for i in range(len(input)):
                    shell("cp {infile} {outfile}".format(infile=input[i], outfile= output[i]))
    else:
        print('"Homer module is not being used."')
except KeyError:
    print("Homer module is not being used and relevant details file was copied.")
    rule copy_for_homerconservation_details_nohomer:
        """Copy the input csv files for homer conservation table"""
        input:
             "cidc_chips/report/using_homer_details.yaml"
        output:
            output_path + "/report/Downstream/01_details.yaml"
        run:
            for i in range(len(input)):
                shell("cp {infile} {outfile}".format(infile=input[i], outfile= output[i]))
###############################################################################
rule copy_for_report:
    """Copy the input csv files to report directory so they can be modified"""
    input:
        output_path + "/contam/contamination.csv",
        output_path + "/align/mapping.csv",
        output_path + "/frips/pbc.csv",
        output_path + "/frips/frips.csv",
        output_path + "/peaks/peakStats.csv",
        output_path + "/ceas/meta.csv",
        output_path + "/ceas/dhs.csv",
        "cidc_chips/report/chips_workflow.png",
        "cidc_chips/report/intro_details.yaml",
        "cidc_chips/report/softwares_details.yaml",
        "cidc_chips/report/ref_assembly.csv",
        #"cidc_chips/report/frip_details.yaml",
        "cidc_chips/report/fragsize_details.yaml",
        "cidc_chips/report/trackView_details.yaml"
    output:
        output_path + "/data/contamination2.csv",
        output_path + "/data/mapped_reads.csv",
        output_path + "/data/pbc.csv",
        output_path + "/data/frip_score.csv",
        output_path + "/data/number_of_peaks.csv",
        output_path + "/data/putative_targets.csv",
        output_path + "/data/dhs.csv",
        output_path + "/report/Overview/01_chips_workflow.png",
        output_path + "/report/Overview/01_details.yaml",
        output_path + "/report/Overview/02_details.yaml",
        output_path + "/report/Overview/03_assembly.csv",
        #output_path + "/report/Peaks_Level_Quality/03_details.yaml",
        output_path + "/report/Reads_Level_Quality/06_details.yaml",
        output_path + "/report/Genome_Track_View/01_details.yaml"
    run:
        for i in range(len(input)):
            shell("cp {infile} {outfile}".format(infile=input[i], outfile= output[i]))
###############################################################################
# Write software version information to a text file for user to download (not displayed in the report)
rule software_versions:
    """Store software versions as text file"""
    #input:

    output:
        output_path + "/report/software_versions_all.txt"
    run:
        shell("conda list > {output}")
###############################################################################
# Write select software version information into file
rule software_versions_datatable_noheader:
    """Store software versions in a data table with no header"""
    input:
        output_path + "/report/software_versions_all.txt"
    output:
        temp(output_path + "/report/software_versions_all_noheader.txt")
    run:
        shell("grep -E 'fastp|fastqc|bwa|macs2|bedtools|homer' {input} > {output}")

###############################################################################
# Display select software version information for report table
rule software_versions_datatable:
    """Store software versions in a data table"""
    input:
        output_path + "/report/software_versions_all_noheader.txt"
    output:
        output_path + "/report/Overview/02_select_software_versions.dt"
    run:
        software_nohead = pd.read_table(output_path + "/report/software_versions_all_noheader.txt", header=None, delim_whitespace=True)
        software_nohead.columns =  ['Name', 'Version', 'Build', 'Channel']
        software_nohead.to_csv(output_path + "/report/Overview/02_select_software_versions.dt")
###############################################################################
# Parse contamination for bar graph where myco are summed
rule parse_contamination_for_graph:
    """Parse the contamination output"""
    input:
        output_path + "/data/contamination2.csv"
    output:
        output_path + "/data/contamination2_parsed.csv"
    run:
        df2 = pd.read_csv(output_path + "/data/contamination2.csv", index_col=0)
        myco_cols = [col for col in df2.columns if 'myco_' in col]
        df2['myco'] = df2[myco_cols].sum(axis=1)
        cols_to_include = [col2 for col2 in df2.columns if 'myco_' not in col2]
        df2[cols_to_include].to_csv(output_path + "/data/contamination2_parsed.csv")
###############################################################################
# Parse PBC output
rule parse_pbc:
    """"Parse the PBC output"""
    input:
        output_path + "/data/pbc.csv"
    output:
        output_path + "/data/pbc_parsed.csv"
    run:
        pbc = pd.read_csv(output_path + "/data/pbc.csv")
        pbc['PBC'] = pbc['N1']/pbc['Nd']
        pd.DataFrame(pbc['PBC']).to_csv(output_path + "/data/pbc_parsed.csv")
###############################################################################
# Parse FRIP output
rule parse_frip:
    """Parse FRIP output"""
    input:
        output_path + "/data/frip_score.csv"
    output:
        output_path + "/data/frip_score_parsed.csv"
    run:
        frip = pd.read_csv(output_path + "/data/frip_score.csv")
        pd.DataFrame(frip['FRiP']).to_csv(output_path + "/data/frip_score_parsed.csv")
###############################################################################
# Parse the peak annotation information
rule parse_peak_annotations:
    """Parse peak annotations"""
    input:
        output_path + "/data/putative_targets.csv"
    output:
        output_path + "/data/putative_targets_parsed.csv"
    run:
        peaks_annotated = pd.read_csv(output_path + "/data/putative_targets.csv")
        peaks_annotated['% peaks in promoters'] = 100*(peaks_annotated['Promoter']/peaks_annotated['Total'])
        peaks_annotated['% peaks in exons'] = 100*(peaks_annotated['Exon']/peaks_annotated['Total'])
        peaks_annotated['% peaks in introns'] = 100*(peaks_annotated['Intron']/peaks_annotated['Total'])
        peaks_annotated['% peaks in intergenic regions'] = 100*(peaks_annotated['Intergenic']/peaks_annotated['Total'])
        pd.DataFrame(peaks_annotated[['% peaks in promoters', '% peaks in exons', '% peaks in introns', '% peaks in intergenic regions']]).to_csv(output_path + "/data/putative_targets_parsed.csv")
###############################################################################
# Parse DHS information
rule parse_dhs:
    """Parse DHS"""
    input:
        output_path + "/data/dhs.csv"
    output:
        output_path + "/data/dhs_parsed.csv"
    run:
        dhs = pd.read_csv(output_path + "/data/dhs.csv")
        dhs['% peaks in DHS'] = 100*(dhs['DHS']/dhs['Total'])
        pd.DataFrame(dhs['% peaks in DHS']).to_csv(output_path + "/data/dhs_parsed.csv")
###############################################################################
# Create Downstream analysis table depending on whether the user is running the homer module
try:
    if config["motif"] == "homer":
        print("Homer module is being run.")
        rule make_conservation_table_only_with_homer:
            input:
                HOMERCONSERV_FILES_CONSERV_PROJECT
            output:
                temp(output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_conserv.csv")
            run:
                print(input)
                try:
                    os.mkdir(output_path + "/report/Downstream/plots_conserv/")
                except OSError as error:
                    print(error)
                for i in input:
                    shell( "cp {i} {path}/report/Downstream/plots_conserv/".format(i = i, path = output_path))
                f=open(output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_conserv.csv",'a', encoding='utf-8')
                w=csv.writer(f, delimiter='\t')
                f.write('Sample\tConservation\n')
                for path, dirs, files in os.walk(output_path + "/report/Downstream/plots_conserv/"):
                    print(files)
                    for filename in sorted(files): #naming structure needs to make it sorted by the runRep name
                        print(filename)
                        w.writerow([filename.split("_conserv_thumb")[0], "img:"+path.split("/report/")[1]+"/"+filename.strip()])
                f.close()
    ###############################################################################
        rule make_logo_table_only_with_homer:
            input:
                HOMERCONSERV_FILES_LOGOS_PROJECT
            output:
                temp(output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_logo_noM.csv")
            run:
                print("logo rule running")
                print(input)
                try:
                    os.mkdir(output_path + "/report/Downstream/plots_logo/")
                except OSError as error:
                    print(error)
                for i in input:
                    oldname = i
                    newname = oldname.replace(r'/', '_')
                    shell( "cp {i} {path}/report/Downstream/plots_logo/{newname}".format(i = i, path = output_path, newname = newname))
                f=open(output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_logo_noM.csv",'a', encoding='utf-8')
                w=csv.writer(f, delimiter='\t')
                f.write('Sample\tHomer Motif Logo\n')
                for path, dirs, files in os.walk(output_path + "/report/Downstream/plots_logo/"):
                    print(files)
                    for filename in sorted(files): #naming structure needs to make it sorted by the runRep name
                        print(filename)
                        sample_name_front = filename.split("_motif_")[1]
                        sample_name = sample_name_front.split("_results_")[0]
                        w.writerow([sample_name,"img:"+path.split(output_path + "/report/")[1]+"/"+filename.strip()])
        ###############################################################################
        rule make_homerconservation_table_conserv_noM:
            """Make the homer conservation table conserv names no M"""
            input:
                output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_conserv.csv"
            output:
                temp(output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_conserv_noM.csv")
            run:
                shell("perl -p -e 's/\r//g' {input} > {output}")
        ###############################################################################
        rule homer_conserv_numerical:
            """Homer conservation table numerical pieces"""
            input:
                LOGOS_TXT_FILES
            output:
                temp(output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_numerical.csv")
            run:
                print(input)
                grand_file = open(output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_numerical.csv","a")
                try:
                    with open(sorted(input)[1], "r") as header_line_prep:
                        header_line = header_line_prep.readlines()[0]
                        grand_file.write(header_line)
                    for i in range(len(input)):
                        with open(sorted(input)[i], "r") as sample_file:
                            to_add = sample_file.readlines()[1]
                            grand_file.write(to_add)
                except:
                    with open(str(input), "r") as sample_file:
                        i = 0
                        while i < 2:
                            add_line = sample_file.readline()
                            grand_file.write(add_line)
                            i +=1
                grand_file.close()
        ###############################################################################
        rule homer_conserv_numerical_noM:
            """Make the homer conservation table sample names no M"""
            input:
                output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_numerical.csv"
            output:
                temp(output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_numerical_noM.csv")
            run:
                shell("perl -p -e 's/\r//g' {input} > {output}")
        ###############################################################################
        rule homer_final:
            """Homer conservation table remove certain columns"""
            input:
                output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_numerical_noM.csv",
                output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_conserv_noM.csv",
                output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_logo_noM.csv"
            output:
                output_path + "/report/Downstream/01_conservation_and_top_motifs.csv"
            run:
                df_1 = pd.read_csv(output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_numerical_noM.csv", sep = '\t')
                df_2 = pd.read_csv(output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_conserv_noM.csv", sep = '\t')
                df_3 = pd.read_csv(output_path + "/report/Downstream/01_conservation_and_top_motifs_partial_logo_noM.csv", sep = '\t')
                df_num_conserv = df_1.join(df_2)
                df_num_conserv_logo = df_num_conserv.join(df_3.set_index('Sample'), on='Sample')
                df_num_conserv_logo['Motif'] = df_num_conserv_logo['Motif Name'].str.split('(').str[0]
                df_num_conserv_logo['Negative Log P-value'] = df_num_conserv_logo['Log P-value']*-1
                df_num_conserv_logo[['Sample', 'Conservation', 'Motif', 'Homer Motif Logo', 'Negative Log P-value']].to_csv(output_path + "/report/Downstream/01_conservation_and_top_motifs.csv")
    else:
        print("Not using homer.")
except KeyError:
    # Make Downstream analysis table without motif logo thumbnails
    rule make_conservation_table_only:
        input:
            HOMERCONSERV_FILES_CONSERV_PROJECT
        output:
            output_path + "/report/Downstream/01_conservation.csv"
        run:
            try:
                os.mkdir(output_path + "/report/Downstream/plots_conserv/")
            except OSError as error:
                print(error)
            for i in input:
                shell( "cp {i} {path}/report/Downstream/plots_conserv/".format(i = i, path = output_path))
            f=open(output_path + "/report/Downstream/01_conservation.csv",'a', encoding='utf-8')
            w=csv.writer(f, delimiter='\t')
            f.write('Sample\tConservation\n')
            for path, dirs, files in os.walk(output_path + "/report/Downstream/plots_conserv/"):
                print(files)
                for filename in sorted(files): #naming structure needs to make it sorted by the runRep name
                    print(filename)
                    w.writerow([filename.split("_conserv_thumb")[0], "img:"+path.split("/report/")[1]+"/"+filename.strip()])
###############################################################################
# Parse fragment files for line plot
rule make_lines_for_frag_hist:
    """Make the input csv files for frag plot"""
    input:
        expand(output_path + "/frag/{sample}/{sample}_frags.txt", sample = list(config["samples"].keys()))
    output:
        expand(output_path + "/report/Reads_Level_Quality/frag_files/{sample}_frags_hist.csv", sample = list(config["samples"].keys()))
    run:
        print(input)
        def parse_frag_normed(frag_files_path):
            frag_files = []
            for file in frag_files_path:
                fname_prior = os.path.basename(file)
                fname = os.path.splitext(fname_prior)[0]
                print(fname)
                a = pd.read_csv(file, header=None)
                counts, bins = np.histogram(a, bins=np.linspace(1,1000,200), density=True)
                a_df = pd.DataFrame(counts).reset_index()
                a_df.columns = ['Bin','Density']
                a_df['Fragment size in bp'] = a_df['Bin']*5
                outname = fname+'_hist.csv'
                print(outname)
                outdir = output_path + "/report/Reads_Level_Quality/frag_files/"
                fullname = os.path.join(outdir, outname)
                a_df[['Fragment size in bp', 'Density']].to_csv(fullname)
            return None
        parse_frag_normed(input)
###############################################################################
# Render the lineplot for fragment length
rule frag_plot_paste:
    """Make fragment plots"""
    input:
        expand(output_path + "/report/Reads_Level_Quality/frag_files/{sample}_frags_hist.csv", sample = list(config["samples"].keys()))
    output:
        output_path + "/report/Reads_Level_Quality/06_fragment_length_line.plotly"
    run:
        dfList = []
        for i in range(len(input)):
            dfList.append(pd.read_csv(input[i], index_col=None, header=0).drop('Unnamed: 0', axis=1))
        dfs = [df.set_index('Fragment size in bp') for df in dfList]
        df_save = pd.concat(dfs, axis=1)
        df_save.columns = list(config["samples"].keys())
        df_save.to_csv(output_path + "/report/Reads_Level_Quality/06_fragment_length_line.plotly")
###############################################################################
# Plot contamination
rule cohort_report_data_quality_plots:
    """Plot contamination"""
    input:
        output_path + "/data/contamination2_parsed.csv"
    output:
        csv=output_path + "/report/Reads_Level_Quality/05_contamination_bar.plotly",
        details=output_path + "/report/Reads_Level_Quality/05_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'The reported values for each species represent the percent of 100,000 reads that map to the reference genome of that species.'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'Percentage of 100,000 reads','value':'Percentage of 100,000 reads'}}}),
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details} &&
        cp {input} {output.csv}"""
###############################################################################
# Plot mapping
rule cohort_report_data_quality_plots2:
    """Plot mapping"""
    input:
        output_path + "/data/mapped_reads.csv"
    output:
        csv=output_path + "/report/Reads_Level_Quality/02_mapped_reads_bar.plotly",
        details=output_path + "/report/Reads_Level_Quality/02_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'Mapped reads refer to the number of reads successfully mapping to the genome, while uniquely mapped reads are the subset of mapped reads mapping only to one genomic location.'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'Number of reads','value':'Number of reads'}}}),
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details} &&
        cp {input} {output.csv}"""
###############################################################################
# Create PBC with sample names in csv
rule copy_for_cohort_report_data_quality_plots4:
    """Create PBC with sample names in csv"""
    input:
        output_path + "/data/pbc_parsed.csv",
        output_path + "/data/mapped_reads.csv"
    output:
        temp(output_path + "/data/pbc_parsed_samplenames.csv")
    run:
        df2 = pd.read_csv(output_path + "/data/mapped_reads.csv", index_col=0)
        pd.read_csv(output_path + "/data/pbc_parsed.csv",index_col=0).set_index(df2.index).to_csv(output_path + "/data/pbc_parsed_samplenames.csv")
###############################################################################
# Render PBC
rule cohort_report_data_quality_plots4:
    """Plot PBC"""
    input:
        output_path + "/data/pbc_parsed_samplenames.csv"
    output:
        csv=output_path + "/report/Reads_Level_Quality/03_pcr_bottleneck_coefficient_bar.plotly",
        details=output_path + "/report/Reads_Level_Quality/03_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'The PCR bottleneck coefficient (PBC) refers to the number of locations with exactly one uniquely mapped read divided by the number of unique genomic locations.'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'PBC score','value':'PBC score'}}}),
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details} &&
        cp {input} {output.csv}"""
###############################################################################
# Create FRIP with sample names in csv
rule copy_for_cohort_report_data_quality_plots5:
    """Create FRIP csv with sample names"""
    input:
        output_path + "/data/frip_score_parsed.csv",
        output_path + "/data/mapped_reads.csv"
    output:
        temp(output_path + "/data/frip_score_parsed_samplenames.csv")
    run:
        df2 = pd.read_csv(output_path + "/data/mapped_reads.csv", index_col=0)
        pd.read_csv(output_path + "/data/frip_score_parsed.csv",index_col=0).set_index(df2.index).to_csv(output_path + "/data/frip_score_parsed_samplenames.csv")
###############################################################################
# Render FRIP plot
rule cohort_report_data_quality_plots5:
    """Render FRIP"""
    input:
        output_path + "/data/frip_score_parsed_samplenames.csv"
    output:
        csv=output_path + "/report/Peaks_Level_Quality/03_fraction_of_reads_in_peaks_bar.plotly",
        details=output_path + "/report/Peaks_Level_Quality/03_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'The fraction of reads in peaks (FRIP) score is the fraction of 4 million subsampled reads that fall within a defined peak region.'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'FRIP score (% of reads)','value':'FRIP score (% of reads)'}}}),
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details} &&
        cp {input} {output.csv}"""
###############################################################################
# Render number of peaks
rule cohort_report_data_quality_plots6:
    """Render number of peaks"""
    input:
        output_path + "/data/number_of_peaks.csv"
    output:
        csv=output_path + "/report/Peaks_Level_Quality/02_number_of_peaks_bar.plotly",
        details=output_path + "/report/Peaks_Level_Quality/02_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'The total peaks called, the peaks with a > 10 fold change (10FC), and the peaks with a > 20 fold change (20FC) for each run are represented here.'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'Number of peaks','value':'Number of peaks'}}}),
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details} &&
        cp {input} {output.csv}"""
###############################################################################
# Create csv with sample names for peak annotations
rule copy_for_cohort_report_data_quality_plots7:
    """Copy the csv files for rendering table of read data"""
    input:
        output_path + "/data/putative_targets_parsed.csv",
        output_path + "/data/mapped_reads.csv"
    output:
        temp(output_path + "/data/putative_targets_parsed_samplenames.csv")
    run:
        df2 = pd.read_csv(output_path + "/data/mapped_reads.csv", index_col=0)
        pd.read_csv(output_path + "/data/putative_targets_parsed.csv",index_col=0).set_index(df2.index).to_csv(output_path + "/data/putative_targets_parsed_samplenames.csv")

###############################################################################
# Render peak annotation table
rule cohort_report_data_quality_plots7:
    input:
        output_path + "/data/putative_targets_parsed_samplenames.csv"
    output:
        csv=output_path + "/report/Peaks_Level_Quality/04_peak_annotations_bar.plotly",
        details=output_path + "/report/Peaks_Level_Quality/04_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'The proportions of peaks for each sample overlapping with the promoters, exons, introns, and intergenic regions are shown here.'""",
        plot_options = yaml_dump({'plotly': {'opacity':1.0, 'orientation': "h", 'labels':{'X':'Percentage of peaks','value':'Percentage of peaks'}}}),
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details} &&
        cp {input} {output.csv}"""
###############################################################################
# Create csv for DHS with sample names
rule copy_for_cohort_report_data_quality_plots9:
    input:
        output_path + "/data/dhs_parsed.csv",
        output_path + "/data/mapped_reads.csv"
    output:
        temp(output_path + "/data/dhs_parsed_samplenames.csv")
    run:
        df2 = pd.read_csv(output_path + "/data/mapped_reads.csv", index_col=0)
        pd.read_csv(output_path + "/data/dhs_parsed.csv",index_col=0).set_index(df2.index).to_csv(output_path + "/data/dhs_parsed_samplenames.csv")
###############################################################################
# Render DHS plot
rule cohort_report_data_quality_plots8:
    input:
        output_path + "/data/dhs_parsed_samplenames.csv"
    output:
        csv=output_path + "/report/Peaks_Level_Quality/05_DNAse_I_hypersensitivity_bar.plotly",
        details=output_path + "/report/Peaks_Level_Quality/05_details.yaml",
    params:
        files = lambda wildcards,input: " -f ".join(input),
        caption="""caption: 'DNAse hypersensitive sites (DHS) may represent highly active regions of the genome. The data below represent the percentage of 4 million subsampled peaks that intersect with DHS peaks as defined by list of known DHS regions (specific to each species).'""",
        plot_options = yaml_dump({'plotly': {'barmode':"overlay",'opacity':1.0, 'orientation': "h", 'labels':{'X':'Percentage of peaks','value':'Percentage of peaks'}}}),
    group: "cohort_report"
    shell:
        """echo "{params.caption}" >> {output.details} &&
        echo "{params.plot_options}" >> {output.details} &&
        cp {input} {output.csv}"""
###############################################################################
# Render contamination table without summing myco
rule report_simple_table1:
    input:
         output_path + "/data/contamination2.csv"
    output:
         csv=output_path + "/report/Reads_Level_Quality/04_contamination_table.dt",
         details=output_path + "/report/Reads_Level_Quality/04_details.yaml"
    params:
         caption="caption: 'Contamination percentages for all reference genomes are included here.' "
    shell:
        """echo "{params.caption}" >> {output.details} &&
        cp {input} {output.csv}"""
###############################################################################
# Create raw csv for read level data
rule copy_for_report_table_2:
    """Copy the csv files for rendering table of read data"""
    input:
        output_path + "/report/Reads_Level_Quality/02_mapped_reads_bar.plotly",
        output_path + "/report/Reads_Level_Quality/03_pcr_bottleneck_coefficient_bar.plotly"
    output:
        temp(output_path + "/report/Reads_Level_Quality/01_read_level_summary_table_for_format.csv")
    run:
        mapping_df = pd.read_csv(output_path + "/report/Reads_Level_Quality/02_mapped_reads_bar.plotly", index_col=0)
        pbc_df = pd.read_csv(output_path + "/report/Reads_Level_Quality/03_pcr_bottleneck_coefficient_bar.plotly").set_index(mapping_df.index)
        result2 = pd.concat([mapping_df, pbc_df], axis=1, sort=False)
        result2['Uniquely Mapped'] = result2['UniquelyMapped']
        result2[['Total', 'Mapped', 'Uniquely Mapped', 'PBC']].to_csv(output_path + "/report/Reads_Level_Quality/01_read_level_summary_table_for_format.csv")
###############################################################################
# Format csv for read level data table
rule copy_for_report_table_2_format_data:
    """Format the data in csv files for rendering table of read data"""
    input:
        output_path + "/report/Reads_Level_Quality/01_read_level_summary_table_for_format.csv"
    output:
        temp(output_path + "/report/Reads_Level_Quality/01_read_level_summary_format.csv")
    run:
        read_to_format_df = pd.read_csv(output_path + "/report/Reads_Level_Quality/01_read_level_summary_table_for_format.csv", index_col=0)
        read_to_format_df['Total (M)'] = np.rint(read_to_format_df['Total']/1000000).astype(int)
        read_to_format_df['Mapped (M)'] = np.rint(read_to_format_df['Mapped']/1000000).astype(int)
        read_to_format_df['Uniquely Mapped (M)'] = np.rint(read_to_format_df['Uniquely Mapped']/1000000).astype(int)
        read_to_format_df['PBC'] = np.around(read_to_format_df['PBC'], decimals=2)
        read_to_format_df[['Total (M)', 'Mapped (M)', 'Uniquely Mapped (M)', 'PBC']].to_csv(output_path + "/report/Reads_Level_Quality/01_read_level_summary_format.csv")
###############################################################################
# Render read level data table
rule report_simple_table2:
    input:
         output_path + "/report/Reads_Level_Quality/01_read_level_summary_format.csv"
    output:
         csv=output_path + "/report/Reads_Level_Quality/01_read_level_summary_table.dt",
         details=output_path + "/report/Reads_Level_Quality/01_details.yaml"
    params:
         caption="caption: 'Abbreviations: M, million; PBC, PCR bottlneck coefficient.' "
    shell:
        """echo "{params.caption}" >> {output.details} &&
        cp {input} {output.csv}"""
###############################################################################
# Create csv for peak level data table
rule copy_for_report_table_3:
    """Copy the csv files for rendering table of peak data"""
    input:
        output_path + "/report/Peaks_Level_Quality/03_fraction_of_reads_in_peaks_bar.plotly",
        output_path + "/report/Peaks_Level_Quality/02_number_of_peaks_bar.plotly",
        output_path + "/report/Peaks_Level_Quality/04_peak_annotations_bar.plotly",
        output_path + "/report/Peaks_Level_Quality/05_DNAse_I_hypersensitivity_bar.plotly"
    output:
        temp(output_path + "/report/Peaks_Level_Quality/01_peak_level_summary_table_for_format.csv")
    run:
        df2 = pd.read_csv(output_path + "/report/Peaks_Level_Quality/02_number_of_peaks_bar.plotly", index_col=0)
        df1 = pd.read_csv(output_path + "/report/Peaks_Level_Quality/03_fraction_of_reads_in_peaks_bar.plotly").set_index(df2.index)
        df3 = pd.read_csv(output_path + "/report/Peaks_Level_Quality/04_peak_annotations_bar.plotly").set_index(df2.index)
        df4 = pd.read_csv(output_path + "/report/Peaks_Level_Quality/05_DNAse_I_hypersensitivity_bar.plotly").set_index(df2.index)
        result = pd.concat([df2, df1, df3, df4], axis=1, sort=False)
        result3c = result[['Total', '10FC', '20FC', 'FRiP', '% peaks in promoters', '% peaks in exons', '% peaks in introns', '% peaks in intergenic regions', '% peaks in DHS']].to_csv(output_path + "/report/Peaks_Level_Quality/01_peak_level_summary_table_for_format.csv")
###############################################################################
# Format csv for peak level data table
rule copy_for_report_table_3_format_data:
    """Format the data in csv files for rendering table of peak data"""
    input:
        output_path + "/report/Peaks_Level_Quality/01_peak_level_summary_table_for_format.csv"
    output:
        temp(output_path + "/report/Peaks_Level_Quality/01_peak_level_summary_format.csv")
    run:
        read_to_format_df2 = pd.read_csv(output_path + "/report/Peaks_Level_Quality/01_peak_level_summary_table_for_format.csv", index_col=0)
        read_to_format_df2['% prom'] = np.around(read_to_format_df2['% peaks in promoters'],decimals=2)
        read_to_format_df2['% exons'] = np.around(read_to_format_df2['% peaks in exons'],decimals=2)
        read_to_format_df2['% introns'] = np.around(read_to_format_df2['% peaks in introns'],decimals=2)
        read_to_format_df2['% inter'] = np.around(read_to_format_df2['% peaks in intergenic regions'],decimals=2)
        read_to_format_df2['% DHS'] = np.around(read_to_format_df2['% peaks in DHS'],decimals=2)
        read_to_format_df2[['Total', '10FC', '20FC', 'FRiP', '% prom', '% exons', '% introns', '% inter', '% DHS']].to_csv(output_path + "/report/Peaks_Level_Quality/01_peak_level_summary_format.csv")
###############################################################################
# Render table for peak level data
rule report_simple_table3_format:
    input:
         output_path + "/report/Peaks_Level_Quality/01_peak_level_summary_format.csv"
    output:
         csv=output_path + "/report/Peaks_Level_Quality/01_peak_level_summary_table.dt",
         details=output_path + "/report/Peaks_Level_Quality/01_details.yaml"
    params:
         caption="caption: 'Abbreviations: 10FC, > 10 fold change; 20FC, > 20 fold change; FRiP, Fraction of reads in peaks; Prom, Promoter; Inter, Intergenic; DHS, DNAseI hypersensitivity sites' "
    shell:
        """echo "{params.caption}" >> {output.details} &&
        cp {input} {output.csv}"""
###############################################################################
# Make genome tracks ini file
rule genome_tracks_init:
    """Make genome track initial file"""
    input:
        GENOMETRACK_FILES
    output:
        output_path + "/report/tracks_all.ini"
    run:
        print(input)
        shell("make_tracks_file --trackFiles {input} -o {output}")
###############################################################################
rule genome_track_bed:
    """Make new bed files for genome track"""
    input:
        ref_bed
    output:
        extend_bed,
        tss_bed
    params:
        up = upstream,
        down = downstream
    run:
        pygenomeTacks_bed(ref_bed, extend_bed, tss_bed, upstream, downstream)

# Modify genome tracks ini file to include vlines
rule genome_tracks_tss_init:
    """Modify genome track initial file"""
    input:
        output_path + "/report/tracks_all.ini",
    output:
        output_path + "/report/tracks_all_vlines.ini"
    run:
        with open(output_path + "/report/tracks_all.ini") as f:
            with open(output_path + "/report/tracks_all_vlines.ini", "w") as f1:
                for line in f:
                    f1.write(line)
                f1.write("\n")
                f1.write("[spacer]\n")
                f1.write("[gene]\n")
                f1.write("file={ref_file_to_use}\n".format(ref_file_to_use = extend_bed))
                f1.write("title=Genes\n")
                f1.write("height=2\n")
                f1.write("fontsize=10\n")
                f1.write("file_type=bed\n")
                f1.write("\n")
                f1.write("[vlines]\n")
                f1.write("file={ref_file_to_use}\n".format(ref_file_to_use = tss_bed))
                f1.write("type=vlines\n")

###############################################################################
# Render genome tracks to display
rule genome_tracks_plot:
    """Make genome track plot"""
    input:
        output_path + "/report/tracks_all_vlines.ini",
        extend_bed,
        tss_bed
    output:
        TRACK_PNG_LIST
    run:
        lookup_coords = pd.read_csv(extend_bed,sep = '\t', header=None, index_col=3).iloc[:,-4]
        for list_num, gene in enumerate(GENE_INCLUDED_LIST):
            print(list_num)
            print(gene)
            try:
                region_plot = lookup_coords[gene]
                print(region_plot)
                shell("pyGenomeTracks --tracks {input} --region {region} --trackLabelFraction 0.2 --width 38 --dpi 130 -o {output}".format(input = output_path + "/report/tracks_all_vlines.ini", region = region_plot, output = output[list_num]))
            except:
                print(gene + ' not found')
###############################################################################
rule report_auto_render:
    """Generalized rule to dynamically generate the report BASED
    on what is in the report directory"""
    input:
        report_targets
    params:
        jinja2_template="cidc_chips/report/index.sample.html",
        output_path = output_path + "/report",
        sections_list=",".join(['Overview', "Reads_Level_Quality", "Peaks_Level_Quality", "Genome_Track_View","Downstream"]), #define sections order here
        title="CHIPs Report",
    output:
        output_path+ "/report/report.html"
    message:
        "REPORT: Generating example report"
    shell:
        """python cidc_chips/report.py -d {params.output_path} -s {params.sections_list} -j {params.jinja2_template} -t "{params.title}" -o {output} && cp -r cidc_chips/report/static {params.output_path}"""
