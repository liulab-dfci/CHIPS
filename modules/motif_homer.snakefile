#MODULE: homer motif module--performs motif analysis and generates motif table
import subprocess
# _logfile=output_path + "/logs/motif.log"
_threads=8
_minPeaks = 500

#NOTE: using the _refs from chips.snakefile
def motif_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for run in config["runs"].keys():
        for rep in _reps[run]:
            runRep = "%s.%s" % (run, rep)
            # ls.append(output_path + "/motif/%s/" % runRep)
            ls.append(output_path + "/motif/%s/results/homerResults.html" % runRep)
            ls.append(output_path + "/motif/%s/results/knownResults.txt" % runRep)
            ls.append(output_path + "/motif/%s/results/knownResults/known1.logo.png"% runRep)
            ls.append(output_path + "/peaks/%s/%s_annotatePeaks.txt" % (runRep,runRep))
            ls.append(output_path + "/peaks/%s/%s_annotatePeaks.tsv" % (runRep,runRep))
            ls.append(output_path + "/peaks/%s/%s_annotatePeaks.csv" % (runRep,runRep))
    #ls.append(output_path + "/motif/motifSummary.csv")
    return ls

rule motif_all:
    input:
        motif_targets

def _createEmptyMotif(motif_html):
    """When the _sorted_5k_summits.bed has too few peaks, or is empty,
    we still want to create an emtpy homerResult.html
    INPUT: output paths of these files
    """
    #CHECK for dir existence:
    _path = "/".join(motif_html.split("/")[:-1])
    if not os.path.exists(_path):
        os.makedirs(_path, exist_ok=True)
    #Create an empty mdseqpos_index.html
    subprocess.call(['touch', motif_html])

rule motif_homer:
    """call HOMER on top 5k summits"""
    input:
        bed = output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_5k_summits.bed"
    output:
        # path=directory(output_path + "/motif/{run}.{rep}/"),
        # results=directory(output_path + "/motif/{run}.{rep}/results"),
        html=output_path + "/motif/{run}.{rep}/results/homerResults.html",
        txt=output_path + "/motif/{run}.{rep}/results/knownResults.txt",
        logo=output_path + "/motif/{run}.{rep}/results/knownResults/known1.logo.png"
    params:
        genome=config['motif_path'],
        size=600,
        results=lambda wildcards: output_path + "/motif/%s.%s/results" % (wildcards.run,wildcards.rep)
    message: "MOTIF: calling HOMER on top 5k summits"
    threads:_threads
    log: output_path + "/logs/motif/{run}.{rep}.log"
    benchmark: output_path + "/Benchmark/{run}.{rep}_motif_homer.benchmark"
    #conda: "../envs/motif/motif.yaml"
    run:
        #check to see if _sorted_5k_summits.bed is valid
        wc = str(subprocess.check_output(['wc', '-l', input.bed]))
        #this returns a byte string--we need to eliminate the b' ... '
        #and then convert the first elm to int
        wc = int(wc[2:-1].split()[0])

        if wc >= _minPeaks:
            try:
            #PASS- run motif scan
                shell("findMotifsGenome.pl {input} {params.genome} {params.results} -size {params.size} -p {threads} -mask -seqlogo -preparsedDir {params.results} >>{log} 2>&1")
            except:
                print("homer package not installed")
        else:
            #FAIL - create empty outputs
            _createEmptyMotif(output.html)
            _createEmptyMotif(output.txt)
            _createEmptyMotif(output.logo)

rule motif_getMotifSummary:
    """Summarize the top hits for each run into a file"""
    input:
        #Generalized INPUT fn defined in chips.snakefile
        _getRepInput(output_path + "/motif/$runRep/results/homerResults.html")
    output:
        output_path + "/motif/motifSummary.csv"
    message: "MOTIF: summarizing motif runs"
    # log: output_path + "/logs/motif/{run}.{rep}.log"
    params:
        files = lambda wildcards, input: [" -m %s" % i for i in input]
    conda: "../envs/motif/motif.yaml"
    # run:
    #     files = " -m ".join(input)
    #     shell("cidc_chips/modules/scripts/motif_homerSummary.py -m {files} -o {output} 2>> {log}")
    shell:
        "cidc_chips/modules/scripts/motif_homerSummary.py {params.files} -o {output} " # 2>> {log}"

rule motif_homerAnnotatePeaks:
    """Annotate peak files.
    NOTE: only for motif_homer modules
    """
    input:
        bed = output_path + "/peaks/{run}.{rep}/{run}.{rep}_sorted_peaks.bed"
    output:
        output_path + "/peaks/{run}.{rep}/{run}.{rep}_annotatePeaks.txt"
    params:
        genome=config['motif_path'],
    message: "MOTIF: homer annotatePeaks"
    log: output_path + "/logs/motif/{run}.{rep}.log"
    benchmark: output_path + "/Benchmark/{run}.{rep}_motif_homerAnnotatePeaks.benchmark"
    #conda: "../envs/motif/motif.yaml"
    shell:
        "annotatePeaks.pl {input} {params.genome} > {output}"

rule motif_homerProcessAnnPeaks:
    """Process peaks/{run}/{run}_annotatePeaks.txt files"""
    input:
        output_path + "/peaks/{run}.{rep}/{run}.{rep}_annotatePeaks.txt"
    output:
        tsv=output_path + "/peaks/{run}.{rep}/{run}.{rep}_annotatePeaks.tsv",
        csv=output_path + "/peaks/{run}.{rep}/{run}.{rep}_annotatePeaks.csv",
    message: "MOTIF: Post-process homer annotatePeaks.txt file"
    log: output_path + "/logs/motif/{run}.{rep}.log"
    #conda: "../envs/motif/motif.yaml"
    shell:
        "cidc_chips/modules/scripts/motif_annPeaksTsvCsv.sh {input} {output.tsv} {output.csv}"
