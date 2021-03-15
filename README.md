# Chips- ChIP-seq analysis pipeline in snakemake

# Introduction to Chips
- relation to ChiLin
- relation to viper

# Table of Contents

# Installing Chips
You will only need to install Chips once, either for your own use, or if you are a system administrator, for the entire system (see **Appendix C**).  In other words, you will only need to perform the steps described in this section only once.  
NOTE: this section ends with **Using Chips** (below)

### Required software
We assume that the following tools are already installed on your system and that you have some basic familiarity in using them:
`git`
`wget`
### Installing Miniconda
Chips uses the [Conda](https://conda.io/docs/intro.html) packaging system to manage and install all of its required software packages.
To install miniconda:

1. download the Miniconda installer:

    ```
    $ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    ```

2. run the installer:

    ```
    $ bash Miniconda3-latest-Linux-x86_64.sh
    ```

3. update channels of conda:

    ```
    $ conda config --add channels defaults
    ```

    ```
    $ conda config --add channels bioconda
    ```

    ```
    $ conda config --add channels conda-forge
    ```

### Installing the Chips conda environments
Conda environments are briefly explained [here](https://conda.io/docs/using/envs.html).  Briefly, if you are familiar with [Python Virtual Environments](http://python-guide-pt-br.readthedocs.io/en/latest/dev/virtualenvs/) or [Docker Containers](https://www.docker.com/what-container) then Conda environments should be a familiar concept.  

If you are **not familiar** with these concepts, then a conda environment is simply a **self-contained package space that is composed of various packages.**  So for example, a **bioinformatics** conda space may include packages such as **R**, **samtools**, **bedtools**, etc.

Chips is dependent on conda environments, *chips*.

0. **clone the chips source code**:

    ```
    git clone git@bitbucket.org:plumbers/cidc_chips.git
    ```  
    ** NOTE: this command will create a directory called 'chips'.  After the next five steps, this directory can be safely deleted as we will explain how to *Setup a Chips Project** below. **

1. **installing chips**:  
    After cloning the git repository, create the chips environment by doing this:

    ```
    $ cd cidc_chips  
    $ conda env create -f environment.yml -n chips
    ```

    Or if you have mamba installed in your base environment, a faster method is:  

    ```
    $ mamba env create -f environment.yml -n chips
    ```

    Activate chips Conda Environment:  

    ```
    $ conda activate chips
    ```

2. **Post installation steps: configuring homer**:
    NOTE: Chips uses the [homer](http://homer.ucsd.edu/homer/motif/index.html) software for motif analysis.  It also has the capability of using the [MDSeqPos](https://github.com/XinDong9511/mdseqpos) motif finder for a similar analysis.  If you are interested in using MDSeqPos for motif analysis, please see **Appendix D**.

    To activate/initialize homer:

    *  Run the configure script:

    ```
    $ perl ~/miniconda3/envs/chips/share/homer/.//configureHomer.pl -install
    ```  

    *  Install the required assemblies:

    For human samples:

    ```
    $ perl ~/miniconda3/envs/chips/share/homer/.//configureHomer.pl -install hg38
    ```

    ```
    $ perl ~/miniconda3/envs/chips/share/homer/.//configureHomer.pl -install hg19
    ```

    For mouse samples:

    ```
    $ perl ~/miniconda3/envs/chips/share/homer/.//configureHomer.pl -install mm9
    ```
### Downloading the Chips static reference files
Chips comes pre-packaged with static reference files (e.g. bwa index, refSeq tables, etc.) for hg38/hg19 and mm9/mm10.  You can download those files [ref_files](http://cistrome.org/~galib/ref_files.tar.gz). Many of these files are commonly used static reference files, but if you would like to use the files that you already have, **OR** if you are interested in sup then please see **Appendix E**.

# Using CHIPs
### Anatomy of a Chips project
All work in Chips is done in a **PROJECT/** directory, which is simply a directory to contain a single Chips analysis run.  **PROJECT/** directories can be named anything (and they usually start with a simple mkdir command, e.g. mkdir chips_for_paper),  but what is CRITICAL about a **PROJECT/** directory is that you fill them with the following core components:
(We first lay out the directory structure and explain each element below)
> PROJECT/  
> ├── cidc_chips/  
> ├── data/  - *optional*  
> ├── config.yaml  
> ├── metasheet.csv  
> ├── ref.yaml -  ***only if you are using chips OTHER THAN hg19 and mm9***   
> └── ref_files/

The 'cidc_chips' directory contains all of the chips source code.  We'll explain how to download that directory below.  The 'data' directory is an optional directory that contains all of your raw data. It is optional because those paths __may__ be fully established in the config.yaml, __however__ it is best practice to gather your raw data within 'data' using [symbolic links](https://www.cyberciti.biz/faq/creating-soft-link-or-symbolic-link/).

The *config.yaml* and *metasheet.csv* are configurations for your VIPER run (also explained below).

The ref.yaml file is explained in **Appendix E**.

After a successful **Chips** run, another 'analysis' folder is generated which contains all of the resulting output files.

### Setting up a Chips project
0. **Create Project Directory**
    As explained above, the **PROJECT** directory is simply a directory to contain an entire Chips run.  **It can be named anything, but for this section, we'll simply call it 'PROJECT'**  
    ```
    $ mkdir PROJECT
    ```

    ```
    $ cd PROJECT
    ```

1. **Create Data Directory**
    As explained above, creating a data directory is a place to gather all of your **raw data files (.fastq, .fastq.gz, .bam)**.  It is optional, but **highly recommended**.
    ```
    $ mkdir data
    ```
    And in 'data', copy over or make symbolic links to your raw data files  

2. **Clone CHIPs Repository**
    In your PROJECT directory:  
    ```
    $ mv cidc_chips/ PROJECT/
    ```

3. **Create config.yaml and metasheet.csv**

    a. **copy chips/config.yaml and chips/metasheet.csv into the PROJECT dir:**

    In the PROJECT directory:

    ```
    $ cp cidc_chips/config.yaml .
    ```

    ```
    $ cp cidc_chips/metasheet.csv .
    ```

    b. **setup config.yaml**
        The config.yaml is where you define Chips run parameters and the ChIP-seq samples for analysis.

    * **genes_to_plot**: If set, genomic region and TSS will be displayed in Genome Trackview figure. Multiple genes should be separated by space (default: GAPDH ACTB TP53).
    * **upstream/downstream**: Upstream and Downstream of the genome region can be extended to have a better view of peaks.
    * **output_path**: Directory to save all the output files (default: analysis).
    * **assembly**: typically hg19/hg38 for human or mm9/mm10 for mouse (default: hg19)
    * **Choose the motif software**: choose either homer or MDSeqPos (default: homer)
    * **Contamination Panel**:The contamination panel is a panel that Chips will check for "cross-species" contamination. Out of the box, the config.yaml has hg19 and mm9 as assemblies to check.  **IF you would like to add other species/assemblies, simply add as many BWA indices as you would like**
    * **cnv_analysis**: Set to 'true' to enable copy number variation analysis
    * **samples**: __The most important part of the config file is to define the samples for Chips analysis.__ Each sample is given an arbitrary name, e.g. MCF7_ER, MCF7_input, etc.  **Sample names, however, can not start with a number, and cannot contain '.', '-' (dash--use underscores instead)** (POSSIBLY others). For each sample, define the path to the raw data file (.fastq, .fastq.gz, .bam). For paired-end samples, simply add another line to define the path to the second pair.

    c. **setup metasheet.csv**:
    The metasheet.csv is where you group the **samples** (defined in config.yaml) into Treatment, Control (and if applicable, replicates).  For Chips, each of these groupings is called a **run**.  

    Open metasheet.csv in Excel or in a text-editor.You will notice the first (uncommented) line is:

    `RunName,Treat1,Cont1,Treat2,Cont2`

    **RunName**- arbitrary name for the run, e.g. *MCF7_ER_run*  
    **Treat1**- The sample name that corresponds to treatment sample.  **It must exactly match the sample name used in config.yaml**  
    **Cont1**- (optional) The input-control sample that should be paired with Treat1.  
    **Treat2**- (optional) if you have replicates, define the treatment of the replicate here.  
    **Cont2**- (optional) the input-control, if available, for Treat2  

4. **Set Up Refs**
    - A pre-built [ref_files](http://cistrome.org/~galib/ref_files.tar.gz) can be downloaded from the link.
    - makesure in config.yaml, ref: "cidc_chips/ref.yaml"
    - linking to static refs.
    - copying ref.yaml

### Running Chips

1. Acitivate the environment

```bash
conda activate chips
```

2. dry run

```
$ snakemake -np  -s cidc_chips/chips.snakefile --rerun-incomplete
```

3. full run  
```
$ nohup snakemake -s cidc_chips/chips.snakefile --rerun-incomplete -j 8 > run.out &
```

More information for using snakemake can be found [here](https://snakemake.readthedocs.io/en/stable/index.html).

### Appendix A: System requirements
### Appendix B: Recommended requirements
### Appendix C: Installing Chips system-wide
###### for system administrator, those who wish to share their Chips installation
### Appendix D: Installing the MDSeqPos motif finder for chips

```
$ conda activate chips
$ cd mdseqpos/lib
$ cp settings.py.example settings.py
```
Modify `settings.py` like below:
```python
#This should be absolute directory where you ref_files folder is.
ASSEMBLY_DIR = '***/***/ref_files'
BUILD_DICT = { "hg19": "hg19/",
               "hg38": "hg38/",
               "mm9":"mm9/",
               "mm10": "mm10/"
               }
```
Then do:

```bash
cd ..
./version_updater.py
python setup.py install
```
At last, type `MDSeqPos.py` to ensure MDSeqPos is installed and check the usage.

### Appendix E: Generating static reference files for Chips
- all of the required files  
- using your own files
- supporting something new
- adding to ref.yaml
