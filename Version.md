# Versions of Tools and Reference Files Used in Chips

## Tools

### environment: chips(stable)

| Software         | Version | Source                | Notes |
|:-----------------|:--------|:----------------------|:------|
| snakemake        | 5.4.5   | bioconda              |       |
| samtools         | 1.10    | bioconda              |       |
| python           | 3.6.12  | conda-forge           |       |
| r                | 3.5.1   | conda-forge           |       |
| numpy            | 1.17.3  | conda                 |       |
| bwa              | 0.7.15  | bioconda              |       |
| bowtie2          | 2.3.4.1 | bioconda              |       |
| picard           | 2.20.0  | bioconda              |       |
| bedtools         | 2.27.1  | bioconda              |       |
| seqtk            | 1.3     | bioconda              |       |
| fastqc           | 0.11.9  | bioconda              |       |
| fastp            | 0.20.1  | bioconda              |       |
| ggplot2          | 3.3.0   | conda-forge r         |       |
| reshape2         | 1.4.4   | conda-forge r         |       |
| git              | 2.26.0  | conda-forge           |       |
| perl             | 5.26.2  | conda-forge           |       |
| homer            | 4.11    | bioconda              |       |
| weblogo          | 2.8.2   | bioconda              |       |
| seqLogo          | 1.50.0  | bioconda bioconductor |       |
| bedgraphtobigwig | 377     | bioconda ucsc         |       |
| bedsort          | 377     | bioconda ucsc         |       |
| qdnaseq          | 1.18.0  | bioconda bioconductor |       |
| seaborn          | 0.11.1  | conda-forge           |       |
| r.utils          | 2.9.2   | conda-forge r         |       |
| pybigwig         | 0.3.17  | bioconda              |       |
| pybedtools       | 0.8.1   | bioconda              |       |
| numpy            | 1.19.5  | conda                 |       |
| cython           | 0.29.21 | conda                 |       |
| jinja2           | 2.11.2  | conda                 |       |
| macs2            | 2.2.7.1 | bioconda              |       |


## Reference

All reference file could be downloaded from [here](http://cistrome.org/~xindong/chips_reference_files/).

### CIDC

##### GDC_hg38

| Reference Key               | Location                                          | Version        | Source                                                                                  | Notes                           |
|:----------------------------|:--------------------------------------------------|:---------------|:----------------------------------------------------------------------------------------|:--------------------------------|
| bwa_index                   | ./ref_files/GDC_hg38/bwa_indices/GRCh38.d1.vd1.fa | GRCh38         | https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files | GDC.h38.d1.vd1 BWA Index Files  |
| geneTable                   | ./ref_files/GDC_hg38/GDC_hg38.refGene             | GENCODE v22    | Download from UCSC genome table browser                                                 | used for calculate ceas and CNV |
| geneBed                     | ./ref_files/hg38/GDC_hg38.bed                     | GENCODE v22    | reformat gtf by cidc_chips/static/scripts/GtfToBed.py                                   | used for calculating RP         |
| conservation                | ./ref_files/GDC_hg38/conservation/                | GRCh38, 100way | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/                        |                                 |
| DHS                         | ./ref_files/GDC_hg38/regions/hg38_cCREs.bed       |                | cCRE regions from ENCODE Project                                                        |                                 |
| exons                       | ./ref_files/GDC_hg38/regions/exon.bed             | GENCODE v22    | extract from geneTable                                                                  |                                 |
| promoters                   | ./ref_files/GDC_hg38/regions/promoter.bed         | GENCODE v22    | extract from geneTable                                                                  |                                 |
| velcro_regions              | MISSING                                           | -              | -                                                                                       | Blacklist Region                |
| chrom_lens                  | ./ref_files/GDC_hg38/regions/chromInfo_hg38.txt   |                | extract from rawgenome by `samtools faidx`                                              |                                 |
| rawgenome (not in ref.yaml) | ./ref_files/GDC_hg38/rawgenome                    |                | split rawgenome into each chromosome                                                    | required by MDSeqPos            |
| masked (not in ref.yaml)    | ./ref_files/GDC_hg38/masked                       |                |                                                                                         | required by MDSeqPos            |

### Cistrome

##### hg38

| Reference Key  | Location                                              | Version      | Source                                                                           | MD5 | Notes                              |
|:---------------|:------------------------------------------------------|:-------------|:---------------------------------------------------------------------------------|:----|:-----------------------------------|
| bwa_index      | ./ref_files/hg38/bwa_indices/hg38/hg38.fa             | hg38         | http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz       |     | use bwa to build index             |
| geneTable      | ./ref_files/hg38/hg38.refGene                         | refseq hg38  | refseq from UCSC table browser                                                   |     | Calculate RP                       |
| geneBed        | ./ref_files/hg38/hg38_refGene.bed                     | refseq hg38  | reformat feature table by cidc_chips/static/scripts/FeatureTableToBed.py         |     |                                    |
| conservation   | ./ref_files/hg38/conservation/hg38.phastCons100way.bw | hg38, 100way | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/                 |     | Previous using 7way in Cistrome DB |
| DHS            | ./ref_files/hg38/regions/DHS_hg38.bed                 |              | Union DHS regions from Cistrome DB                                               |     |                                    |
| exons          | ./ref_files/hg38/regions/exon.bed                     |              | extract from geneTable                                                           |     |                                    |
| promoters      | ./ref_files/hg38/regions/promoter.bed                 |              | extract from geneTable                                                           |     |                                    |
| velcro_regions | MISSING                                               | -            | -                                                                                |     | Blacklist Region                   |
| chrom_lens     | ./ref_files/hg38/regions/chromInfo_hg38.txt           |              | UCSC table browser                                                               |     |                                    |
| rawgenome      | ./ref_files/hg38/rawgenome/                           |              | http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz       |     | required by MDSeqPos               |
| masked         | ./ref_files/hg38/masked/                              |              | http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFaMasked.tar.gz |     | required by MDSeqPos               |


##### mm10

| Reference Key  | Location                                              | Version     | Source                                                                      | MD5 | Notes                                               |
|:---------------|:------------------------------------------------------|:------------|:----------------------------------------------------------------------------|:----|:----------------------------------------------------|
| bwa_index      | ./ref_files/mm10/bwa_indices/mm10.fa                  | mm10        | http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/chromFa.tar.gz       |     | use bwa to build index                              |
| geneTable      | ./ref_files/mm10/mm10.refGene                         | refseq mm10 | refseq from UCSC table browser                                              |     | Calculate RP                                        |
| geneBed        | ./ref_files/mm10/mm10_refGene.be                      | refseq mm10 | reformat feature table by cidc_chips/static/scripts/FeatureTableToBed.py    |     |                                                     |
| conservation   | ./ref_files/mm10/conservation/mm10.60way.phastCons.bw | mm10, 60way |                                                                             |     | Previous using 60way in Cistrome DB                 |
| DHS            | ./ref_files/mm10/regions/mm10.DHS.bed                 |             | Union DHS regions from Cistrome DB                                          |     | merging all the peaks of DNase-seq data from ENCODE |
| exons          | ./ref_files/mm10/regions/mm10.exon.bed                |             | extract from geneTable                                                      |     |                                                     |
| promoters      | ./ref_files/mm10/regions/mm10.promoter.bed            |             | extract from geneTable                                                      |     |                                                     |
| velcro_regions | MISSING                                               | -           | -                                                                           |     | Blacklist Region                                    |
| chrom_lens     | ./ref_files/mm10/regions/mm10.len                     |             | UCSC table browser                                                          |     |                                                     |
| rawgenome      | ./ref_files/hg38/rawgenome/                           |             | http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/chromFa.tar.gz       |     | required by MDSeqPos                                |
| masked         | ./ref_files/hg38/masked/                              |             | http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/chromFaMasked.tar.gz |     | required by MDSeqPos                                |

#### Contamination Panel

| Reference Key  | Location                                                                          | Version | Source | MD5 | Notes |
|:---------------|:----------------------------------------------------------------------------------|:--------|:-------|:----|:------|
| hg19           | ./ref_files/contam_panel/hg19/hg19.fa                                             |         |        |     |       |
| mm9            | ./ref_files/contam_panel/mm9/mm9.fa                                               |         |        |     |       |
| dm3            | ./ref_files/contam_panel/dm3/dm3.fa                                               |         |        |     |       |
| S_cerevisiae   | ./ref_files/contam_panel/S_cerevisiae/S_cerevisiae.fa                             |         |        |     |       |
| e_coli         | ./ref_files/contam_panel/e_coli/e_coli.fasta                                      |         |        |     |       |
| myco_PG-8A     | ./ref_files/contam_panel/mycoplasma/GCF_000018785.1_ASM1878v1/myco_PG-8A.fna      |         |        |     |       |
| myco_ATCC23114 | ./ref_files/contam_panel/mycoplasma/GCF_000085865.1_ASM8586v1/myco_ATCC23114.fna  |         |        |     |       |
| myco_m64       | ./ref_files/contam_panel/mycoplasma/GCF_000186005.1_ASM18600v1/myco_m64.fna       |         |        |     |       |
| myco_SK76      | ./ref_files/contam_panel/mycoplasma/GCF_000313635.1_ASM31363v1/myco_SK76.fna      |         |        |     |       |
| myco_ATCC23714 | ./ref_files/contam_panel/mycoplasma/GCF_000420105.1_ASM42010v1/myco_ATCC23714.fna |         |        |     |       |
| myco_ATCC23064 | ./ref_files/contam_panel/mycoplasma/GCF_000485555.1_ASM48555v1/myco_ATCC23064.fna |         |        |     |       |
| myco_HAZ145_1  | ./ref_files/contam_panel/mycoplasma/GCF_001547975.1_ASM154797v1/myco_HAZ145_1.fna |         |        |     |       |
| myco_ATCC29342 | ./ref_files/contam_panel/mycoplasma/GCF_000027345.1_ASM2734v1/myco_ATCC29342.fna  |         |        |     |       |
| myco_R         | ./ref_files/contam_panel/mycoplasma/GCF_000092585.1_ASM9258v1/myco_R.fna          |         |        |     |       |
| myco_WVU1853   | ./ref_files/contam_panel/mycoplasma/GCF_000969765.1_ASM96976v1/myco_WVU1853.fna   |         |        |     |       |
