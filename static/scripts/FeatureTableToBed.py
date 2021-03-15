#!/usr/bin/env python
# -*- coding: utf-8 -*-


import pandas as pd
import argparse

def add_tss(bed_table):
    bed_table.loc[bed_table.strand == '+', 'TSS'] = bed_table.loc[:,"start"]
    bed_table.loc[bed_table.strand == '-', 'TSS'] = bed_table.loc[:,"end"]
    bed_table["TSS"] = bed_table.loc[:,"TSS"].astype(int)
    return bed_table


def change_format(feature):
    feature_bed = feature.loc[:,['chromosome','start','end','product_accession','strand','symbol']]
    feature_bed.loc[:,'chromosome'] = ['chr'+str(x) for x in feature_bed.loc[:,'chromosome'].tolist()]
    # gtf is 1-based system, but bed is 0-based
    feature_bed["start"] = feature_bed.loc[:,"start"] - 1
    feature_bed["end"] = feature_bed.loc[:,"end"] - 1
    feature_bed = add_tss(feature_bed)
    feature_bed['coordinate'] = [x[0]+':'+str(x[1])+'-'+str(x[2]) for x in feature_bed.values.tolist()]
    feature_bed = feature_bed.loc[:,['chromosome','start','end','coordinate','product_accession','strand','symbol','TSS']]
    return feature_bed

def getnew(refGene_path,feature_path,full_table_path):
    #filter:  mrna chromosome protein_coding
    #Read in annotation file
    feature = pd.read_csv(feature_path,sep='\t',low_memory=False) 
    feature = feature[feature['seq_type'] == 'chromosome']
    protein_codings = feature[feature['class']=='protein_coding'].loc[:,'symbol'].tolist()
    new_feature = feature[(feature['# feature']=='mRNA') & (feature['symbol'].isin(protein_codings))]
    new_feature_bed = change_format(new_feature)
    if refGene_path:
        refGene = pd.read_csv(refGene_path,sep = '\t',names = ['bin_num','product_accession','chromosome','strand','start','end','cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','score','symbol','cdsStartStat','cdsEndStat','exonFrames'], low_memory=False)
        refGene = refGene[refGene['chromosome'].isin(['chr'+str(i) for i in range(1, 23)]+['chrX', 'chrY'])]
        feature = feature[feature['# feature']=='gene']
        all_feature_bed = change_format(feature)
        # get the symbol in refGene but not new feature
        refgene_symbol = list(set(refGene.loc[:,"symbol"].tolist()) - set(new_feature_bed.loc[:,"symbol"].tolist()))
        gene_ann = all_feature_bed[all_feature_bed.loc[:,"symbol"].isin(refgene_symbol)]
        # Final full table
        full_table = pd.concat([new_feature_bed, gene_ann],axis=0,sort=True)
    else:
        full_table = new_feature_bed
    # full_table=add_tss(full_table)
    full_table = full_table.loc[:,['chromosome','start','end','coordinate','product_accession','strand','symbol','TSS']]
    full_table.to_csv(full_table_path, index = None, sep = '\t')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare annotation for RP calculation')
    parser.add_argument("-r", "--refGene", dest="refGene", default= "" , help="refGene file that downloaded from USCS table")
    parser.add_argument("-f", "--feature", dest="feature", required=True, help="Feature table for updating your annotation")    
    parser.add_argument("-o", "--output", dest="output", required=True, help='prefix for the output file')
    args = parser.parse_args()

    refGene_path=args.refGene
    feature_path=args.feature
    full_table_path=args.output

    getnew(refGene_path,feature_path,full_table_path)




