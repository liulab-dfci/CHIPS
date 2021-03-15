#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pyranges as pr
import pandas as pd
import argparse

def gtfToBed(gtf,output):
    gr = pr.read_gtf(gtf)
    df = gr.df
    geneAndTrans = df[(df["Feature"] == "gene" )|(df["Feature"] == "transcript")]
    AnnoBed = geneAndTrans.loc[:,["Chromosome","Start","End","gene_name","Strand","gene_id"]]
    AnnoBed = AnnoBed.drop_duplicates()
    AnnoBed = AnnoBed.rename(columns={'Chromosome': 'chromosome', 'Start': 'start', 'End':'end' , 'gene_name':'symbol', 
                                      'Strand':'strand', 'gene_id':'product_accession'})
    AnnoBed.loc[:,"start"] = AnnoBed.loc[:,"start"] - 1
    AnnoBed.loc[:,"end"] = AnnoBed.loc[:,"end"] - 1
    AnnoBed.loc[AnnoBed.strand == '+', 'TSS'] = AnnoBed.start
    AnnoBed.loc[AnnoBed.strand == '-', 'TSS'] = AnnoBed.end
    AnnoBed.TSS = AnnoBed.TSS.astype(int)
    AnnoBed["coordinate"] = [x[0]+':'+str(x[1])+'-'+str(x[2]) for x in AnnoBed.values.tolist()]
    AnnoBed = AnnoBed[['chromosome','start','end','coordinate','product_accession','strand','symbol','TSS']]
    AnnoBed.to_csv(output, index = None, sep = '\t')

def main():
    parser = argparse.ArgumentParser(description='Prepare annotation for RP calculation')
    parser.add_argument("-g", "--gtf", dest="gtf", required=True, help="gtf table")    
    parser.add_argument("-o", "--output", dest="output", required=True, help='output file')
    args = parser.parse_args()

    gtf=args.gtf
    output=args.output
    gtfToBed(gtf,output)
    
if __name__ == "__main__":
    main()
