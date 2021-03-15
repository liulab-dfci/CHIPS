import sys
import pandas as pd
import argparse
import time
import re
import os
import math
import collections
    
def Info(infoStr):
    print("[%s] %s" % (time.strftime('%H:%M:%S'), infoStr))    

def readPeaks(file_path):
    Info("Start reading peaks bed file.")
    peaks = []  # chr summit peak None
    peaki = 0

    with open(file_path, "r") as pfile:
        for line in pfile:
            if line.startswith("#"):
                continue
            else:
                try:
                    peaks.append((line.rstrip().split("\t")[0],
                              (int(line.rstrip().split("\t")[1]) + int(line.rstrip().split("\t")[2])) / 2.0, 0,
                              "peak_%s" % peaki))
                    peaki += 1
                except:
                    continue
    peaks_top = peaks[:]
    Info("Finish reading peaks bed file. %s peaks in total. " % peaki)
    return(peaks, peaks_top)


def readCoordinate(new_ref_path):
    genome = new_ref_path
    coordinate = []
    gene_ann = collections.defaultdict()
    with open(genome, "r") as gfile:
        next(gfile)
        for line in gfile.readlines():
            line = line.rstrip().split('\t')
            if line[4] != '':
                # append chromosome, TSS, 1, refSeq id
                coordinate.append((line[0], int(line[7]), 1, line[4]))
                gene_ann[line[4]] = [line[0], line[1], line[2], line[4], line[5], line[6]]
    return(coordinate, gene_ann)

def propCalc(genes_list, peaks_list):
    PEAK = 0
    GENE = 1
    prop_dict = {}
    for top_peak_number in [2000, 5000, 10000, 1000000]:
        w = genes_list + peaks_list[:top_peak_number]
        D = {}
        A = {}
        w.sort()
        for elem in w:
            if elem[2] == PEAK:
                A[elem[-1]] = [elem[0], elem[1]]
                D[elem[-1]] = float("inf")
            else:
                for peak in list(A.keys()):
                    p_elem = A[peak]
                    if p_elem[0] == elem[0]:
                        D[peak] = elem[1] - p_elem[1]
                        A.pop(peak)
                    else:
                        A.pop(peak)                  
        w.reverse()
        for elem in w:
            if elem[2] == PEAK:
                A[elem[-1]] = [elem[0], elem[1]]
            else:
                for peak in list(A.keys()):
                    p_elem = A[peak]
                    if p_elem[0] == elem[0]:
                        D[peak] = min(p_elem[1] - elem[1], D[peak])
                        A.pop(peak)
                    else:
                        A.pop(peak)
        D_values= D.values()
        D_len = float(len(D_values))
        peak_within_1k = [i for i in D.values() if i <=1000]
        prop_dict[top_peak_number] = len(peak_within_1k) / D_len
        
    return(prop_dict)


# Get peaks with the distance (15*decay) of gene's TSS
def peaksInRange(genes_list, peaks_list, decay):
    # genes_list [chr, tss, 1, name]
    # peaks_list [chr, summit/center, 0, None]
    PEAK = 0
    GENE = 1
    padding = decay * 15
    
    w = genes_list + peaks_list

    D = {}
    A = {}
    D2 = {}
    A2 = {}

    w.sort()
    for elem in w:
        if elem[2] == GENE:
            A[elem[-1]] = [elem[0], elem[1]]
            D[elem[-1]] = []
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                if (g[0] != elem[0]) or ((elem[1] - g[1]) > padding):
                    dlist.append(gene_name)
                else:
                    A[gene_name].append(elem[1] - g[1])  # peak in distance will calculate the distance
            for gene_name in dlist:
                D[gene_name] += A.pop(gene_name)[2:]

    w.reverse()
    for elem in w:
        if elem[2] == GENE:
            A[elem[-1]] = [elem[0], elem[1]]
        else:
            dlist = []
            for gene_name in list(A.keys()):
                g = A[gene_name]
                if (g[0] != elem[0]) or (-(elem[1] - g[1]) > padding):
                    dlist.append(gene_name)
                else:
                    A[gene_name].append(-(elem[1] - g[1]))
            for gene_name in dlist:
                    D[gene_name] += A.pop(gene_name)[2:]

    for gene_name in list(A.keys()):
        D[gene_name] += A.pop(gene_name)[2:]


    return(D)

def calcRPscore(peaks_file, genes_file, decay):
    peaks_list, peaks_list_top = readPeaks(peaks_file)
    genes_list, gene_ann = readCoordinate(genes_file) # gene_ann for retreve gene inforamtion
    if decay == "auto":
        prop_dict = propCalc(genes_list, peaks_list)
        if max(prop_dict.values()) >= 0.2:
            Info(str(max(prop_dict.values())) + """ of peaks are located in the promotor region.
                 This is more than the 20% promoter-type threshold so the the decay distance is set to 1.0kb,
                 appropriate for promotor-type analysis. The half decay distance can be specified in the parameters.""")
            decay = 1000
        else:
            Info(str(max(prop_dict.values())) + """ of peaks are located in the promotor region.
                 This is less than the 20% promoter-type threshold so the the decay distance is set to 10.0kb,
                 appropriate for enhancer-type analysis. The half decay distance can be specified in the parameters.""")
            decay = 10000
    else:
        decay = float(decay) #float(str(decay)[:-1]) * 1000  #remove k
    # D, peaks_in_range = peaksInRange(genes_list, peaks_list_top, decay)
    D = peaksInRange(genes_list, peaks_list_top, decay)
    RP = []
    for gene_name in D: # gene_name here is refseq id
        score = sum(map(lambda x: 2.0 ** x, map(lambda x: -x / decay, D[gene_name])))
        ann = gene_ann[gene_name] # retreve gene annotation

        # output format: chr, start, end, refseq, score, strand, symbol
        entry = ann[0:4] + [score] + ann[4:] 
        RP.append(entry)
    # RP_df = pd.DataFrame(list(RP.items()))
    return(RP, decay)


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return(arg)  # return an open file handle


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate RP score.')
    parser.add_argument("-p", "--peak",dest="peak_file", required=True,
                    help="input peak file", metavar="FILE",
                    type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-a", "--annotation",dest="annotation_file", required=True,
                    help="genome annotation file with format: chromosome, start, end, label, product_accession, strand, symbol, TSS ", metavar="FILE",
                    type=lambda x: is_valid_file(parser, x))    
    parser.add_argument("-d", dest="decay", type=int,
                        help="Half decay distance to use. Select from auto, or a integer (base pair). By default is automatically use 10000 or 1000")
    parser.add_argument("-n", dest="prefix", type=str,
        help='prefix for the output file')
    
    args = parser.parse_args()
    peaks_file = args.peak_file
    genome = args.annotation_file
    decay = 'auto' if not args.decay else args.decay
    output = args.prefix

    #====== run rp caculator =====
    RP_result, decay_used = calcRPscore(peaks_file, genome, decay) # return a list
    RP_result.sort(key=lambda x:x[-3], reverse=True)

    #======= write out rp ========
    if not output:
        output = 'gene_score.txt'

    out = open(output, 'w')

    
    output1 = output.split('/')[-1]

    opts_string = "# Argument List:\n" +\
                           "# Name = %s\n" % output1 +\
                           "# peak file = %s\n" % peaks_file +\
                           "# decay = %d bp\n" % decay_used +\
                           "# genome = %s\n" % genome
    out.write(opts_string)
    out.write('#chrom\ttxStart\ttxEnd\trefseq\tscore\tstrand\tsymbol\n')
    for line in RP_result:
        out.write('%s\t%s\t%s\t%s\t%.3f\t%s\t%s\n'%(
                       line[0], line[1], line[2], line[3], line[4], line[5], line[6]))
    out.close()
    Info("Finished! result output to <%s>"%output)


