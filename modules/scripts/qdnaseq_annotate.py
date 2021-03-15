#!/usr/bin/env python
"""
Copyright 2014 Len Taing

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------------------
UPDATED 2017-07-19: ported to python3 AND used to process segmented.igv, a 
qdnaseq output file

Outputs: qdnaseq_genes.igv
-------------------------------------------------------------------------------
A program to classify genomic sites into one of four categories:
promoter, exon, intron, intergenic.

Input: a bed file of genomic sites (note: the conversion of bed regions to
sites should be done by the user), we use chr and start (i.e. col0 and col1)

refGenes geneTable- expected format is: see UCSC hg19 refGenes table as example
bin name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts 
exonEnds score name2 cdsStartStat cdsEndStat exonFrames

Output:
{Promoter: %sites in promoter,
 Exon: %sites in exon,
 Intron: %sites in intron,
 Intergenic: %sites in intergenic
}

NOTE: intergenic is the 'other' category--meaning if it's neither promoter,
exon, or intron, then the default classification is "intergenic".  
THEREFORE: the percents will = 100%

UPDATE: will also output 4 files, XXX_promoter.bed, XXX_exon.bed, 
XXX_intron.bed, XXX_intergenic.bed
"""
import os
import sys
from optparse import OptionParser

import bisect
import math
import copy

def inRegion(reg, site):
    """given a region tuple, (start, end) 
    returns True if site is >= start <= end"""
    return site >= reg[0] and site <= reg[1]

class Gene():
    """Gene object--will help us store gene information and perform some fns"""
    #should make these command line parameters!
    upstream = 0 #gene starts at TSS
    downstream = 0 #gene ends at TTS

    def __init__(self, line):
        """Given a line/row from the UCSC geneTable, grab: 
        name, name2, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exons
        NOTE: refGenes locations are oriented according to the 5' orientation
        ensuring that *Start will always be < *End even though the semantics
        are reversed for - strand genes. 
        **WE will keep this 5' orientation**
        """
        tmp = line.strip().split("\t")
        self.name = tmp[1] #UCSC ID
        self.name2 = tmp[-4] #human readable gene name
        self.chrom = tmp[2]
        self.strand = tmp[3] #can only be "+" or "-"
        self.txStart = int(tmp[4])
        self.txEnd = int(tmp[5])
        self.cdsStart = int(tmp[6])
        self.cdsEnd = int(tmp[7])
        exStarts = [int(i) for i in tmp[9].split(",") if i]
        exEnds = [int(i) for i in tmp[10].split(",") if i]
        self.exons = zip(exStarts, exEnds)

    def chunk(self):
        """a method that will return the start and end sites of the gene,
        e.g. the 'chunk' of chromosome it occupies.
        NOTE: accounts for promoter region, which uses the upstream!
        returns (start, end)

        this is used to see if a site falls within a gene 'chunk'
        """
        if self.strand == '+':
            #check for negatives?
            return (self.txStart - Gene.upstream, self.txEnd + Gene.downstream)
        else: #negative strand
            #check for negative?
            return (self.txStart - Gene.downstream, self.txEnd + Gene.upstream)
    
    def siteInGene(self, site):
        """Given a site, e.g. 5246835, returns True if the site falls within
        the gene chunk, else False
        """
        c = self.chunk()
        #return site >= c[0] and site <= c[1]
        return inRegion(c, site)

    def __str__(self):        
        ls = ["Gene: %s\%s" % (self.name,self.name2),
              "Loc: %s:%s-%s" % (self.chrom, self.txStart, self.txEnd)]
        return "\n".join(ls)

    #Comparisons based on genomic location--assumming other is on same chr
    def __eq__(self, other):
        return self.chunk() == other.chunk()
    def __ne__(self, other):
        return self.chunk() != other.chunk()  ## discuss with len 201486
    #order by chunk()[0]s
    def __lt__(self, other):
        return self.chunk()[0] < other.chunk()[0]
    def __le__(self, other):
        return self.chunk()[0] <= other.chunk()[0]
    def __gt__(self, other):
        """order by chunk()[0]s-Return true if self is DOWNSTREAM from other"""
        return self.chunk()[0] > other.chunk()[0]
    def __ge__(self, other):
        return self.chunk()[0] >= other.chunk()[0]

class Chromosome(list):
    """list of SORTED Gene objects (sorted by gene.chunk()[0])"""
    #gene ids makes sure we only take the 1st isoform of the gene we encounter
    #NO redundancies!
    geneIds = None
    
    def __init__(self, *args):
        list.__init__(self, *args)
        self.geneIds = {}

    def add(self, gene):
        """adds the gene into the sorted list"""
        if gene.name not in self.geneIds:
            #NEW gene to add
            self.geneIds[gene.name] = gene.name2
            i = bisect.bisect(self, gene)
            self[i:i] = [gene]

    def len(self):
        return len(self)

    def findGene(self, site, lo, hi):
        """Given a site, e.g. 5246835, will return the gene chunk that the site
        falls into, otherwise None"""
        #binary search
        if hi < lo or lo >= len(self) or hi < 0: 
            return None
        else:
            mid = int((hi + lo)/2)
            #print self[mid]
            #print hi
            #print lo
            if self[mid].siteInGene(site): #Found!
                return self[mid]
            else:
                if self[mid].chunk()[0] < site:
                    return self.findGene(site, mid+1, hi)
                else:
                    return self.findGene(site, lo, mid-1)

#NOTE: this class is a subclass of GENES (overriding it's init) and reads
#in qdnaseq regions
class Qdnaseq(Gene):
    ct = 0

    #PRIMARY Constructor
    def __init__(self, line):
        #NOTE: hdr = chromosome\tstart\tend\tfeature\tCNV_scores--
        tmp = line.strip().split("\t")
        #NAME these regions by their current line-count in the file
        Qdnaseq.ct += 1
        self.name = self.name2 = Qdnaseq.ct
        self.chrom = "chr%s" % tmp[0]
        self.strand = "+" #all + strand
        self.txStart = int(tmp[1])
        self.txEnd = int(tmp[2])
        self.CNV = [float(s) for s in tmp[4:]]
        #NOTE: cds start is not used so it doesn't matter what we set it to
        self.cdsStart = int(tmp[1])
        self.cdsEnd = int(tmp[2])
        self.exons=[]
        self.weight = 1

    def getCNV(self):
        return self.CNV

    #def __str__(self):        
    #    ls = ["Loc: %s:%s-%s" % (self.chrom, self.txStart, self.txEnd)]
    #    return "\n".join(ls)

def main():
    optparser = OptionParser()
    optparser.add_option("-g", "--gt", help="UCSC geneTable")
    optparser.add_option("-i", "--igv", help="qdnaseq_segmented.igv file")
    optparser.add_option("-o", "--outpath", help="output path", default="./")
    #NOTE: for qdnaseq we will not consider promoter regions --upstream = 0
    optparser.add_option("-u", "--up", help="How many bps should we consider upstream of TSS? default: 0", default=0)
    optparser.add_option("-d", "--down", help="How many bps should we consider downstream of TTS? default: 0", default=0)
    optparser.add_option("-e", "--exon", help="How many bps should we consider downstream of TTS? default: 0", default=0)
    optparser.add_option("-t", "--gene", help="How many bps should we consider downstream of TTS? default: 0", default=0)

    (options, args) = optparser.parse_args(sys.argv)

    Gene.upstream = int(options.up)
    Gene.downstream = int(options.down)

    if not options.gt or not options.igv:
        print("USAGE: bedAnnotate.py -g [UCSC geneTable] -i [qdnaseq igv file] -u [upstream of TSS- default:0 (optional)] -d [downstream of TTS- default:0 (optional)]")
        sys.exit(-1)

    #l = "625	NM_000518	chr11	-	5246695	5248301	5246827	5248251	3	5246695,5247806,5248159,	5246956,5248029,5248301,	0	HBB	cmpl	cmpl	0,2,0,"
    outpath = options.outpath

    #READ in the segmented igv and the sampleNames (last set of lines in hdr)
    sampleNames = []
    #process the qdnaseq_segmented.igv file
    #NOTE: we assume the regions "tile" or partition the entire genome
    f = open(options.igv)
    igv_set = {}
    for l in f:
        #copy the headers back into the out_file
        if l.startswith("#"):
            continue
        if l.startswith("chromosome"): #get the sample names
            sampleNames = l.strip().split('\t')[4:]
            continue

        #MAIN body of qdnaseq_segmented.igv
        q = Qdnaseq(l)
        if q.chrom not in igv_set:
            igv_set[q.chrom] = Chromosome()
        igv_set[q.chrom].add(q)
    f.close()

    #take a peak and autoinfer the tiling size
    chr1 = igv_set['chr1']
    _tileSize = abs(chr1[0].txEnd - chr1[0].txStart) + 1

    out_file = open(os.path.join(outpath,"qdnaseq_genes.txt"),"w")
    out_file.write("Gene\t%s\n" % "\t".join(sampleNames))
    igv_file = open(os.path.join(outpath,"qdnaseq_genes.igv"),"w")
    #compose the igv header
    igv_file.write("#type=COPY_NUMBER\n#track coords=1\n")
    igv_file.write("%s\t" % "\t".join(["chromosome","start","end","feature"]))
    #ADD sample names
    igv_file.write("%s\n" % "\t".join(sampleNames))

    #NOW we have to go through each of the genes and try to find where it fits
    #in the igv_set
    f = open(options.gt)
    for l in f:
        if l.startswith("#"):
            continue
        gene = Gene(l)
        gene_len = gene.txEnd - gene.txStart

        #GENE is on a chromosome not represented in segmented.igv!--skip
        if gene.chrom not in igv_set:
            continue

        chrom = igv_set[gene.chrom]
        hi = len(chrom) - 1

        _tiles = Chromosome()
        start = gene.txStart
        while start < gene.txEnd:
            #FIND the IGV segment of the current start
            t = chrom.findGene(start, 0, hi)

            if not t: #skip
                start += _tileSize
                continue

            #ADD the tile to _tiles, but make sure we don't clip 
            #the end of the gene
            newT = copy.copy(t)
            #Change the start
            newT.txStart = start
            if newT.txEnd >= gene.txEnd:
                newT.txEnd = gene.txEnd
            #calc the weight
            newT.weight = float(newT.txEnd - newT.txStart) / gene_len
                    
            _tiles.add(newT)

            #ADVANCE:
            start = newT.txEnd + 1
            
        #now weight the CNVS:
        weightedCNV = [list(map(lambda c: t.weight*c, t.CNV)) for t in _tiles]
        #for t in _tiles:
        #    print(t.weight,t.CNV)
        #for w in weightedCNV:
        #    print(w)

        #now SUM the weights to find the value for each sample = scores
        scores = [0 for s in sampleNames]

        for (col, name) in enumerate(sampleNames):
            s = 0
            for row in range(len(_tiles)):
                s += weightedCNV[row][col]
            scores[col] = s
        #print(scores)

        #OUTPUT to txt
        out_cnvs = ["%.3f" % s for s in scores]
        #write to output
        out_line = [gene.name, gene.name2]
        out_line.extend(out_cnvs)
        #print("%s\n" % "\t".join(out_line))
        out_file.write("%s\n" % "\t".join(out_line))

        #OUTPUT to igv
        igv_line = [gene.chrom[3:], str(gene.txStart), str(gene.txEnd), "%s:%s" % (gene.name2, gene.name)]
        igv_line.extend(out_cnvs)
        igv_file.write("%s\n" % "\t".join(igv_line))

    f.close()
    out_file.close()
    igv_file.close()

if __name__=='__main__':
    main()

