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
UPDATED: 2018-02-23: re-geared to generate exon.bed and promoter.bed ref_files
-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
UPDATED 2016-10-13: ported to python3 AND to split the input into 
_exon.bed, _promoter.bed, etc.
-------------------------------------------------------------------------------
A program to output the sites of exons and promoters (given a refGene table)

Input: 
refGenes geneTable- expected format is: see UCSC hg19 refGenes table as example
bin name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts 
exonEnds score name2 cdsStartStat cdsEndStat exonFrames

Output:
exon.bed
promoter.bed

NOTE: intergenic/intronic are ignored b/c not needed as ref file

NOTE: intergenic is the 'other' category--meaning if it's neither promoter,
exon, or intron, then the default classification is "intergenic".  
THEREFORE: the percents will = 100%
"""
import os
import sys
from optparse import OptionParser

import bisect

def inRegion(reg, site):
    """given a region tuple, (start, end) 
    returns True if site is >= start <= end"""
    return site >= reg[0] and site <= reg[1]

class Gene():
    """Gene object--will help us store gene information and perform some fns"""
    #should make these command line parameters!
    upstream = 2000 #promoter = 2kb upstream of TSS
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
        self.name = tmp[1]
        self.name2 = tmp[-4]
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

    def classifySite(self, site):
        """Given a site, returns 
        "Intergenic" if it is OUTSIDE of gene--
        NOTE: since we don't know about other genes, "Intergenic" might not 
        be the true classification of the site

        "Promoter" if it is in upstream(2kb) - TSS
        "Exon" if it is in one of the exon regions
        "Intron" otherwise
        """
        def inPromoter(s):
            #VERSION 1-WHERE 5'UTR is NOT considered in promoter
            if self.strand == "+":
                p = (self.txStart - Gene.upstream, self.txStart) 
            else:
                p = (self.txEnd, self.txEnd + Gene.upstream)
            
            #VERSION 2- Where 5'UTR is in promoter
            #if self.strand == "+":
            #    p = (self.txStart - Gene.upstream, self.cdsStart) 
            #else:
            #    p = (self.cdsEnd, self.txEnd + Gene.upstream)

            #return s >= p[0] and s <= p[1]
            return inRegion(p, s)

        def inExons(s):
            for e in self.exons:
                #if s >= e[0] and s <= e[1]:
                if inRegion(e, s):
                    return True
            return False

        if not self.siteInGene(site):  ## discuss with len 201486
            return "Intergenic"
        elif inPromoter(site):
            return "Promoter"
        elif inExons(site):
            return "Exon"
        else:
            return "Intron"

    def __str__(self):        
        ls = ["Gene: %s\%s" % (self.name,self.name2),
              "Loc: %s:%s-%s" % (self.chrom, self.txStart, self.txEnd),
              "Strand: %s" % self.strand,
              "Coding: %s-%s" % (self.cdsStart,self.cdsEnd),
              "Exons: %s" % self.exons]
        return "\n".join(ls)

    #Comparisons based on genomic location--assumming other is on same chr
    def __eq__(self, other):
        return self.chunk() == other.chunk()
    def __ne__(self, other):
        return self.chunk() == other.chunk()  ## discuss with len 201486
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
    def __init__(self, *args):
        list.__init__(self, *args)

    def add(self, gene):
        """adds the gene into the sorted list"""
        i = bisect.bisect(self, gene)
        self[i:i] = [gene]

    def len(self):
        return len(self)

    def findGene(self, site, lo, hi):
        """Given a site, e.g. 5246835, will return the gene chunk that the site
        falls into, otherwise None"""
        #binary search
        #if hi < lo or lo >= self.len or hi < 0: ## need discuss with len 201486
        if hi < lo or lo >= len(self) or hi < 0: ## need discuss with len 201486
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

def main():
    optparser = OptionParser()
    optparser.add_option("-g", "--gt", help="UCSC geneTable")
    optparser.add_option("-l", "--lengths", help="lengths file") #used to cull out random stuff from refGene
    optparser.add_option("-n", "--basename", help="basename for the output")
    optparser.add_option("-o", "--outpath", help="output path", default="./")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.gt:
        print("USAGE: parseRefGene.py -g [UCSC geneTable] -l [lengths file] -o output path -n basename")
        sys.exit(-1)

    outpath = options.outpath
    basename = "%s." % options.basename if options.basename else ""

    #read the genetable
    f = open(options.gt)
    genome = {}
    for l in f:
        if l.startswith("#"):
            continue
        g = Gene(l)
        if g.chrom not in genome:
            genome[g.chrom] = Chromosome()
        genome[g.chrom].add(g)
    f.close()

    #IF lengths file is given, then use that to cull the genome of random stuff
    validChrs = list(genome.keys())
    if (options.lengths):
        f = open(options.lengths)
        validChrs = [l.strip().split("\t")[0] for l in f]
        f.close()
        #print(validChrs)
        #CULL geneTable:
        for k in list(genome.keys()):
            if k not in validChrs:
                del genome[k]

    #OUTPUT:
    out_exons = open(os.path.join(outpath, "%sexon.bed" % basename), "w")
    out_prom = open(os.path.join(outpath, "%spromoter.bed" % basename), "w")
    #ITERATE over entire genome --in order of validKeys, if given
    for chrom in validChrs:
        if chrom in genome:
            for gene in genome[chrom]:
                #PROMOTER--IF +, go upstrem (minus); else downstream
                if gene.strand == "+":
                    (start, end) = (gene.txStart - 2000, gene.txStart)
                else:
                    (start, end) = (gene.txStart, gene.txStart + 2000)
                prom_l = "\t".join([chrom, str(start), str(end), "%s" % gene.name,"0", gene.strand])
                out_prom.write("%s\n" % prom_l)
                #EXONS
                for (i, exon) in enumerate(gene.exons):
                    exon_l = "\t".join([chrom, str(exon[0]), str(exon[1]), "%s_exon%s" % (gene.name, str(i+1)),"0", gene.strand])
                    out_exons.write("%s\n" % exon_l)

if __name__=='__main__':
    main()

