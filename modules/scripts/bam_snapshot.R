#Gviz bam coverage plot test
#ref: http://www.sthda.com/english/wiki/gviz-visualize-genomic-data#alignmentstrack

suppressMessages(library("Gviz"))
coverage_plot <- function(bam_f, paired, sample_name, gene, genome, zoom_factor, png_out) {
    #print(c(bam_f, paired, sample_name, gene, genome, png_out))
    alTrack <- AlignmentsTrack(bam_f, isPaired=paired, name=sample_name);
    bmt <- BiomartGeneRegionTrack(genome = genome, symbol=gene,
                                  filter = list(with_ox_refseq_mrna = TRUE),
                                  name="");
    
    start = attr(bmt,"start");
    end = attr(bmt,"end");
    chr= attr(bmt,"chromosome");
    #HACK: ACTB and TMPRSS2 genes are not using the correct start/end
    #so we need to correct for them
    if (gene== "ACTB") { #**ERROR: need to correct for mouse location; non-hg19
        start=5566779
        end=5570232
        chr="chr7"
    } else if (gene=="TMPRSS2") {
        start=42836478
        end=42879992
        chr="chr21"
    }
    
    len=end-start;
    #print(c(start,end,len))

    #NOTE: we want 2x zoom out, so we adjust the start and end by 50%
    #zoom_f = floor(0.50*len);
    #NOTE: we want 3x zoom out, so we adjust the start and end by 100%
    #GENERALLY: zoom_f = (zoom_factor - 1)/2 * len e.g. for 2x: (2-1)/2 *len
    zoom_f = floor((zoom_factor - 1.00)/2*len);
    #print(zoom_f);
    #checking the start location is in bounts
    z_start = if( start - zoom_f > 0) start-zoom_f else 1;
    #print(z_start)
    
    png(file = png_out, width = 900, height = 600, units="px", pointsize = 8);
    plotTracks(c(alTrack, bmt), from = z_start, to = end + zoom_f, chromosome = chr, type = "coverage", main=gene, sizes=c(280, 20));
    dev.off();
}

args <- commandArgs( trailingOnly = TRUE )
arg_bam = args[1]
arg_isPaired = args[2] == TRUE; #convert string to boolean
arg_sample_name = args[3]
arg_gene = args[4]
arg_genome = args[5]
arg_zoom_factor = as.numeric(args[6])
arg_out = args[7]


coverage_plot(arg_bam, arg_isPaired, arg_sample_name, arg_gene, arg_genome, arg_zoom_factor, arg_out);

