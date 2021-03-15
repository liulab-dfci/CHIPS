#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))

plotGC_f <- function(gc_in, gc_out, gc_thumb_out) {
    data <- read.table(gc_in, sep="\t", header=TRUE, check.names=F )
    #print(data)

    x <- data.frame(data)
    colnames(x) <- c('GC','Count')
    gc <- as.vector(x['GC'])[[1]] #we need to access the first elm of list
    count <- as.vector(x['Count'])[[1]] #to make it a proper vector!

    #Get median of frequency dist.--
    #dist <- rep(gc, count)
    #print(median(dist))
    #HERE: Is when I realize that this is already calculated in fastqc.stat
    #and we don't need to do it again, but I'm keeping this for future ref!
    
    #full
    png(gc_out, height=8,width=8, unit='in', res=300)
    plot(x, type='l', col=rainbow(1)[1])
    #mark the median, 50 mark
    mark_v <- 50
    abline(v=mark_v, col="black")

    #NOTE: if you do this before  next plot it will create an Rplot.pdf file!!
    #DON'T do this!!
    #junk <- dev.off() 
    
    #thumbnail
    png(gc_thumb_out, height=85,width=150, unit='px')
    par(mar=c(0,0,0,0))
    plot(x, type='l',col=rainbow(1)[1],bty='n', xaxs='i', yaxs='i', ann=FALSE, xaxt='n', yaxt='n', bty='n')
    #mark the median, 50 mark
    abline(v=mark_v, col="black")
    junk <- dev.off()
    
}
args <- commandArgs( trailingOnly = TRUE )
arg_in = args[1]
arg_full_out = args[2]
arg_thumb_out = args[3]

plotGC_f(arg_in, arg_full_out, arg_thumb_out)

