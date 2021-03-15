#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))

plot_frag_dist <- function(frag_in, frag_out, sample_name) {
    #CHECK for empty files
    info = file.info(frag_in)
    empty = (info$size == 0)

    if (! empty) {
        data <- read.table(frag_in, sep=",", header=FALSE, check.names=F )
        colnames(data) <- c('FragLen')

        #calc median
        med <- median(data$FragLen)

        #PLOT png
        png(frag_out, width = 8, height = 8, unit="in",res=300)
        hist(data$FragLen, breaks=50, xlab=paste0("Fragment Length\n median = ", med), main=paste0("Fragment Length Distribution for ", sample_name))
    
        abline(v=med, col="red", lwd=3)
    
        junk <- dev.off()
    } else {
        #ERROR, no data--just produce empty plot
        png(frag_out, width = 8, height = 8, unit="in",res=300)
        plot.new()
        title(main="ERROR: no fragment length info available")
        junk <- dev.off()
    }    
}
args <- commandArgs( trailingOnly = TRUE )
arg_in = args[1]
arg_out = args[2]
arg_name = args[3]
#...
plot_frag_dist(arg_in, arg_out, arg_name)

