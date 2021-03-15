#!/usr/bin/env Rscript

#NOTE: this was highly customized % graph--
#use plot_pbc.R instead for plotting RAW values
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))

plot_pbc_f <- function(nonChrM_in, png_out) {
    data <- read.csv(nonChrM_in, sep=",", header=TRUE, check.names=F )
    #print(data)
    #rownames(data) <- data[,1]
    x <- data.frame(data)
    colnames(x) <- c('Sample', 'Total_reads', 'nonChrM_reads')
    
    #print(x)
    ##NORMALIZE the values
    x$nonChrM_reads <- x$nonChrM_reads/x$Total_reads
    #drop TOTAL READS b/c now we normalized
    x <- subset(x, select=c("Sample", "nonChrM_reads"))
    x1 <- melt(x, id.var="Sample")
    #print(x1)

    png(png_out, width = 8, height = 8, unit="in",res=300)

    limits<-seq(0, 1.0, length.out=9)
    cust_labels<-limits

    #NOTE: lightblue-ish #91b6d4
    colors <- c(nonChrM_reads="steelblue", Total_reads="Grey")
    #KEY directive to make the barplot horizontal: coord_flip()
    ggplot(x1,aes(x=Sample, y=value, fill=variable)) +
        geom_bar(stat = "identity", position="identity") +
	geom_text(aes(label= paste0(round(value*100,2),"%"), y=value), size=3) +
        scale_y_continuous("",limits=c(0,1.0), labels=cust_labels, breaks=limits) +
        scale_fill_manual(values=colors) +
        labs(title="% of Non-ChrM Reads\n\n", x = "", y="") +
        guides(fill=guide_legend(title=NULL)) +
        theme_bw() +
        theme(axis.text.y = element_text(angle=0, hjust = 1, vjust=0.5, size=10), legend.position="top") +
        coord_flip()

    #dev.off()

}
args <- commandArgs( trailingOnly = TRUE )
arg_in = args[1]
arg_out = args[2]
#...
plot_pbc_f(arg_in, arg_out)

