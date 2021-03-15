#Script to run qdnaseq analysis
# params:
#   sampleName- abitrary name to give this run
#   dataDir- path to where the bam files are
#   binFilePath- path to qdnaseq annot file
#   outputDir- dir to output
# outputs:
#

library("QDNAseq")
qdnaseq_analysis <- function(sampleName, dataDir, binFilePath, outputDir) {
    #readCounts
    bins <- readRDS(binFilePath)
    readCounts <- binReadCounts(bins, path=dataDir)

    #Filter readCounts
    readCountsFiltered <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE, chromosomes = c("Y", "MT"))
    readCountsFiltered <- estimateCorrection(readCountsFiltered)

    #GET copy numbers
    copyNumbers <- correctBins(readCountsFiltered)
    copyNumbersNormalized <- normalizeBins(copyNumbers)
    copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

    #Write out copy numbers
    exportBins(copyNumbersSmooth, file=paste0(outputDir,sampleName,".txt"))
    exportBins(copyNumbersSmooth, file=paste0(outputDir,sampleName,".igv"),format="igv")
    exportBins(copyNumbersSmooth, file=paste0(outputDir,sampleName,".bed"), format="bed")

    #SEGMENT
    copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
    copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
        
    #CALL copyNumbers
    copyNumbersCalled <- callBins(copyNumbersSegmented)

    #SAVE SEGMENTED and CALLED copy numbers
    exportBins(copyNumbersSegmented, file=paste0(outputDir,sampleName,"_segmented.igv"), format="igv", type="segments")
    exportBins(copyNumbersCalled, file=paste0(outputDir,sampleName,"_calls.igv"), format="igv", logTransform=FALSE, type="calls")

    #PLOT
    pdf(file=paste0(outputDir,sampleName,".pdf"))
    plot(copyNumbersSmooth)
    plot(copyNumbersSegmented)
    plot(copyNumbersCalled)
    dev.off()
    
    gcts.copynumber <- as.data.frame(copyNumbersCalled@assayData$copynumber)
    exportSegments <- function(copyNumbersSegmented, sampleName){
        fileName=paste0(sampleName,".igv")
    }
}

args <- commandArgs( trailingOnly = TRUE )
arg_sName = args[1]
arg_dataDir = args[2]
arg_bin = args[3]
arg_out = args[4]

qdnaseq_analysis(arg_sName, arg_dataDir, arg_bin, arg_out);


