library(ChIPQC)

#source("https://bioconductor.org/biocLite.R")
#biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")

setwd("/mnt/projects/fiona/results/")

samples <- read.delim("/mnt/projects/fiona/results/sample_key_chipseq_exp2.csv")
qc <- ChIPQC(samples, annotation="hg19", chromosomes="1")
qc

# the following plots don't work (error message from ChIPQC package)

coveragehistogram(qc)
plotCoverageHist(qc, facetBy=c("SampleID"))
plotPrincomp(qc)
plotPeakProfile(qc)
