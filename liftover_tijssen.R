# source("https://bioconductor.org/biocLite.R") ; biocLite("rtracklayer")
library(rtracklayer)
library(reshape2)
library(tidyr)
library(GenomicRanges)

# Tijssen peak positions from TableS1 --> GenomicRanges
t <- read.delim("/mnt/projects/fiona/data/Tijssen_2011_TableS1_peak_coordinates.txt", stringsAsFactors = F, check.names = F)
t <- melt(t, measure.vars = names(t), na.rm = T, variable.name = "factor")
t <- t[t$value != "",]
t <- separate(t, col=value, into = c("chr", "start", "end"), sep="_")
t$start <- as.numeric(t$start)
t$end <- as.numeric(t$end)
t <- makeGRangesFromDataFrame(t, keep.extra.columns = T)

# liftover hg18 to hg19
chain <- import.chain("/mnt/projects/generic/data/ucsc/hg18ToHg19.over.chain")
t.hg19 <- liftOver(t, chain)@unlistData

write.table(as.data.frame(t.hg19), "/mnt/projects/fiona/results/Tijssen_2011_TableS1_peak_coordinates.hg19.txt", col.names = T, row.names = F, quote = F, sep="\t")
