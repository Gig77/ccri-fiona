#source("https://bioconductor.org/biocLite.R") ; biocLite("DiffBind")
library(DiffBind)

samples <- read.csv("/mnt/projects/fiona/data/sample-key-chipseq.csv")
samples <- samples[samples$SampleID != "35124_NALM6_RHD",]

db          <- dba(sampleSheet=samples, peakCaller="macs", peakFormat = "macs")
db.count    <- dba.count(db, minOverlap = 1, bRemoveDuplicates = T)
db.contrast <- dba.contrast(db.count, group1=samples$Factor=="ER", name1="ER", name2="RUNX1")
db.analyze  <- dba.analyze(db.contrast, bFullLibrarySize=TRUE, method=DBA_DESEQ2)
db.report   <- dba.report(db.analyze, th=0.05, bUsePval=T, method=DBA_DESEQ2)

db.report[db.report$Fold > 0,]
