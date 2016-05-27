#source("https://bioconductor.org/biocLite.R") ; biocLite("DiffBind")
library(DiffBind)

samples <- read.csv("/mnt/projects/fiona/data/sample-key-chipseq.csv")
samples <- samples[samples$SampleID != "35124_NALM6_RHD",]

# setup using sample sheet
db <- dba(sampleSheet=samples, peakCaller="macs", peakFormat = "macs")

# calc/load counts
countfile <- "/mnt/projects/fiona/results/diffbind.counts.R"
if(file.exists(countfile)) {
  load(countfile)
} else {
  db.count    <- dba.count(db, minOverlap = 1, bRemoveDuplicates = T)
  save(db.count, file=countfile)
}

# get diff peaks
db.contrast <- dba.contrast(db.count, group1=samples$Factor=="ER", name1="ER", name2="RUNX1")
db.analyze  <- dba.analyze(db.contrast, bFullLibrarySize=TRUE, method=DBA_DESEQ2)
db.report   <- dba.report(db.analyze, th=1, bUsePval=T, method=DBA_DESEQ2)

# filter, sort
#diffpeak <- as.data.frame(db.report[db.report$Fold > 0,])
diffpeak <- as.data.frame(db.report)
diffpeak <- diffpeak[order(diffpeak$Fold, decreasing = T),]
diffpeak <- diffpeak[!grepl("^GL", diffpeak$seqnames),]

# find best original peak overlapping with merged peak region
library(GenomicRanges)

diffbind.gr <- makeGRangesFromDataFrame(diffpeak, keep.extra.columns = T)

reh <- read.delim("/mnt/projects/fiona/results/macs/ChIP24_REH_ER_peaks.bed") ; names(reh) <- c("seqnames", "start", "end", "id", "score")
summit <- read.delim("/mnt/projects/fiona/results/macs/ChIP24_REH_ER_summits.bed", check.names = F, stringsAsFactors = F, header = F) ; names(summit) <- c("seqnames", "summit_pos", "end", "id", "score")
summit$peak_id <- gsub("MACS_summit", "MACS_peak", summit$id)
reh <- merge(reh, summit[,c("peak_id", "summit_pos")], by.x="id", by.y="peak_id", all.x=T)
reh.gr <- makeGRangesFromDataFrame(reh, keep.extra.columns = T)

at2 <- read.delim("/mnt/projects/fiona/results/macs/ChIP24_AT2_ER_peaks.bed") ; names(at2) <- c("seqnames", "start", "end", "id", "score")
summit <- read.delim("/mnt/projects/fiona/results/macs/ChIP24_AT2_ER_summits.bed", check.names = F, stringsAsFactors = F, header = F) ; names(summit) <- c("seqnames", "summit_pos", "end", "id", "score")
summit$peak_id <- gsub("MACS_summit", "MACS_peak", summit$id)
at2 <- merge(at2, summit[,c("peak_id", "summit_pos")], by.x="id", by.y="peak_id", all.x=T)
at2.gr <- makeGRangesFromDataFrame(at2, keep.extra.columns = T)

nalm6er <- read.delim("/mnt/projects/fiona/results/macs/ChIP23_NALM6_ER_peaks.bed") ; names(nalm6er) <- c("seqnames", "start", "end", "id", "score")
summit <- read.delim("/mnt/projects/fiona/results/macs/ChIP23_NALM6_ER_summits.bed", check.names = F, stringsAsFactors = F, header = F) ; names(summit) <- c("seqnames", "summit_pos", "end", "id", "score")
summit$peak_id <- gsub("MACS_summit", "MACS_peak", summit$id)
nalm6er <- merge(nalm6er, summit[,c("peak_id", "summit_pos")], by.x="id", by.y="peak_id", all.x=T)
nalm6er.gr <- makeGRangesFromDataFrame(nalm6er, keep.extra.columns = T)

nalm6runx1 <- read.delim("/mnt/projects/fiona/results/macs/ChIP22_NALM6_RUNX1_peaks.bed") ; names(nalm6runx1) <- c("seqnames", "start", "end", "id", "score")
summit <- read.delim("/mnt/projects/fiona/results/macs/ChIP22_NALM6_RUNX1_summits.bed", check.names = F, stringsAsFactors = F, header = F) ; names(summit) <- c("seqnames", "summit_pos", "end", "id", "score")
summit$peak_id <- gsub("MACS_summit", "MACS_peak", summit$id)
nalm6runx1 <- merge(nalm6runx1, summit[,c("peak_id", "summit_pos")], by.x="id", by.y="peak_id", all.x=T)
nalm6runx1.gr <- makeGRangesFromDataFrame(nalm6runx1, keep.extra.columns = T)

o.reh <- findOverlaps(diffbind.gr, reh.gr, minoverlap = 100)
o.at2 <- findOverlaps(diffbind.gr, at2.gr, minoverlap = 100)
o.nalm6er <- findOverlaps(diffbind.gr, nalm6er.gr, minoverlap = 100)
o.nalm6runx1 <- findOverlaps(diffbind.gr, nalm6runx1.gr, minoverlap = 100)

o.combined <- data.frame(diffbind=o.reh@queryHits, cl="REH", id=reh$id[o.reh@subjectHits], seqnames=reh$seqnames[o.reh@subjectHits], start=reh$start[o.reh@subjectHits], end=reh$end[o.reh@subjectHits], summit=reh$summit_pos[o.reh@subjectHits], score=reh$score[o.reh@subjectHits])
o.combined <- rbind(o.combined, data.frame(diffbind=o.at2@queryHits, cl="AT2", id=at2$id[o.at2@subjectHits], seqnames=at2$seqnames[o.at2@subjectHits], start=at2$start[o.at2@subjectHits], end=at2$end[o.at2@subjectHits], summit=at2$summit_pos[o.at2@subjectHits], score=at2$score[o.at2@subjectHits]))
o.combined <- rbind(o.combined, data.frame(diffbind=o.nalm6er@queryHits, cl="NALM6", id=nalm6er$id[o.nalm6er@subjectHits], seqnames=nalm6er$seqnames[o.nalm6er@subjectHits], start=nalm6er$start[o.nalm6er@subjectHits], end=nalm6er$end[o.nalm6er@subjectHits], summit=nalm6er$summit_pos[o.nalm6er@subjectHits], score=nalm6er$score[o.nalm6er@subjectHits]))
o.combined <- rbind(o.combined, data.frame(diffbind=o.nalm6runx1@queryHits, cl="RUNX1", id=nalm6runx1$id[o.nalm6runx1@subjectHits], seqnames=nalm6runx1$seqnames[o.nalm6runx1@subjectHits], start=nalm6runx1$start[o.nalm6runx1@subjectHits], end=nalm6runx1$end[o.nalm6runx1@subjectHits], summit=nalm6runx1$summit_pos[o.nalm6runx1@subjectHits], score=nalm6runx1$score[o.nalm6runx1@subjectHits]))
o.combined <- o.combined[order(o.combined$diffbind, -o.combined$score),]
o.combined <- o.combined[!duplicated(o.combined$diffbind),]

# write peaks and summits file
peaks <- data.frame(seqnames=o.combined$seqnames, start=o.combined$start, end=o.combined$end, id=paste0("MACS_peak_", 1:nrow(o.combined)), fold=as.numeric(sprintf("%.2f", diffpeak$Fold[o.combined$diffbind])), FDR=sprintf("%.2g", diffpeak$FDR[o.combined$diffbind]))
write.table(peaks, "/mnt/projects/fiona/results/macs/diffbind_ER_vs_RUNX_peaks.bed.part", col.names = F, row.names = F, quote = F, sep = "\t")

summits <- data.frame(seqnames=o.combined$seqnames, start=o.combined$summit, end=o.combined$summit+1, id=paste0("MACS_summit_", 1:nrow(o.combined)), fold=as.numeric(sprintf("%.2f", diffpeak$Fold[o.combined$diffbind])), FDR=sprintf("%.2g", diffpeak$FDR[o.combined$diffbind]))
write.table(summits, "/mnt/projects/fiona/results/macs/diffbind_ER_vs_RUNX_summits.bed", col.names = F, row.names = F, quote = F, sep = "\t")
