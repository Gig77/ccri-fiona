options(warn=1)
library(optparse)

option_list <- list(
  make_option("--peak-file", type="character", help="Annotated ChIP-seq peaks from Homer"),
  make_option("--summit-file", type="character", help="Peak summits (MACS output, BED format)"),
  make_option("--out-file", type="character", help="Output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

#opt <- data.frame('peak-file'= "/mnt/projects/fiona/results/homer/runx1_peaks.annotated.tsv", 'summit-file'= "/mnt/projects/fiona/results/macs/runx1_summits.bed", 'out-file'= "/mnt/projects/fiona/results/homer/runx1_peaks.annotated.with-expr.tsv", stringsAsFactors=F, check.names=F)
#opt <- data.frame('peak-file'= "/mnt/projects/fiona/results/homer/rhd_peaks.annotated.tsv", 'summit-file'= "/mnt/projects/fiona/results/macs/rhd_summits.bed", 'out-file'= "/mnt/projects/fiona/results/homer/rhd_peaks.annotated.with-expr.tsv", stringsAsFactors=F, check.names=F)
stopifnot(!is.null(opt$'peak-file'))
stopifnot(!is.null(opt$'out-file'))
stopifnot(!is.null(opt$'summit-file'))

peaks <- read.delim(opt$'peak-file', check.names = F, stringsAsFactors = F)

# merge summit coordinate
summit <- read.delim(opt$'summit-file', check.names = F, stringsAsFactors = F, header = F)
names(summit) <- c("chr", "summit_pos", "end", "id", "score")
summit$peak_id <- gsub("MACS_summit", "MACS_peak", summit$id)
peaks <- merge(peaks, summit[,c("peak_id", "summit_pos")], by.x=1, by.y=1, all.x=T)

# convert motif distances from distances relative to peak start to distances relative to peak summit
convertMotifDist <- function(colname) {
  motif.dists <- strsplit(as.character(peaks[,colname]), "\\(.*?\\),?", perl=T)
  for (i in 1:length(motif.dists)) {
    d <- as.numeric(motif.dists[[i]])
    if (length(d) > 0) {
      peaks[i,colname] <- paste0(peaks$Start[i]+d-peaks$summit_pos[i]-1, collapse=",")
    } else {
      peaks[i,colname] <- NA
    }
  }
  return(peaks)
}

peaks <- convertMotifDist('RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Peak(sequence,strand,conservation)')
peaks$'RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs' <- ifelse(is.na(peaks$'RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Peak(sequence,strand,conservation)'), NA, sapply(strsplit(peaks$'RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Peak(sequence,strand,conservation)', ","), length))

peaks <- convertMotifDist('ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Peak(sequence,strand,conservation)')
peaks$'ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs' <- ifelse(is.na(peaks$'ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Peak(sequence,strand,conservation)'), NA, sapply(strsplit(peaks$'ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Peak(sequence,strand,conservation)', ","), length))

peaks <- convertMotifDist('EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Peak(sequence,strand,conservation)')
peaks$'EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs' <- ifelse(is.na(peaks$'EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Peak(sequence,strand,conservation)'), NA, sapply(strsplit(peaks$'EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Peak(sequence,strand,conservation)', ","), length))

peaks <- convertMotifDist('GATA3(Zf)/iTreg-Gata3-ChIP-Seq(GSE20898)/Homer Distance From Peak(sequence,strand,conservation)')
peaks$'GATA3(Zf)/iTreg-Gata3-ChIP-Seq(GSE20898)/Homer No. motifs' <- ifelse(is.na(peaks$'GATA3(Zf)/iTreg-Gata3-ChIP-Seq(GSE20898)/Homer Distance From Peak(sequence,strand,conservation)'), NA, sapply(strsplit(peaks$'GATA3(Zf)/iTreg-Gata3-ChIP-Seq(GSE20898)/Homer Distance From Peak(sequence,strand,conservation)', ","), length))

names(peaks) <- gsub("Distance From Peak(sequence,strand,conservation)", "Distance From Summit", names(peaks), fixed=T)
  


# ER vs. empty

ERvsEmpty <- read.delim("/mnt/projects/fiona/results/anduril/execute/deseqAnnotated_oeERvsEmpty/table.csv")
ERvsEmpty <- ERvsEmpty[!is.na(ERvsEmpty$Gene),c("Gene", "fc", "q")]
ERvsEmpty <- ERvsEmpty[order(ERvsEmpty$q),]
ERvsEmpty <- ERvsEmpty[!duplicated(ERvsEmpty$Gene),]
names(ERvsEmpty) <- c("Gene", "fcERvsEmpty", "qERvsEmpty")
peaks.ann <- merge(peaks, ERvsEmpty, by.x = "Gene Name", by.y = "Gene", all.x = T)

# RHD vs. empty
RHDvsEmpty <- read.delim("/mnt/projects/fiona/results/anduril/execute/deseqAnnotated_oeRHDvsEmpty/table.csv")
RHDvsEmpty <- RHDvsEmpty[!is.na(RHDvsEmpty$Gene),c("Gene", "fc", "q")]
RHDvsEmpty <- RHDvsEmpty[order(RHDvsEmpty$q),]
RHDvsEmpty <- RHDvsEmpty[!duplicated(RHDvsEmpty$Gene),]
names(RHDvsEmpty) <- c("Gene", "fcRHDvsEmpty", "qRHDvsEmpty")
peaks.ann <- merge(peaks.ann, RHDvsEmpty, by.x = "Gene Name", by.y = "Gene", all.x = T)

# annotate enhancers
library(GenomicRanges)
peaks.gr <- GRanges(seqnames=peaks.ann$Chr, ranges=IRanges(start=peaks.ann$summit_pos, end=peaks.ann$summit_pos), peak=peaks.ann[,2])

peaks.ann$overlaps_enhancer_in_celllines <- NA
celllines <- c("Gm12878", "H1hesc", "Helas3", "Hepg2", "Huvec", "K562")
for (cl in celllines) {
  segmentation.bed <- read.delim(paste0("/mnt/projects/fiona/data/enhancers/wgEncodeAwgSegmentationCombined", cl, ".bed"), header = F, stringsAsFactors = F)
  segmentation.bed <- segmentation.bed[segmentation.bed$V4 %in% c("E"),] # we are only interested in enhancers here
  segmentation.gr <- GRanges(seqnames=segmentation.bed$V1, ranges=IRanges(start=segmentation.bed$V2, end=segmentation.bed$V3))
  o <- findOverlaps(peaks.gr, segmentation.gr)
  peaks.ann$overlaps_enhancer_in_celllines[o@queryHits] <- ifelse(is.na(peaks.ann$overlaps_enhancer_in_celllines[o@queryHits]), cl, paste(peaks.ann$overlaps_enhancer_in_celllines[o@queryHits], cl, sep=","))
}

# sort and write output
peaks.ann <- peaks.ann[order(peaks.ann$'Peak Score', decreasing = T),]
write.table(peaks.ann, opt$'out-file', col.names=T, row.names=F, sep="\t", quote=F, na="")
