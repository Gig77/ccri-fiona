options(warn=1)
library(optparse)

option_list <- list(
  make_option("--homer-peak-file", type="character", help="Annotated ChIP-seq peaks from Homer"),
  make_option("--macs-peak-file", type="character", help="(MACS) BED peak file with additional columns to add to output"),
  make_option("--summit-file", type="character", help="Peak summits (MACS output, BED format)"),
  make_option("--out-file", type="character", help="Output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

#opt <- data.frame('homer-peak-file'= "/mnt/projects/fiona/results/homer/runx1_peaks.annotated.tsv", 'summit-file'= "/mnt/projects/fiona/results/macs/runx1_summits.bed", 'out-file'= "/mnt/projects/fiona/results/homer/runx1_peaks.annotated.with-expr.tsv", stringsAsFactors=F, check.names=F)
#opt <- data.frame('homer-peak-file'= "/mnt/projects/fiona/results/homer/ChIP22_NALM6_RUNX1_peaks.annotated.tsv", 'summit-file'= "/mnt/projects/fiona/results/macs/ChIP22_NALM6_RUNX1_summits.bed", 'out-file'= "/mnt/projects/fiona/results/homer/ChIP22_NALM6_RUNX1_peaks.annotated.with-expr.tsv", stringsAsFactors=F, check.names=F)
#opt <- data.frame('homer-peak-file'= "/mnt/projects/fiona/results/homer/diffbind_ER_vs_RUNX_peaks.annotated.tsv", 'macs-peak-file'= "/mnt/projects/fiona/results/macs/diffbind_ER_vs_RUNX_peaks.bed", 'summit-file'= "/mnt/projects/fiona/results/macs/diffbind_ER_vs_RUNX_summits.bed", 'out-file'= "/mnt/projects/fiona/results/homer/diffbind_ER_vs_RUNX_peaks.annotated.with-expr.tsv", stringsAsFactors=F, check.names=F)
stopifnot(!is.null(opt$'homer-peak-file'))
stopifnot(!is.null(opt$'out-file'))
stopifnot(!is.null(opt$'summit-file'))

peaks <- read.delim(opt$'homer-peak-file', check.names = F, stringsAsFactors = F)

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

peaks <- convertMotifDist('ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer Distance From Peak(sequence,strand,conservation)')
peaks$'ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs' <- ifelse(is.na(peaks$'ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer Distance From Peak(sequence,strand,conservation)'), NA, sapply(strsplit(peaks$'ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer Distance From Peak(sequence,strand,conservation)', ","), length))

names(peaks) <- gsub("Distance From Peak(sequence,strand,conservation)", "Distance From Summit", names(peaks), fixed=T)

peaks.ann <- peaks

#-----------------------------------------------------------------
# annotate enhancers
#-----------------------------------------------------------------

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

#-----------------------------------------------------------------
# annotate constitutive peaks (i.e. overlapping with RUNX1 peak) vs. de novo peaks (i.e. not overlapping with RUNX1 peak)
#-----------------------------------------------------------------

peaks.gr <- GRanges(seqnames=peaks.ann$Chr, ranges=IRanges(start=peaks.ann$Start, end=peaks.ann$End))
nalm6.runx1 <- read.delim("/mnt/projects/fiona/results/homer/ChIP22_NALM6_RUNX1_peaks.annotated.tsv", stringsAsFactors = F, check.names = F)
nalm6.runx1.gr <- GRanges(seqnames = nalm6.runx1$Chr, ranges=IRanges(nalm6.runx1$Start, nalm6.runx1$End))
o.nalm6.runx1 <- findOverlaps(peaks.gr, nalm6.runx1.gr, minoverlap=100, ignore.strand=TRUE)
better_peak <- peaks.ann$`Peak Score`[o.nalm6.runx1@queryHits] > nalm6.runx1$`Peak Score`[o.nalm6.runx1@subjectHits]
peaks.ann$runx1_overlap <- "de novo"
peaks.ann$runx1_overlap[unique(o.nalm6.runx1@queryHits[better_peak])] <- "constitutive_better"
peaks.ann$runx1_overlap[unique(o.nalm6.runx1@queryHits[!better_peak])] <- "constitutive_worse"

#-----------------------------------------------------------------
# annotate expression datasets
#-----------------------------------------------------------------

# ER vs. empty (experiment 1)

ERvsEmptyB1 <- read.delim("/mnt/projects/fiona/results/anduril/execute/deseqAnnotated_oeERvsEmptyB1/table.csv")
ERvsEmptyB1 <- ERvsEmptyB1[!is.na(ERvsEmptyB1$Gene),c("Gene", "fc", "q")]
ERvsEmptyB1 <- ERvsEmptyB1[order(ERvsEmptyB1$q),]
ERvsEmptyB1 <- ERvsEmptyB1[!duplicated(ERvsEmptyB1$Gene),]
names(ERvsEmptyB1) <- c("Gene", "fcERvsEmptyB1", "qERvsEmptyB1")
peaks.ann <- merge(peaks.ann, ERvsEmptyB1, by.x = "Gene Name", by.y = "Gene", all.x = T)

# RHD vs. empty (experiment 1)
RHDvsEmptyB1 <- read.delim("/mnt/projects/fiona/results/anduril/execute/deseqAnnotated_oeRHDvsEmptyB1/table.csv")
RHDvsEmptyB1 <- RHDvsEmptyB1[!is.na(RHDvsEmptyB1$Gene),c("Gene", "fc", "q")]
RHDvsEmptyB1 <- RHDvsEmptyB1[order(RHDvsEmptyB1$q),]
RHDvsEmptyB1 <- RHDvsEmptyB1[!duplicated(RHDvsEmptyB1$Gene),]
names(RHDvsEmptyB1) <- c("Gene", "fcRHDvsEmptyB1", "qRHDvsEmptyB1")
peaks.ann <- merge(peaks.ann, RHDvsEmptyB1, by.x = "Gene Name", by.y = "Gene", all.x = T)

# ER vs. empty (experiment 2)

ERvsEmptyB2 <- read.delim("/mnt/projects/fiona/results/anduril/execute/deseqAnnotated_oeERvsEmptyB2/table.csv")
ERvsEmptyB2 <- ERvsEmptyB2[!is.na(ERvsEmptyB2$Gene),c("Gene", "fc", "q")]
ERvsEmptyB2 <- ERvsEmptyB2[order(ERvsEmptyB2$q),]
ERvsEmptyB2 <- ERvsEmptyB2[!duplicated(ERvsEmptyB2$Gene),]
names(ERvsEmptyB2) <- c("Gene", "fcERvsEmptyB2", "qERvsEmptyB2")
peaks.ann <- merge(peaks.ann, ERvsEmptyB2, by.x = "Gene Name", by.y = "Gene", all.x = T)

# RHD vs. empty (experiment 2)
RHDvsEmptyB2 <- read.delim("/mnt/projects/fiona/results/anduril/execute/deseqAnnotated_oeRHDvsEmptyB2/table.csv")
RHDvsEmptyB2 <- RHDvsEmptyB2[!is.na(RHDvsEmptyB2$Gene),c("Gene", "fc", "q")]
RHDvsEmptyB2 <- RHDvsEmptyB2[order(RHDvsEmptyB2$q),]
RHDvsEmptyB2 <- RHDvsEmptyB2[!duplicated(RHDvsEmptyB2$Gene),]
names(RHDvsEmptyB2) <- c("Gene", "fcRHDvsEmptyB2", "qRHDvsEmptyB2")
peaks.ann <- merge(peaks.ann, RHDvsEmptyB2, by.x = "Gene Name", by.y = "Gene", all.x = T)

fuka.d20plus <- read.delim("/mnt/projects/chrisi/results/fuka/matAnn.telamlKD.REHandAT2.esetnsF.REH.AT2.balanced.annot.tsv")
fuka.d20plus <- fuka.d20plus[,c("syms", "Padj", "logFC")]
fuka.d20plus <- fuka.d20plus[order(fuka.d20plus$Padj),]
fuka.d20plus <- fuka.d20plus[!duplicated(fuka.d20plus$syms),]
names(fuka.d20plus) <- c("Gene", "fuka.kdER.REHandAT2.d20plus.padj", "fuka.kdER.REHandAT2.d20plus.logfc")
peaks.ann <- merge(peaks.ann, fuka.d20plus, by.x = "Gene Name", by.y = "Gene", all.x = T)

fuka.d20plus.REH <- read.delim("/mnt/projects/chrisi/results/fuka/telamlKD.REHandAT2.esetnsF.onlyG_late_REH.annot.tsv")
fuka.d20plus.REH <- fuka.d20plus.REH[,c("syms", "Padj.onlyG_late_REH", "logFC.onlyG_late_REH")]
fuka.d20plus.REH <- fuka.d20plus.REH[order(fuka.d20plus.REH$Padj.onlyG_late_REH),]
fuka.d20plus.REH <- fuka.d20plus.REH[!duplicated(fuka.d20plus.REH$syms),]
names(fuka.d20plus.REH) <- c("Gene", "fuka.kdER.REH.d20plus.padj", "fuka.kdER.REH.d20plus.logfc")
peaks.ann <- merge(peaks.ann, fuka.d20plus.REH, by.x = "Gene Name", by.y = "Gene", all.x = T)

fuka.d20plus.AT2 <- read.delim("/mnt/projects/chrisi/results/fuka/telamlKD.REHandAT2.esetnsF.onlyG_late_AT2.annot.tsv")
fuka.d20plus.AT2 <- fuka.d20plus.AT2[,c("syms", "Padj.onlyG_late_AT2", "logFC.onlyG_late_AT2")]
fuka.d20plus.AT2 <- fuka.d20plus.AT2[order(fuka.d20plus.AT2$Padj.onlyG_late_AT2),]
fuka.d20plus.AT2 <- fuka.d20plus.AT2[!duplicated(fuka.d20plus.AT2$syms),]
names(fuka.d20plus.AT2) <- c("Gene", "fuka.kdER.AT2.d20plus.padj", "fuka.kdER.AT2.d20plus.logfc")
peaks.ann <- merge(peaks.ann, fuka.d20plus.AT2, by.x = "Gene Name", by.y = "Gene", all.x = T)

boer.TA.vs.noTall <- read.delim("/mnt/projects/chrisi/data/RossBoer/NordischALL.esetnsF.annot.txt", check.names=F)
boer.TA.vs.noTall <- boer.TA.vs.noTall[,c("syms", "adjPval.TAvs.mean.noTall", "TAvs.mean.noTall")]
boer.TA.vs.noTall <- boer.TA.vs.noTall[order(boer.TA.vs.noTall$adjPval.TAvs.mean.noTall),]
boer.TA.vs.noTall <- boer.TA.vs.noTall[!duplicated(boer.TA.vs.noTall$syms),]
names(boer.TA.vs.noTall) <- c("Gene", "boer.TA.vs.noTall.padj", "boer.TA.vs.noTall.logfc")
peaks.ann <- merge(peaks.ann, boer.TA.vs.noTall, by.x = "Gene Name", by.y = "Gene", all.x = T)

boer.TA.vs.rest <- read.delim("/mnt/projects/chrisi/data/RossBoer/matAnn.GSE13351_BOER.eset_zfilt_th3_nsF.tsv", check.names=F)
boer.TA.vs.rest <- boer.TA.vs.rest[,c("syms", "adjP.TA_vs_rest", "TA_vs_rest")]
boer.TA.vs.rest <- boer.TA.vs.rest[order(boer.TA.vs.rest$adjP.TA_vs_rest),]
boer.TA.vs.rest <- boer.TA.vs.rest[!duplicated(boer.TA.vs.rest$syms),]
names(boer.TA.vs.rest) <- c("Gene", "boer.TA.vs.rest.padj", "boer.TA.vs.rest.logfc")
peaks.ann <- merge(peaks.ann, boer.TA.vs.rest, by.x = "Gene Name", by.y = "Gene", all.x = T)

ross <- read.delim("/mnt/projects/chrisi/data/RossBoer/ROSS2.2003.esetnsF.annot.txt", check.names=F)
ross <- ross[,c("syms", "adjPval.TAvs.mean_noTALL", "TAvs.mean_noTALL")]
ross <- ross[order(ross$adjPval.TAvs.mean_noTALL),]
ross <- ross[!duplicated(ross$syms),]
names(ross) <- c("Gene", "ross.TA.vs.noTall.padj", "ross.TA.vs.noTall.logfc")
peaks.ann <- merge(peaks.ann, ross, by.x = "Gene Name", by.y = "Gene", all.x = T)

chrisi <- read.delim("/mnt/projects/chrisi/results/deseq/zuber+strobl-expressed-vs-notexpressed.deseq2.chipseq-annotated.tsv", check.names=F)
chrisi <- chrisi[,c("hgnc_symbol", "padj", "log2FoldChange")]
chrisi <- chrisi[order(chrisi$padj),]
chrisi <- chrisi[!duplicated(chrisi$hgnc_symbol),]
names(chrisi) <- c("Gene", "chrisi.oeER.padj", "chrisi.oeER.logfc")
peaks.ann <- merge(peaks.ann, chrisi, by.x = "Gene Name", by.y = "Gene", all.x = T)

veronika.E1 <- read.delim("/mnt/projects/helena_veronika/results/anduril/execute/deseqAnnotated_shG1vsNT/table.csv")
veronika.E1 <- veronika.E1[,c("Gene", "q", "fc")]
veronika.E1 <- veronika.E1[order(veronika.E1$q),]
veronika.E1 <- veronika.E1[!duplicated(veronika.E1$Gene),]
names(veronika.E1) <- c("Gene", "veronika.E1.kdER.vs.empty.padj", "veronika.E1.kdER.vs.empty.logfc")
peaks.ann <- merge(peaks.ann, veronika.E1, by.x = "Gene Name", by.y = "Gene", all.x = T)

veronika.E2.d3 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseqAnnotated_ERvsNTd3/table.csv")
veronika.E2.d3 <- veronika.E2.d3[,c("Gene", "q", "fc")]
veronika.E2.d3 <- veronika.E2.d3[order(veronika.E2.d3$q),]
veronika.E2.d3 <- veronika.E2.d3[!duplicated(veronika.E2.d3$Gene),]
names(veronika.E2.d3) <- c("Gene", "veronika.E2.kdER.vs.empty.D3.padj", "veronika.E2.kdER.vs.empty.D3.logfc")
peaks.ann <- merge(peaks.ann, veronika.E2.d3, by.x = "Gene Name", by.y = "Gene", all.x = T)

veronika.E2.d8 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseqAnnotated_ERvsNTd8/table.csv")
veronika.E2.d8 <- veronika.E2.d8[,c("Gene", "q", "fc")]
veronika.E2.d8 <- veronika.E2.d8[order(veronika.E2.d8$q),]
veronika.E2.d8 <- veronika.E2.d8[!duplicated(veronika.E2.d8$Gene),]
names(veronika.E2.d8) <- c("Gene", "veronika.E2.kdER.vs.empty.D8.padj", "veronika.E2.kdER.vs.empty.D8.logfc")
peaks.ann <- merge(peaks.ann, veronika.E2.d8, by.x = "Gene Name", by.y = "Gene", all.x = T)

veronika.E2.d15 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseqAnnotated_ERvsNTd15/table.csv")
veronika.E2.d15 <- veronika.E2.d15[,c("Gene", "q", "fc")]
veronika.E2.d15 <- veronika.E2.d15[order(veronika.E2.d15$q),]
veronika.E2.d15 <- veronika.E2.d15[!duplicated(veronika.E2.d15$Gene),]
names(veronika.E2.d15) <- c("Gene", "veronika.E2.kdER.vs.empty.D15.padj", "veronika.E2.kdER.vs.empty.D15.logfc")
peaks.ann <- merge(peaks.ann, veronika.E2.d15, by.x = "Gene Name", by.y = "Gene", all.x = T)

helena <- read.delim("/mnt/projects/helena_veronika/results/anduril/execute/deseqAnnotated_oeERvsEmpty/table.csv")
helena <- helena[,c("Gene", "q", "fc")]
helena <- helena[order(helena$q),]
helena <- helena[!duplicated(helena$Gene),]
names(helena) <- c("Gene", "helena.oeER.vs.empty.padj", "helena.oeER.vs.empty.logfc")
peaks.ann <- merge(peaks.ann, helena, by.x = "Gene Name", by.y = "Gene", all.x = T)

kasia.L24vsC24 <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_L24vsC24/table.csv")
kasia.L24vsC24 <- kasia.L24vsC24[order(kasia.L24vsC24$q),]
kasia.L24vsC24 <- kasia.L24vsC24[!duplicated(kasia.L24vsC24$Gene),]
kasia.L24vsC24 <- kasia.L24vsC24[,c("Gene", "q", "fc")]
names(kasia.L24vsC24) <- c("Gene", "kasia.HDACi.L24vsC24.padj", "kasia.HDACi.L24vsC24.logfc")
peaks.ann <- merge(peaks.ann, kasia.L24vsC24, by.x = "Gene Name", by.y = "Gene", all.x = T)

kasia.H24vsC24 <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_H24vsC24/table.csv")
kasia.H24vsC24 <- kasia.H24vsC24[order(kasia.H24vsC24$q),]
kasia.H24vsC24 <- kasia.H24vsC24[!duplicated(kasia.H24vsC24$Gene),]
kasia.H24vsC24 <- kasia.H24vsC24[,c("Gene", "q", "fc")]
names(kasia.H24vsC24) <- c("Gene", "kasia.HDACi.H24vsC24.padj", "kasia.HDACi.H24vsC24.logfc")
peaks.ann <- merge(peaks.ann, kasia.H24vsC24, by.x = "Gene Name", by.y = "Gene", all.x = T)

kasia.L48vsC48 <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_L48vsC48/table.csv")
kasia.L48vsC48 <- kasia.L48vsC48[order(kasia.L48vsC48$q),]
kasia.L48vsC48 <- kasia.L48vsC48[!duplicated(kasia.L48vsC48$Gene),]
kasia.L48vsC48 <- kasia.L48vsC48[,c("Gene", "q", "fc")]
names(kasia.L48vsC48) <- c("Gene", "kasia.HDACi.L48vsC48.padj", "kasia.HDACi.L48vsC48.logfc")
peaks.ann <- merge(peaks.ann, kasia.L48vsC48, by.x = "Gene Name", by.y = "Gene", all.x = T)

kasia.H48vsC48 <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_H48vsC48/table.csv")
kasia.H48vsC48 <- kasia.H48vsC48[order(kasia.H48vsC48$q),]
kasia.H48vsC48 <- kasia.H48vsC48[!duplicated(kasia.H48vsC48$Gene),]
kasia.H48vsC48 <- kasia.H48vsC48[,c("Gene", "q", "fc")]
names(kasia.H48vsC48) <- c("Gene", "kasia.HDACi.H48vsC48.padj", "kasia.HDACi.H48vsC48.logfc")
peaks.ann <- merge(peaks.ann, kasia.H48vsC48, by.x = "Gene Name", by.y = "Gene", all.x = T)

kasia.H24vsL24 <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_H24vsL24/table.csv")
kasia.H24vsL24 <- kasia.H24vsL24[order(kasia.H24vsL24$q),]
kasia.H24vsL24 <- kasia.H24vsL24[!duplicated(kasia.H24vsL24$Gene),]
kasia.H24vsL24 <- kasia.H24vsL24[,c("Gene", "q", "fc")]
names(kasia.H24vsL24) <- c("Gene", "kasia.HDACi.H24vsL24.padj", "kasia.HDACi.H24vsL24.logfc")
peaks.ann <- merge(peaks.ann, kasia.H24vsL24, by.x = "Gene Name", by.y = "Gene", all.x = T)

kasia.H48vsL48 <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_H48vsL48/table.csv")
kasia.H48vsL48 <- kasia.H48vsL48[order(kasia.H48vsL48$q),]
kasia.H48vsL48 <- kasia.H48vsL48[!duplicated(kasia.H48vsL48$Gene),]
kasia.H48vsL48 <- kasia.H48vsL48[,c("Gene", "q", "fc")]
names(kasia.H48vsL48) <- c("Gene", "kasia.HDACi.H48vsL48.padj", "kasia.HDACi.H48vsL48.logfc")
peaks.ann <- merge(peaks.ann, kasia.H48vsL48, by.x = "Gene Name", by.y = "Gene", all.x = T)

kasia.H48vsH24 <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_H48vsH24/table.csv")
kasia.H48vsH24 <- kasia.H48vsH24[order(kasia.H48vsH24$q),]
kasia.H48vsH24 <- kasia.H48vsH24[!duplicated(kasia.H48vsH24$Gene),]
kasia.H48vsH24 <- kasia.H48vsH24[,c("Gene", "q", "fc")]
names(kasia.H48vsH24) <- c("Gene", "kasia.HDACi.H48vsH24.padj", "kasia.HDACi.H48vsH24.logfc")
peaks.ann <- merge(peaks.ann, kasia.H48vsH24, by.x = "Gene Name", by.y = "Gene", all.x = T)

kasia.L48vsL24 <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_L48vsL24/table.csv")
kasia.L48vsL24 <- kasia.L48vsL24[order(kasia.L48vsL24$q),]
kasia.L48vsL24 <- kasia.L48vsL24[!duplicated(kasia.L48vsL24$Gene),]
kasia.L48vsL24 <- kasia.L48vsL24[,c("Gene", "q", "fc")]
names(kasia.L48vsL24) <- c("Gene", "kasia.HDACi.L48vsL24.padj", "kasia.HDACi.L48vsL24.logfc")
peaks.ann <- merge(peaks.ann, kasia.L48vsL24, by.x = "Gene Name", by.y = "Gene", all.x = T)

kasia.LvsC <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_LvsC/table.csv")
kasia.LvsC <- kasia.LvsC[order(kasia.LvsC$q),]
kasia.LvsC <- kasia.LvsC[!duplicated(kasia.LvsC$Gene),]
kasia.LvsC <- kasia.LvsC[,c("Gene", "q", "fc")]
names(kasia.LvsC) <- c("Gene", "kasia.HDACi.LvsC.padj", "kasia.HDACi.LvsC.logfc")
peaks.ann <- merge(peaks.ann, kasia.LvsC, by.x = "Gene Name", by.y = "Gene", all.x = T)

kasia.HvsC <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_HvsC/table.csv")
kasia.HvsC <- kasia.HvsC[order(kasia.HvsC$q),]
kasia.HvsC <- kasia.HvsC[!duplicated(kasia.HvsC$Gene),]
kasia.HvsC <- kasia.HvsC[,c("Gene", "q", "fc")]
names(kasia.HvsC) <- c("Gene", "kasia.HDACi.HvsC.padj", "kasia.HDACi.HvsC.logfc")
peaks.ann <- merge(peaks.ann, kasia.HvsC, by.x = "Gene Name", by.y = "Gene", all.x = T)

kasia.HvsL <- read.delim("/mnt/projects/kasia/results/anduril/execute/deseqAnnotated_HvsL/table.csv")
kasia.HvsL <- kasia.HvsL[order(kasia.HvsL$q),]
kasia.HvsL <- kasia.HvsL[!duplicated(kasia.HvsL$Gene),]
kasia.HvsL <- kasia.HvsL[,c("Gene", "q", "fc")]
names(kasia.HvsL) <- c("Gene", "kasia.HDACi.HvsL.padj", "kasia.HDACi.HvsL.logfc")
peaks.ann <- merge(peaks.ann, kasia.HvsL, by.x = "Gene Name", by.y = "Gene", all.x = T)

#-----------------------------------------------------------------
# annotate other ChIP-seq datasets
#-----------------------------------------------------------------

# ChIP-seq Tijssen et al. 2011 (http://www.ncbi.nlm.nih.gov/pubmed/21571218) - primary human megakaryocytes
tijssen <- read.delim("/mnt/projects/fiona/results/Tijssen_2011_TableS1_peak_coordinates.hg19.txt")
tijssen <- makeGRangesFromDataFrame(tijssen, keep.extra.columns = T)
o.tijssen <- findOverlaps(peaks.gr, tijssen, minoverlap=100, ignore.strand=TRUE)
if (length(o.tijssen) > 0) {
  o.tijssen <- data.frame(peak=peaks[o.tijssen@queryHits,1], factor=tijssen[o.tijssen@subjectHits]$factor)
  o.tijssen <- o.tijssen[!duplicated(o.tijssen),]
  o.tijssen <- aggregate(factor~peak, paste, collapse="_", data=o.tijssen)
  o.tijssen$factor <- sapply(o.tijssen$factor, function(x) paste(sort(unique(strsplit(x, "_")[[1]])), collapse=","))
  peaks.ann <- merge(peaks.ann, o.tijssen, by.x = 2, by.y = 1, all.x = T)
  colnames(peaks.ann)[length(peaks.ann)] <- "Tijssen2011.chipseq.human.megakaryocytes"
} else {
  peaks.ann$Tijssen2011.chipseq.human.megakaryocytes <- NA
}

# ChIP-seq Wilson et al. 2010 (http://www.ncbi.nlm.nih.gov/pubmed/20887958) - mosue HPC-7 cells
wilson <- read.csv("/mnt/projects/chrisi/results/chipseq/Wilson_Gottgens_ChIPseq.txt",  stringsAsFactors=F, sep="\t", header=T, fill=T)
wilson <- wilson[!is.na(wilson$Runx1) & wilson$Runx1 != "",]
library(biomaRt)
human = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice" , dataset="hsapiens_gene_ensembl") # GRCh37, v75
mouse = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", dataset="mmusculus_gene_ensembl") # GRCm38, v75
humOrt <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values=wilson$Runx1, mart = mouse, attributesL = c("hgnc_symbol", "entrezgene"), martL = human)
wilson <- humOrt[!is.na(humOrt$EntrezGene.ID), c("HGNC.symbol", "MGI.symbol")]
wilson <- wilson[!duplicated(wilson),]
colnames(wilson) <- c("Gene", "wilson2010.chipseq.runx1.mouse")
wilson <- aggregate(wilson2010.chipseq.runx1.mouse~Gene, paste, collapse="|", data=wilson)
wilson <- wilson[wilson$Gene != "",]
peaks.ann <- merge(peaks.ann, wilson, by.x = "Gene Name", by.y = "Gene", all.x = T)

# ChIP-seq Niebuhr et. al 2013 (http://www.ncbi.nlm.nih.gov/pubmed/23704093) - BMiFLT3 (15-3) cells (proB cells from a murine BCP-ALL)
niebuhr <-  read.csv("/mnt/projects/chrisi/results/chipseq/Niebuhr_TableS3_Runx1 Peaks Called in ProB-Cells.txt", stringsAsFactors=F, sep="\t", header=T, fill=T)
#niebuhr <- niebuhr[niebuhr$dist_tss > -5000 & niebuhr$dist_tss < 1000,]
#niebuhr <- niebuhr[niebuhr$score >= 100,]
humOrt <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values=niebuhr$nearest.gene, mart = mouse, attributesL = c("hgnc_symbol", "entrezgene"), martL = human)
niebuhr <- humOrt[!is.na(humOrt$EntrezGene.ID), c("HGNC.symbol", "MGI.symbol")]
niebuhr <- niebuhr[!duplicated(niebuhr),]
colnames(niebuhr) <- c("Gene", "niebuhr2013.chipseq.runx1.mouse")
niebuhr <- aggregate(niebuhr2013.chipseq.runx1.mouse~Gene, paste, collapse="|", data=niebuhr)
niebuhr <- niebuhr[niebuhr$Gene != "",]
peaks.ann <- merge(peaks.ann, niebuhr, by.x = "Gene Name", by.y = "Gene", all.x = T)

# add additional columns present in initial (pre-annotation) BED file (didn't find a way to make HOMER preserve such columns when annotating peaks)
if (!is.null(opt$`macs-peak-file`)) {
  macs.peaks <- read.delim(opt$'macs-peak-file', check.names = F, stringsAsFactors = F, header = F)
  macs.peaks$V1 <- paste0("chr", macs.peaks$V1)
  macs.peaks$V2 <- macs.peaks$V2 + 1
  names(macs.peaks) <- paste0(opt$`macs-peak-file`, ".", names(macs.peaks))
  peaks.ann <- merge(peaks.ann, macs.peaks, by.x=c(2, 3, 4, 5, 7), by.y=c(4, 1, 2, 3, 5), all.x = T)
}

# sort and write output
peaks.ann <- peaks.ann[order(peaks.ann$'Peak Score', decreasing = T),]
write.table(peaks.ann, opt$'out-file', col.names=T, row.names=F, sep="\t", quote=F, na="")
