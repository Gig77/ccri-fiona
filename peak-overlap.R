library(GenomicRanges)
library(Vennerable)
library(gridExtra)

minscore <- 2
minoverlap <- 100

at2 <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
at2 <- at2[!grepl("^GL", at2$Chr),]
at2 <- at2[at2$`Peak Score` >= minscore,]
at2.gr <- GRanges(seqnames = at2$Chr, ranges=IRanges(at2$Start, at2$End), mcols=data.frame(score=at2$`Peak Score`))

reh <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_REH_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
reh <- reh[!grepl("^GL", reh$Chr),]
reh <- reh[reh$`Peak Score` >= minscore,]
reh.gr <- GRanges(seqnames = reh$Chr, ranges=IRanges(reh$Start, reh$End), mcols=data.frame(score=reh$`Peak Score`))

nalm6.runx1 <- read.delim("/mnt/projects/fiona/results/homer/ChIP22_NALM6_RUNX1_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.runx1 <- nalm6.runx1[!grepl("^GL", nalm6.runx1$Chr),]
nalm6.runx1 <- nalm6.runx1[nalm6.runx1$`Peak Score` >= minscore,]
nalm6.runx1.gr <- GRanges(seqnames = nalm6.runx1$Chr, ranges=IRanges(nalm6.runx1$Start, nalm6.runx1$End), mcols=data.frame(score=nalm6.runx1$`Peak Score`))

nalm6.er <- read.delim("/mnt/projects/fiona/results/homer/ChIP23_NALM6_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.er <- nalm6.er[!grepl("^GL", nalm6.er$Chr),]
nalm6.er <- nalm6.er[nalm6.er$`Peak Score` >= minscore,]
nalm6.er.gr <- GRanges(seqnames = nalm6.er$Chr, ranges=IRanges(nalm6.er$Start, nalm6.er$End), mcols=data.frame(score=nalm6.er$`Peak Score`))

# find overlaps

o.at2.reh <- findOverlaps(at2.gr, reh.gr, minoverlap=minoverlap)
o.nalm6.runx1.er <- findOverlaps(nalm6.runx1.gr, nalm6.er.gr, minoverlap=minoverlap)
o.nalm6.at2 <- findOverlaps(nalm6.er.gr, at2.gr, minoverlap=minoverlap)
o.nalm6.reh <- findOverlaps(nalm6.er.gr, reh.gr, minoverlap=minoverlap)

pdf("/mnt/projects/fiona/results/peak-overlaps.pdf", height=10)

# venn diagram AT2 vs REH

at2.and.reh <- round((length(unique(o.at2.reh@queryHits)) + length(unique(o.at2.reh@subjectHits)))/2)
at2.only <- length(at2.gr) - length(unique(o.at2.reh@queryHits))
reh.only <- length(reh.gr) - length(unique(o.at2.reh@subjectHits))
plot(Venn(SetNames = c("AT2", "REH"), Weight = c("10" = at2.only, "11" = at2.and.reh, "01" = reh.only)), doWeights=TRUE, type="circles")
grid.text(sprintf("ChIP-seq peak overlap (min. score = %d, min. overlap = %d bp)", minscore, minoverlap), vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

# AT2 vs REH peak score distributions

plot(density(log(at2$`Peak Score`, 2)), xlim=c(1,10), col="blue", lwd=3, xlab="Peak score (log2)", main="AT2 and REH peak score distribution")
lines(density(log(reh$`Peak Score`, 2)), xlim=c(1,10), col="red", lwd=3)
legend("topright", c("AT2", "REH"), fill=c("blue", "red"))

# AT2 vs REH peak width distributions

plot(density(reh$End-reh$Start, adjust=2), xlim=c(1,1500), col="red", lwd=3, xlab="Peak width (bp)", main="AT2 and REH peak width distribution")
lines(density(at2$End-at2$Start, adjust=2), col="blue", lwd=3)
legend("topright", c("AT2", "REH"), fill=c("blue", "red"))

# AT2 and REH peak scores shared vs unique

par(mfrow=c(2,1))

hist(log(at2.gr$mcols.score[-unique(o.at2.reh@queryHits)], 2), breaks=30, col=rgb(0,0,1,0.5), ylim=c(0,500), xlim=c(1,10), xlab="Peak score (log2)", main="AT2 peak score distribution")
hist(log(at2.gr$mcols.score[unique(o.at2.reh@queryHits)], 2), breaks=30, col=rgb(1,0,0,0.5), ylim=c(0,500), xlim=c(1,10), add=T)
legend("topright", c("shared with REH", "only AT2"), fill=c("red", "blue"))

hist(log(reh.gr$mcols.score[-unique(o.at2.reh@subjectHits)], 2), breaks=30, col=rgb(0,0,1,0.5), ylim=c(0,500), xlim=c(1,10), xlab="Peak score (log2)", main="REH peak score distribution")
hist(log(reh.gr$mcols.score[unique(o.at2.reh@subjectHits)], 2), breaks=30, col=rgb(1,0,0,0.5), ylim=c(0,500), xlim=c(1,10), add=T)
legend("topright", c("shared with AT2", "only REH"), fill=c("red", "blue"))

par(mfrow=c(1,1))

# venn diagram NALM6 RUNX1 vs. E/R

nalm6.runx1.and.er <- round((length(unique(o.nalm6.runx1.er@queryHits)) + length(unique(o.nalm6.runx1.er@subjectHits)))/2)
nalm6.runx1.only <- length(nalm6.runx1.gr) - length(unique(o.nalm6.runx1.er@queryHits))
nalm6.er.only <- length(nalm6.er.gr) - length(unique(o.nalm6.runx1.er@subjectHits))
plot(Venn(SetNames = c("NALM6 RUNX1", "NALM6 ER"), Weight = c("10" = nalm6.runx1.only, "11" = nalm6.runx1.and.er, "01" = nalm6.er.only)), doWeights=TRUE, type="circles")
grid.text(sprintf("ChIP-seq peak overlap (min. score = %d, min. overlap = %d bp)", minscore, minoverlap), vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

# RUNX1 vs E/R peak score distributions

plot(density(log(nalm6.er$`Peak Score`, 2), adjust=2), xlim=c(1,10), col="blue", lwd=3, xlab="Peak score (log2)", main="NALM6 RUNX1 and E/R peak score distribution")
lines(density(log(nalm6.runx1$`Peak Score`, 2), adjust=2), xlim=c(1,10), col="red", lwd=3)
legend("topright", c("E/R", "RUNX1"), fill=c("blue", "red"))

# RUNX1 vs E/R peak width distributions

plot(density(nalm6.er$End-nalm6.er$Start, adjust=2), xlim=c(1,2500), col="blue", lwd=3, xlab="Peak width (bp)", main="NALM6 RUNX1 and E/R peak width distribution")
lines(density(nalm6.runx1$End-nalm6.runx1$Start, adjust=2), col="red", lwd=3)
legend("topright", c("E/R", "RUNX1"), fill=c("blue", "red"))

# NALM6 RUNX1 and E/R peak scores shared vs unique

par(mfrow=c(2,1))

hist(log(nalm6.runx1.gr$mcols.score[-unique(o.nalm6.runx1.er@queryHits)], 2), breaks=30, col=rgb(0,0,1,0.5), ylim=c(0,3000), xlim=c(1,10), xlab="Peak score (log2)", main="NALM6 RUNX1 peak score distribution")
hist(log(nalm6.runx1.gr$mcols.score[unique(o.nalm6.runx1.er@queryHits)], 2), breaks=40, col=rgb(1,0,0,0.5), ylim=c(0,3000), xlim=c(1,10), add=T)
legend("topright", c("shared with E/R", "only RUNX1"), fill=c("red", "blue"))

hist(log(nalm6.er.gr$mcols.score[-unique(o.nalm6.runx1.er@subjectHits)], 2), breaks=40, col=rgb(0,0,1,0.5), ylim=c(0,1500), xlim=c(1,10), xlab="Peak score (log2)", main="NALM6 E/R peak score distribution")
hist(log(nalm6.er.gr$mcols.score[unique(o.nalm6.runx1.er@subjectHits)], 2), breaks=40, col=rgb(1,0,0,0.5), ylim=c(0,1500), xlim=c(1,10), add=T)
legend("topright", c("shared with RUNX1", "only E/R"), fill=c("red", "blue"))

par(mfrow=c(1,1))

# venn diagram NALM6 E/R vs AT2 vs REH

#source("https://bioconductor.org/biocLite.R") ; biocLite("ChIPpeakAnno")
library(ChIPpeakAnno)
o.nalm6.at2.reh <- findOverlapsOfPeaks(nalm6.er.gr, at2.gr, reh.gr, minoverlap=minoverlap, ignore.strand=TRUE, connectedPeaks = "merge")
venn.counts <- o.nalm6.at2.reh$venn_cnt[,4] 
names(venn.counts) <- apply(o.nalm6.at2.reh$venn_cnt[,1:3], 1, function(x) paste0(x, collapse="")) 
venn <- compute.Venn(Venn(SetNames = c("NALM6 ER", "AT2", "REH"), Weight = venn.counts), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text(sprintf("ChIP-seq peak overlap (min. score = %d, min. overlap = %d bp)", minscore, minoverlap), vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

# score distributions NALM6 E/R vs AT2 vs REH

par(mfrow=c(3,1))

nalm6.er.unique <- as.integer(gsub("nalm6.er.gr__", "", unlist(o.nalm6.at2.reh$peaklist$nalm6.er.gr$peakNames)))
hist(log(nalm6.er.gr$mcols.score[nalm6.er.unique], 2), breaks=50, col=rgb(0,0,1,0.5), xlim=c(1,10), xlab="Peak score (log2)", main="NALM6 E/R peak score distribution")
hist(log(nalm6.er.gr$mcols.score[-nalm6.er.unique], 2), breaks=50, col=rgb(1,0,0,0.5), add=T)
legend("topright", c("only NALM6", "shared with AT or REH"), fill=c("blue", "red"))

at2.unique <- as.integer(gsub("at2.gr__", "", c(unlist(o.nalm6.at2.reh$peaklist$at2.gr$peakNames), grep("at2.gr", unlist(unlist(o.nalm6.at2.reh$peaklist$`at2.gr///reh.gr`$peakNames)), value=T))))
hist(log(at2.gr$mcols.score[at2.unique], 2), breaks=50, col=rgb(0,0,1,0.5), xlim=c(1,10), ylim=c(0, 450), xlab="Peak score (log2)", main="AT2 peak score distribution")
hist(log(at2.gr$mcols.score[-at2.unique], 2), breaks=50, col=rgb(1,0,0,0.5), add=T)
legend("topright", c("only AT2", "shared with NALM6 E/R"), fill=c("blue", "red"))

reh.unique <- as.integer(gsub("reh.gr__", "", c(unlist(o.nalm6.at2.reh$peaklist$reh.gr$peakNames), grep("reh.gr", unlist(unlist(o.nalm6.at2.reh$peaklist$`at2.gr///reh.gr`$peakNames)), value=T))))
hist(log(reh.gr$mcols.score[reh.unique], 2), breaks=50, col=rgb(0,0,1,0.5), xlim=c(1,10), ylim=c(0, 600), xlab="Peak score (log2)", main="REH peak score distribution")
hist(log(reh.gr$mcols.score[-reh.unique], 2), breaks=50, col=rgb(1,0,0,0.5), add=T)
legend("topright", c("only REH", "shared with NALM6 E/R"), fill=c("blue", "red"))

par(mfrow=c(1,1))

# scatter plot peak scores

o.df <- data.frame(runx1.hits=o.nalm6.runx1.er@queryHits, 
                   er.hits=o.nalm6.runx1.er@subjectHits, 
                   runx1.score=nalm6.runx1.gr$mcols.score[o.nalm6.runx1.er@queryHits],
                   er.score=nalm6.er.gr$mcols.score[o.nalm6.runx1.er@subjectHits])
o.df.best <- do.call(rbind, by(o.df, o.df$runx1.hits, function(x) x[which.max(x$er.score),]))
o.df.best <- do.call(rbind, by(o.df.best, o.df.best$er.hits, function(x) x[which.max(x$runx1.score),]))

#plot(o.df.best$runx1.score, o.df.best$er.score, log="xy", cex=0.3, col=rgb(0,0,0,0.3))

library(ggplot2)
ggplot(data=o.df.best, aes(log(runx1.score), log(er.score))) +
  geom_point(alpha=0.3, size=1) + 
  stat_density2d(aes(fill=..level.., alpha=..level..), geom = "polygon", colour = "black") +
  scale_fill_continuous(low="blue", high="red") +
#  geom_smooth(method=lm, linetype=2, colour="red", se=F) +
  guides(alpha="none") +
  labs(color="Density", fill="Density", x="Peak score NALM6 RUNX1 (log scale)", y="Peak score NALM6 E/R (log scale)") + 
  theme_bw() +
  coord_fixed() +
  theme(legend.position=c(0,1), legend.justification=c(0,1)) +
  ggtitle("Comparison of RUNX1 and E/R binding affinity")

#spear <- cor.test(log(o.df.best$er.score), log(o.df.best$runx1.score), method="spearman")
#text(2, 800, sprintf("Spearman correlation = %.2f", spear$estimate), adj=0)


# distance to TSS

plot(density(at2$`Distance to TSS`/1000), xlim=c(-100,100), col="blue", lwd=3, xlab="Distance to nearest TSS in kb", main="REH+AT2 peak distance to nearest TSS")
lines(density(reh$`Distance to TSS`/1000), xlim=c(-100,100), col="red", lwd=3)
legend("topright", c("AT2", "REH"), fill=c("blue", "red"))

plot(density(at2$`Distance to TSS`[unique(o.at2.reh@queryHits)]/1000), xlim=c(-100,100), col="blue", lwd=3, xlab="Distance to nearest TSS in kb", main="AT2 peak distance to nearest TSS")
lines(density(at2$`Distance to TSS`[-unique(o.at2.reh@queryHits)]/1000), xlim=c(-100,100), col="red", lwd=3)
legend("topright", c("shared with REH", "only AT2"), fill=c("blue", "red"))

plot(density(reh$`Distance to TSS`[unique(o.at2.reh@subjectHits)]/1000), xlim=c(-100,100), col="blue", lwd=3, xlab="Distance to nearest TSS in kb", main="REH peak distance to nearest TSS")
lines(density(reh$`Distance to TSS`[-unique(o.at2.reh@subjectHits)]/1000), xlim=c(-100,100), col="red", lwd=3)
legend("topright", c("shared with AT2", "only REH"), fill=c("blue", "red"))

plot(density(at2$`Distance to TSS`[at2$`Peak Score`<=30]/1000), xlim=c(-100,100), col="blue", lwd=3, xlab="Distance to nearest TSS in kb", main="AT2 peak distance to nearest TSS")
lines(density(at2$`Distance to TSS`[at2$`Peak Score`>30]/1000), xlim=c(-100,100), col="red", lwd=3)
legend("topright", c("Peak score <= 30", "Peak score > 30"), fill=c("blue", "red"))

plot(density(reh$`Distance to TSS`[reh$`Peak Score`<=30]/1000), xlim=c(-100,100), col="blue", lwd=3, xlab="Distance to nearest TSS in kb", main="REH peak distance to nearest TSS")
lines(density(reh$`Distance to TSS`[reh$`Peak Score`>30]/1000), xlim=c(-100,100), col="red", lwd=3)
legend("topright", c("Peak score <= 30", "Peak score > 30"), fill=c("blue", "red"))

plot(density(nalm6.er$`Distance to TSS`/1000), xlim=c(-100,100), col="blue", lwd=3, xlab="Distance to nearest TSS in kb", main="NALM6 peak distance to nearest TSS")
lines(density(nalm6.runx1$`Distance to TSS`/1000), xlim=c(-100,100), col="red", lwd=3)
legend("topright", c("E/R", "RUNX1"), fill=c("blue", "red"))

plot(density(nalm6.runx1$`Distance to TSS`[unique(o.nalm6.runx1.er@queryHits)]/1000), xlim=c(-100,100), col="blue", lwd=3, xlab="Distance to nearest TSS in kb", main="NALM6 RUNX1 peak distance to nearest TSS")
lines(density(nalm6.runx1$`Distance to TSS`[-unique(o.nalm6.runx1.er@queryHits)]/1000), xlim=c(-100,100), col="red", lwd=3)
legend("topright", c("shared with E/R", "only RUNX1"), fill=c("blue", "red"))

plot(density(nalm6.er$`Distance to TSS`[unique(o.nalm6.runx1.er@subjectHits)]/1000), xlim=c(-100,100), col="blue", lwd=3, xlab="Distance to nearest TSS in kb", main="NALM6 E/R peak distance to nearest TSS")
lines(density(nalm6.er$`Distance to TSS`[-unique(o.nalm6.runx1.er@subjectHits)]/1000), xlim=c(-100,100), col="red", lwd=3)
legend("topright", c("shared with RUNX1", "only E/R"), fill=c("blue", "red"))

plot(density(nalm6.runx1$`Distance to TSS`[nalm6.runx1$`Peak Score`>30]/1000), xlim=c(-100,100), col="red", lwd=3, xlab="Distance to nearest TSS in kb", main="NALM6 RUNX1 peak distance to nearest TSS")
lines(density(nalm6.runx1$`Distance to TSS`[nalm6.runx1$`Peak Score`<=30]/1000), xlim=c(-100,100), col="blue", lwd=3)
legend("topright", c("Peak score <= 30", "Peak score > 30"), fill=c("blue", "red"))

plot(density(nalm6.er$`Distance to TSS`[nalm6.er$`Peak Score`>30]/1000), xlim=c(-100,100), col="red", lwd=3, xlab="Distance to nearest TSS in kb", main="NALM6 E/R peak distance to nearest TSS")
lines(density(nalm6.er$`Distance to TSS`[nalm6.er$`Peak Score`<=30]/1000), xlim=c(-100,100), col="blue", lwd=3)
legend("topright", c("Peak score <= 30", "Peak score > 30"), fill=c("blue", "red"))

# promoter overlap

library(vcd)
t <- as.table(matrix(c(sum(grepl("promoter", at2$`Detailed Annotation`)), 
                       sum(!grepl("promoter", at2$`Detailed Annotation`)),
                       sum(grepl("promoter", reh$`Detailed Annotation`)), 
                       sum(!grepl("promoter", reh$`Detailed Annotation`))),
                     nrow = 2,
                     dimnames = list('Region'=c("Promoter+", "Promoter-"), 
                                     'Cell line'=c("AT2", "REH")))
)
mosaic(t, pop=F, main=sprintf("p=%.2g", fisher.test(t)$p.value))
labeling_cells(text=t)(t)

t <- as.table(matrix(c(sum(grepl("promoter", nalm6.runx1$`Detailed Annotation`)), 
                       sum(!grepl("promoter", nalm6.runx1$`Detailed Annotation`)),
                       sum(grepl("promoter", nalm6.er$`Detailed Annotation`)), 
                       sum(!grepl("promoter", nalm6.er$`Detailed Annotation`))),
                     nrow = 2,
                     dimnames = list('Region'=c("Promoter+", "Promoter-"), 
                                     'Cell line'=c("NALM6 RUNX1", "NALM6 E/R")))
)
mosaic(t, pop=F, main=sprintf("p=%.2g", fisher.test(t)$p.value))
labeling_cells(text=t)(t)

# enhancer overlap

t <- as.table(matrix(c(sum(at2$overlaps_enhancer_in_celllines!=""), 
                       sum(at2$overlaps_enhancer_in_celllines==""),
                       sum(reh$overlaps_enhancer_in_celllines!=""), 
                       sum(reh$overlaps_enhancer_in_celllines=="")),
                     nrow = 2,
                     dimnames = list('Region'=c("Enhancer+", "Enhancer-"), 
                                     'Cell line'=c("AT2", "REH")))
)
mosaic(t, pop=F, main=sprintf("p=%.2g", fisher.test(t)$p.value))
labeling_cells(text=t)(t)

t <- as.table(matrix(c(sum(nalm6.runx1$overlaps_enhancer_in_celllines!=""), 
                       sum(nalm6.runx1$overlaps_enhancer_in_celllines==""),
                       sum(nalm6.er$overlaps_enhancer_in_celllines!=""), 
                       sum(nalm6.er$overlaps_enhancer_in_celllines=="")),
                     nrow = 2,
                     dimnames = list('Region'=c("Enhancer+", "Enhancer-"), 
                                     'Cell line'=c("NALM6 RUNX1", "NALM6 E/R")))
)
mosaic(t, pop=F, main=sprintf("p=%.2g", fisher.test(t)$p.value))
labeling_cells(text=t)(t)

dev.off()

# ----------------------------------------------------
# test for independence of ChIP binding and differential expression
# ----------------------------------------------------

pdf("/mnt/projects/fiona/results/chip-vs-expression.pdf", height=10)

minPadj <- 0.3
minLogfc <- 0.5
max.dist.TSS.us <- 5000
max.dist.TSS.ds <- 2000

#=========================
# AT2 
#=========================

fuka.d20plus.AT2 <- read.delim("/mnt/projects/chrisi/results/fuka/telamlKD.REHandAT2.esetnsF.onlyG_late_AT2.annot.tsv")
fuka.d20plus.AT2 <- fuka.d20plus.AT2[,c("syms", "Padj.onlyG_late_AT2", "logFC.onlyG_late_AT2")]
names(fuka.d20plus.AT2) <- c("Gene", "fuka.kdER.d20plus.AT2.padj", "fuka.kdER.d20plus.AT2.logfc")

at2.er.fuka <- at2[order(at2$`Peak Score`, decreasing = T),c("Gene Name", "Peak Score", "Distance to TSS")]
at2.er.fuka <- at2.er.fuka[at2.er.fuka$`Distance to TSS` >= -max.dist.TSS.us & at2.er.fuka$`Distance to TSS` <= max.dist.TSS.ds,]
at2.er.fuka <- at2.er.fuka[!duplicated(at2.er.fuka$`Gene Name`),]
at2.er.fuka <- merge(at2.er.fuka, fuka.d20plus.AT2, by.x="Gene Name", by.y="Gene", all=T)
names(at2.er.fuka) <- c("Gene", "PeakScore", "DistTSS", "padj", "logfc")

at2.er.fuka$PeakScore[is.na(at2.er.fuka$PeakScore)] <- 0
at2.er.fuka$padj[is.na(at2.er.fuka$padj)] <- 1
at2.er.fuka$logfc[is.na(at2.er.fuka$logfc)] <- 0
at2.er.fuka$DistTSS[is.na(at2.er.fuka$DistTSS)] <- 0

# AT2 down-regulated

t.at2.downreg <- with(at2.er.fuka, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore > 0 & (padj > minPadj | logfc > -minLogfc)),
  sum(PeakScore == 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore == 0 & (padj > minPadj | logfc > -minLogfc))),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.ER.AT2'           = c("Peak+", "Peak-"),
                  'fuka.kdER.d20plus.AT2' = c("down-reg.", "other")))))

mosaic(t.at2.downreg, pop=F, 
       main="AT2 E/R ChIP vs. Fuka KD", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.at2.downreg)$estimate, fisher.test(t.at2.downreg)$p.value))
labeling_cells(text=t.at2.downreg, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.at2.downreg)

# AT2 up-regulated

t.at2.upreg <- with(at2.er.fuka, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & logfc >= minLogfc),
  sum(PeakScore > 0 & (padj > minPadj | logfc < minLogfc)),
  sum(PeakScore == 0 & padj <= minPadj & logfc >= minLogfc),
  sum(PeakScore == 0 & (padj > minPadj | logfc < minLogfc))),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.ER.AT2'           = c("Peak+", "Peak-"),
                  'fuka.kdER.d20plus.AT2' = c("up-reg.", "other")))))

mosaic(t.at2.upreg, pop=F, 
       main="AT2 E/R ChIP vs. Fuka KD", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.at2.upreg)$estimate, fisher.test(t.at2.upreg)$p.value))
labeling_cells(text=t.at2.upreg, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.at2.upreg)

# AT2 up- or down-regulated

t.at2.diffreg <- with(at2.er.fuka, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & abs(logfc) >= minLogfc),
  sum(PeakScore > 0 & (padj > minPadj | abs(logfc) < minLogfc)),
  sum(PeakScore == 0 & padj <= minPadj & abs(logfc) >= minLogfc),
  sum(PeakScore == 0 & (padj > minPadj | abs(logfc) < minLogfc))),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.ER.AT2'           = c("Peak+", "Peak-"),
                  'fuka.kdER.d20plus.AT2' = c("diff-reg.", "other")))))

mosaic(t.at2.diffreg, pop=F, 
       main="AT2 E/R ChIP vs. Fuka KD", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.at2.diffreg)$estimate, fisher.test(t.at2.diffreg)$p.value))
labeling_cells(text=t.at2.diffreg, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.at2.diffreg)

# AT2 up- vs. down-regulated

t.at2.upvsdn <- with(at2.er.fuka, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore > 0 & padj <= minPadj & logfc >= minLogfc),
  sum(PeakScore == 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore == 0 & padj <= minPadj & logfc >= minLogfc)),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.ER.AT2'           = c("Peak+", "Peak-"),
                  'fuka.kdER.d20plus.AT2' = c("down-reg.", "up-reg.")))))

mosaic(t.at2.upvsdn, pop=F, 
       main="AT2 E/R ChIP vs. Fuka KD", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.at2.upvsdn)$estimate, fisher.test(t.at2.upvsdn)$p.value))
labeling_cells(text=t.at2.upvsdn, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.at2.upvsdn)

#x <- abs(at2.er.fuka$logfc[at2.er.fuka$padj <= minPadj & at2.er.fuka$PeakScore > 0])
#y <- log(at2.er.fuka$PeakScore[at2.er.fuka$padj <= minPadj & at2.er.fuka$PeakScore > 0], 2)
#plot(x, y)
#abline(lm(y~x), col="red", lwd="3")

#=========================
# REH
#=========================

fuka.d20plus.REH <- read.delim("/mnt/projects/chrisi/results/fuka/telamlKD.REHandAT2.esetnsF.onlyG_late_REH.annot.tsv")
fuka.d20plus.REH <- fuka.d20plus.REH[,c("syms", "Padj.onlyG_late_REH", "logFC.onlyG_late_REH")]
names(fuka.d20plus.REH) <- c("Gene", "fuka.kdER.d20plus.REH.padj", "fuka.kdER.d20plus.REH.logfc")

reh.er.fuka <- reh[order(reh$`Peak Score`, decreasing = T),c("Gene Name", "Peak Score", "Distance to TSS")]
reh.er.fuka <- reh.er.fuka[reh.er.fuka$`Distance to TSS` >= -max.dist.TSS.us & reh.er.fuka$`Distance to TSS` <= max.dist.TSS.ds,]
reh.er.fuka <- reh.er.fuka[!duplicated(reh.er.fuka$`Gene Name`),]
reh.er.fuka <- merge(reh.er.fuka, fuka.d20plus.REH, by.x="Gene Name", by.y="Gene", all=T)
names(reh.er.fuka) <- c("Gene", "PeakScore", "DistTSS", "padj", "logfc")

reh.er.fuka$PeakScore[is.na(reh.er.fuka$PeakScore)] <- 0
reh.er.fuka$padj[is.na(reh.er.fuka$padj)] <- 1
reh.er.fuka$logfc[is.na(reh.er.fuka$logfc)] <- 0
reh.er.fuka$DistTSS[is.na(reh.er.fuka$DistTSS)] <- 0

# REH down-regulated

t.reh.downreg <- with(reh.er.fuka, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore > 0 & (padj > minPadj | logfc > -minLogfc)),
  sum(PeakScore == 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore == 0 & (padj > minPadj | logfc > -minLogfc))),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.ER.REH'           = c("Peak+", "Peak-"),
                  'fuka.kdER.d20plus.REH' = c("down-reg.", "other")))))

mosaic(t.reh.downreg, pop=F, 
       main="REH E/R ChIP vs. Fuka KD", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.reh.downreg)$estimate, fisher.test(t.reh.downreg)$p.value))
labeling_cells(text=t.reh.downreg, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.reh.downreg)

# REH up-regulated

t.reh.upreg <- with(reh.er.fuka, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & logfc >= minLogfc),
  sum(PeakScore > 0 & (padj > minPadj | logfc < minLogfc)),
  sum(PeakScore == 0 & padj <= minPadj & logfc >= minLogfc),
  sum(PeakScore == 0 & (padj > minPadj | logfc < minLogfc))),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.ER.REH'           = c("Peak+", "Peak-"),
                  'fuka.kdER.d20plus.REH' = c("up-reg.", "other")))))

mosaic(t.reh.upreg, pop=F, 
       main="REH E/R ChIP vs. Fuka KD", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.reh.upreg)$estimate, fisher.test(t.reh.upreg)$p.value))
labeling_cells(text=t.reh.upreg, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.reh.upreg)

# REH up- or down-regulated

t.reh.diffreg <- with(reh.er.fuka, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & abs(logfc) >= minLogfc),
  sum(PeakScore > 0 & (padj > minPadj | abs(logfc) < minLogfc)),
  sum(PeakScore == 0 & padj <= minPadj & abs(logfc) >= minLogfc),
  sum(PeakScore == 0 & (padj > minPadj | abs(logfc) < minLogfc))),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.ER.REH'           = c("Peak+", "Peak-"),
                  'fuka.kdER.d20plus.REH' = c("diff-reg.", "other")))))

mosaic(t.reh.diffreg, pop=F, 
       main="REH E/R ChIP vs. Fuka KD", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.reh.diffreg)$estimate, fisher.test(t.reh.diffreg)$p.value))
labeling_cells(text=t.reh.diffreg, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.reh.diffreg)

# REH up- vs. down-regulated

t.reh.upvsdn <- with(reh.er.fuka, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore > 0 & padj <= minPadj & logfc >= minLogfc),
  sum(PeakScore == 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore == 0 & padj <= minPadj & logfc >= minLogfc)),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.ER.REH'           = c("Peak+", "Peak-"),
                  'fuka.kdER.d20plus.REH' = c("down-reg.", "up-reg.")))))

mosaic(t.reh.upvsdn, pop=F, 
       main="REH E/R ChIP vs. Fuka KD", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.reh.upvsdn)$estimate, fisher.test(t.reh.upvsdn)$p.value))
labeling_cells(text=t.reh.upvsdn, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.reh.upvsdn)

#=========================
# AT2+REH combined 
#=========================

fuka.d20plus <- read.delim("/mnt/projects/chrisi/results/fuka/matAnn.telamlKD.REHandAT2.esetnsF.REH.AT2.balanced.annot.tsv")
fuka.d20plus <- fuka.d20plus[,c("syms", "Padj", "logFC")]
names(fuka.d20plus) <- c("Gene", "fuka.kdER.d20plus.padj", "fuka.kdER.d20plus.logfc")

at2.reh.er.fuka <- at2[o.at2.reh@queryHits,]
at2.reh.er.fuka <- at2.reh.er.fuka[order(at2.reh.er.fuka$`Peak Score`, decreasing = T),c("Gene Name", "Peak Score", "Distance to TSS", "RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit")]
at2.reh.er.fuka <- at2.reh.er.fuka[at2.reh.er.fuka$`Distance to TSS` >= -max.dist.TSS.us & at2.reh.er.fuka$`Distance to TSS` <= max.dist.TSS.ds,]
at2.reh.er.fuka <- at2.reh.er.fuka[!duplicated(at2.reh.er.fuka$`Gene Name`),]
at2.reh.er.fuka <- merge(at2.reh.er.fuka, fuka.d20plus, by.x="Gene Name", by.y="Gene", all=T)
names(at2.reh.er.fuka) <- c("Gene", "PeakScore", "DistTSS", "runx1.motif.dist", "padj", "logfc")

at2.reh.er.fuka$PeakScore[is.na(at2.reh.er.fuka$PeakScore)] <- 0
at2.reh.er.fuka$padj[is.na(at2.reh.er.fuka$padj)] <- 1
at2.reh.er.fuka$logfc[is.na(at2.reh.er.fuka$logfc)] <- 0
at2.reh.er.fuka$DistTSS[is.na(at2.reh.er.fuka$DistTSS)] <- 0
at2.reh.er.fuka$runx1.motif.dist[is.na(at2.reh.er.fuka$runx1.motif.dist)] <- ""

# AT2+REH down-regulated

t.at2.reh.downreg <- with(at2.reh.er.fuka, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore > 0 & (padj > minPadj | logfc > -minLogfc)),
  sum(PeakScore == 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore == 0 & (padj > minPadj | logfc > -minLogfc))),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.ER.AT2+REH'           = c("Peak+", "Peak-"),
                  'fuka.kdER.d20plus.AT2+REH' = c("down-reg.", "other")))))

mosaic(t.at2.reh.downreg, pop=F, 
       main="AT2+REH E/R ChIP vs. Fuka KD", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.at2.reh.downreg)$estimate, fisher.test(t.at2.reh.downreg)$p.value))
labeling_cells(text=t.at2.reh.downreg, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.at2.reh.downreg)

# AT2+REH up-regulated

t.at2.reh.upreg <- with(at2.reh.er.fuka, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & logfc >= minLogfc),
  sum(PeakScore > 0 & (padj > minPadj | logfc < minLogfc)),
  sum(PeakScore == 0 & padj <= minPadj & logfc >= minLogfc),
  sum(PeakScore == 0 & (padj > minPadj | logfc < minLogfc))),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.ER.AT2+REH'           = c("Peak+", "Peak-"),
                  'fuka.kdER.d20plus.AT2+REH' = c("up-reg.", "other")))))

mosaic(t.at2.reh.upreg, pop=F, 
       main="AT2+REH E/R ChIP vs. Fuka KD", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.at2.reh.upreg)$estimate, fisher.test(t.at2.reh.upreg)$p.value))
labeling_cells(text=t.at2.reh.upreg, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.at2.reh.upreg)

# AT2+REH up- or down-regulated

t.at2.reh.diffreg <- with(at2.reh.er.fuka, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & abs(logfc) >= minLogfc),
  sum(PeakScore > 0 & (padj > minPadj | abs(logfc) < minLogfc)),
  sum(PeakScore == 0 & padj <= minPadj & abs(logfc) >= minLogfc),
  sum(PeakScore == 0 & (padj > minPadj | abs(logfc) < minLogfc))),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.ER.AT2+REH'           = c("Peak+", "Peak-"),
                  'fuka.kdER.d20plus.AT2+REH' = c("diff-reg.", "other")))))

mosaic(t.at2.reh.diffreg, pop=F, 
       main="AT2+REH E/R ChIP vs. Fuka KD", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.at2.reh.diffreg)$estimate, fisher.test(t.at2.reh.diffreg)$p.value))
labeling_cells(text=t.at2.reh.diffreg, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.at2.reh.diffreg)

# AT2+REH up- vs. down-regulated

t.at2.reh.upvsdn <- with(at2.reh.er.fuka, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore > 0 & padj <= minPadj & logfc >= minLogfc),
  sum(PeakScore == 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore == 0 & padj <= minPadj & logfc >= minLogfc)),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.ER.AT2+REH'           = c("Peak+", "Peak-"),
                  'fuka.kdER.d20plus.AT2+REH' = c("down-reg.", "up-reg.")))))

mosaic(t.at2.reh.upvsdn, pop=F, 
       main="AT2+REH E/R ChIP vs. Fuka KD", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.at2.reh.upvsdn)$estimate, fisher.test(t.at2.reh.upvsdn)$p.value))
labeling_cells(text=t.at2.reh.upvsdn, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.at2.reh.upvsdn)

# AT2+REH with RUNX1 motif, down-regulated

t.at2.reh.wRUNX1motif.downreg <- with(at2.reh.er.fuka, as.table(matrix(c(
  sum(PeakScore > 0 & runx1.motif.dist != "" & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore > 0 & runx1.motif.dist != "" & (padj > minPadj | logfc > -minLogfc)),
  sum((PeakScore == 0 | runx1.motif.dist == "") & padj <= minPadj & logfc <= -minLogfc),
  sum((PeakScore == 0 | runx1.motif.dist == "") & (padj > minPadj | logfc > -minLogfc))),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.ER.AT2+REH'           = c("Peak+", "Peak-"),
                  'fuka.kdER.d20plus.AT2+REH' = c("down-reg.", "other")))))

mosaic(t.at2.reh.wRUNX1motif.downreg, pop=F, 
       main="AT2+REH E/R ChIP w/ RUNX1 motif vs. Fuka KD", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.at2.reh.wRUNX1motif.downreg)$estimate, fisher.test(t.at2.reh.wRUNX1motif.downreg)$p.value))
labeling_cells(text=t.at2.reh.wRUNX1motif.downreg, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.at2.reh.wRUNX1motif.downreg)

# AT2+REH with RUNX1 motif, up-regulated

t.at2.reh.wRUNX1motif.upreg <- with(at2.reh.er.fuka, as.table(matrix(c(
  sum(PeakScore > 0 & runx1.motif.dist != "" & padj <= minPadj & logfc >= minLogfc),
  sum(PeakScore > 0 & runx1.motif.dist != "" & (padj > minPadj | logfc < minLogfc)),
  sum((PeakScore == 0 | runx1.motif.dist == "") & padj <= minPadj & logfc >= minLogfc),
  sum((PeakScore == 0 | runx1.motif.dist == "") & (padj > minPadj | logfc < minLogfc))),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.ER.AT2+REH'           = c("Peak+", "Peak-"),
                  'fuka.kdER.d20plus.AT2+REH' = c("up-reg.", "other")))))

mosaic(t.at2.reh.wRUNX1motif.upreg, pop=F, 
       main="AT2+REH E/R ChIP w/ RUNX1 motif vs. Fuka KD", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.at2.reh.wRUNX1motif.upreg)$estimate, fisher.test(t.at2.reh.wRUNX1motif.upreg)$p.value))
labeling_cells(text=t.at2.reh.wRUNX1motif.upreg, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.at2.reh.wRUNX1motif.upreg)

#=========================
# NALM6
#=========================

fiona.nalm6.oeERvsEmptyB2 <- read.delim("/mnt/projects/fiona/results/anduril/execute/deseqAnnotated_oeERvsEmptyB2/table.csv", check.names = F, stringsAsFactors = F)
fiona.nalm6.oeERvsEmptyB2 <- fiona.nalm6.oeERvsEmptyB2[order(fiona.nalm6.oeERvsEmptyB2$q),c("Gene", "q", "fc")]
fiona.nalm6.oeERvsEmptyB2 <- fiona.nalm6.oeERvsEmptyB2[!is.na(fiona.nalm6.oeERvsEmptyB2$Gene) & !is.na(fiona.nalm6.oeERvsEmptyB2$q) & !duplicated(fiona.nalm6.oeERvsEmptyB2$Gene),]

nalm6.er.expr <- nalm6.er[order(nalm6.er$`Peak Score`, decreasing = T),c("Gene Name", "Peak Score", "Distance to TSS")]
nalm6.er.expr <- nalm6.er.expr[nalm6.er.expr$`Distance to TSS` >= -max.dist.TSS.us & nalm6.er.expr$`Distance to TSS` <= max.dist.TSS.ds,]
nalm6.er.expr <- nalm6.er.expr[!duplicated(nalm6.er.expr$`Gene Name`),]
nalm6.er.expr <- merge(nalm6.er.expr, fiona.nalm6.oeERvsEmptyB2, by.x="Gene Name", by.y="Gene", all=T)
names(nalm6.er.expr) <- c("Gene", "PeakScore", "DistTSS", "padj", "logfc")

nalm6.er.expr$PeakScore[is.na(nalm6.er.expr$PeakScore)] <- 0
nalm6.er.expr$padj[is.na(nalm6.er.expr$padj)] <- 1
nalm6.er.expr$logfc[is.na(nalm6.er.expr$logfc)] <- 0
nalm6.er.expr$DistTSS[is.na(nalm6.er.expr$DistTSS)] <- 0

# NALM-6 down-regulated genes

t.nalm6.downreg <- with(nalm6.er.expr, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore > 0 & (padj > minPadj | logfc > -minLogfc)),
  sum(PeakScore == 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore == 0 & (padj > minPadj | logfc > -minLogfc))),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.oeER.NALM6'   = c("Peak+", "Peak-"),
                  'RNAseq.oeER.NALM6' = c("down-reg.", "other")))))

mosaic(t.nalm6.downreg, pop=F, 
       main="NALM-6 E/R ChIP vs. RNA-seq", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.nalm6.downreg)$estimate, fisher.test(t.nalm6.downreg)$p.value))
labeling_cells(text=t.nalm6.downreg, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.nalm6.downreg)

# NALM-6 up-regulated genes

t.nalm6.upreg <- with(nalm6.er.expr, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & logfc >= minLogfc),
  sum(PeakScore > 0 & (padj > minPadj | logfc < minLogfc)),
  sum(PeakScore == 0 & padj <= minPadj & logfc >= minLogfc),
  sum(PeakScore == 0 & (padj > minPadj | logfc < minLogfc))),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.oeER.NALM6'   = c("Peak+", "Peak-"),
                  'RNAseq.oeER.NALM6' = c("up-reg.", "other")))))

mosaic(t.nalm6.upreg, pop=F, 
       main="NALM-6 E/R ChIP vs. RNA-seq", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.nalm6.upreg)$estimate, fisher.test(t.nalm6.upreg)$p.value))
labeling_cells(text=t.nalm6.upreg, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.nalm6.upreg)

# NALM-6 up- or down-regulated genes

t.nalm6.diffreg <- with(nalm6.er.expr, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & abs(logfc) >= minLogfc),
  sum(PeakScore > 0 & (padj > minPadj | abs(logfc) < minLogfc)),
  sum(PeakScore == 0 & padj <= minPadj & abs(logfc) >= minLogfc),
  sum(PeakScore == 0 & (padj > minPadj | abs(logfc) < minLogfc))),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.oeER.NALM6'   = c("Peak+", "Peak-"),
                  'RNAseq.oeER.NALM6' = c("diff-reg.", "other")))))

mosaic(t.nalm6.diffreg, pop=F, 
       main="NALM-6 E/R ChIP vs. RNA-seq", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.nalm6.diffreg)$estimate, fisher.test(t.nalm6.diffreg)$p.value))
labeling_cells(text=t.nalm6.diffreg, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.nalm6.diffreg)

# NALM-6 up- vs. down-regulated genes

t.nalm6.upvsdn <- with(nalm6.er.expr, as.table(matrix(c(
  sum(PeakScore > 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore > 0 & padj <= minPadj & logfc >= minLogfc),
  sum(PeakScore == 0 & padj <= minPadj & logfc <= -minLogfc),
  sum(PeakScore == 0 & padj <= minPadj & logfc >= minLogfc)),
  nrow = 2,
  byrow = T,
  dimnames = list('ChIP.oeER.NALM6'   = c("Peak+", "Peak-"),
                  'RNAseq.oeER.NALM6' = c("down-reg.", "up-reg.")))))

mosaic(t.nalm6.upvsdn, pop=F, 
       main="NALM-6 E/R ChIP vs. RNA-seq", 
       sub=sprintf("OR=%.2f p=%.2g", fisher.test(t.nalm6.upvsdn)$estimate, fisher.test(t.nalm6.upvsdn)$p.value))
labeling_cells(text=t.nalm6.upvsdn, gp_text = gpar(cex=1), rot = 0, margin = unit(0, "lines"))(t.nalm6.upvsdn)

#=============================
# enrichment summary box plot
#=============================

data.at2.reh.fuka <- data.frame("CL" = character(0), "direction" = character(0), "OR" = numeric(0), "CI.low" = numeric(0), "CI.high" = numeric(0), "p" = numeric(0), stringsAsFactors = F)

t <- fisher.test(t.at2.upreg) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("AT2", "E/R repressed", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))
t <- fisher.test(t.at2.downreg) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("AT2", "E/R activated", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))
#t <- fisher.test(t.at2.diffreg) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("AT2", "up- or down-regulated", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))
#t <- fisher.test(t.at2.upvsdn) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("AT2", "up- vs. down-regulated", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))

t <- fisher.test(t.reh.upreg) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("REH", "E/R repressed", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))
t <- fisher.test(t.reh.downreg) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("REH", "E/R activated", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))
#t <- fisher.test(t.reh.diffreg) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("REH", "up- or down-regulated", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))
#t <- fisher.test(t.reh.upvsdn) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("REH", "up- vs. down-regulated", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))

t <- fisher.test(t.at2.reh.upreg) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("AT2+REH", "E/R repressed", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))
t <- fisher.test(t.at2.reh.downreg) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("AT2+REH", "E/R activated", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))
#t <- fisher.test(t.at2.reh.diffreg) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("AT2+REH", "up- or down-regulated", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))
#t <- fisher.test(t.at2.reh.upvsdn) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("AT2+REH", "up- vs. down-regulated", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))

t <- fisher.test(t.at2.reh.wRUNX1motif.upreg) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("AT2+REH\nw/ motif", "E/R repressed", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))
t <- fisher.test(t.at2.reh.wRUNX1motif.downreg) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("AT2+REH\nw/ motif", "E/R activated", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))

t <- fisher.test(t.nalm6.upreg) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("NALM-6", "E/R activated", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))
t <- fisher.test(t.nalm6.downreg) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("NALM-6", "E/R repressed", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))
#t <- fisher.test(t.nalm6.diffreg) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("NALM-6", "up- or down-regulated", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))
#t <- fisher.test(t.nalm6.upvsdn) ; data.at2.reh.fuka = rbind(data.at2.reh.fuka, setNames(data.frame("NALM-6", "up- vs. down-regulated", t$estimate, t$conf.int[1], t$conf.int[2], t$p.value, stringsAsFactors = F), names(data.at2.reh.fuka)))

data.at2.reh.fuka$CL <- factor(data.at2.reh.fuka$CL, levels=c("NALM-6", "AT2", "REH", "AT2+REH", "AT2+REH\nw/ motif"))
data.at2.reh.fuka$direction <- factor(data.at2.reh.fuka$direction, levels=c("E/R activated", "E/R repressed", "up- or down-regulated", "up- vs. down-regulated"))

ggplot(data = data.at2.reh.fuka, aes(x = CL, y = log(OR, 2), fill = direction)) +   
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymax = log(CI.high, 2), ymin = log(CI.low, 2)), colour = "black", width = 0.3, position=position_dodge(.9)) +
  labs(fill="Direction", x="Cell line", y="Log2 Enrichment") + 
  ggtitle("Peak enrichment in differentially regulated genes") +
  theme_bw() +
  coord_fixed()

dev.off()