library(GenomicRanges)
library(Vennerable)
library(gridExtra)

minscore <- 10
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

pdf("/mnt/projects/fiona/results/peak-overlaps.pdf", height=10)

# venn diagram AT2 vs REH

at2.and.reh <- round((length(unique(o.at2.reh@queryHits)) + length(unique(o.at2.reh@subjectHits)))/2)
at2.only <- length(at2.gr) - length(unique(o.at2.reh@queryHits))
reh.only <- length(reh.gr) - length(unique(o.at2.reh@subjectHits))
plot(Venn(SetNames = c("AT2", "REH"), Weight = c("10" = at2.only, "11" = at2.and.reh, "01" = reh.only)), doWeights=TRUE, type="circles")
grid.text(sprintf("ChIP-seq peak overlap (min. score = %d, min. overlap = %d bp)", minscore, minoverlap), vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

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

# NALM6 RUNX1 and E/R peak scores shared vs unique

par(mfrow=c(2,1))

hist(log(nalm6.runx1.gr$mcols.score[-unique(o.nalm6.runx1.er@queryHits)], 2), breaks=30, col=rgb(0,0,1,0.5), ylim=c(0,3000), xlim=c(1,10), xlab="Peak score (log2)", main="NALM6 RUNX1 peak score distribution")
hist(log(nalm6.runx1.gr$mcols.score[unique(o.nalm6.runx1.er@queryHits)], 2), breaks=40, col=rgb(1,0,0,0.5), ylim=c(0,3000), xlim=c(1,10), add=T)
legend("topright", c("shared with E/R", "only RUNX1"), fill=c("red", "blue"))

hist(log(nalm6.er.gr$mcols.score[-unique(o.nalm6.runx1.er@subjectHits)], 2), breaks=40, col=rgb(0,0,1,0.5), ylim=c(0,1500), xlim=c(1,10), xlab="Peak score (log2)", main="NALM6 E/R peak score distribution")
hist(log(nalm6.er.gr$mcols.score[unique(o.nalm6.runx1.er@subjectHits)], 2), breaks=40, col=rgb(1,0,0,0.5), ylim=c(0,1500), xlim=c(1,10), add=T)
legend("topright", c("shared with RUNX1", "only E/R"), fill=c("red", "blue"))

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

dev.off()
