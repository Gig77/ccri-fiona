# read data
source("/mnt/projects/fiona/scripts/read-peaks.R")

# AT2

at2.motifs <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_AT2_ER.motif_hits.tsv", stringsAsFactors = F, check.names = F)
at2.motifs <- at2.motifs[order(at2.motifs$PositionID, at2.motifs$`Motif Name`, at2.motifs$MotifScore, -abs(at2.motifs$Offset), decreasing = T),]
at2.motifs.best <- at2.motifs[!duplicated(at2.motifs[,c("PositionID", "Motif Name")]),]

# AT2 RUNX

at2 <- merge(at2, at2.motifs.best[at2.motifs.best$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(at2)[names(at2) %in% c("Sequence", "MotifScore")] <- c("Sequence.RUNX", "MotifScore.RUNX")
at2$MotifScoreBin.RUNX <- "absent"
at2$MotifScoreBin.RUNX[!is.na(at2$MotifScore.RUNX)] <- "intermediate"
at2$MotifScoreBin.RUNX[at2$MotifScore.RUNX>=10] <- "high"
at2$MotifScoreBin.RUNX[at2$MotifScore.RUNX<=7] <- "low"
at2$MotifScoreBin.RUNX <- factor(at2$MotifScoreBin.RUNX, levels=c("absent", "low", "intermediate", "high"))

# AT2 ETS

at2 <- merge(at2, at2.motifs.best[at2.motifs.best$`Motif Name`=="ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(at2)[names(at2) %in% c("Sequence", "MotifScore")] <- c("Sequence.ETS", "MotifScore.ETS")
at2$MotifScoreBin.ETS <- "absent"
at2$MotifScoreBin.ETS[!is.na(at2$MotifScore.ETS)] <- "intermediate"
at2$MotifScoreBin.ETS[at2$MotifScore.ETS>=10] <- "high"
at2$MotifScoreBin.ETS[at2$MotifScore.ETS<=8] <- "low"
at2$MotifScoreBin.ETS <- factor(at2$MotifScoreBin.ETS, levels=c("absent", "low", "intermediate", "high"))

# AT2 EBF

at2 <- merge(at2, at2.motifs.best[at2.motifs.best$`Motif Name`=="EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(at2)[names(at2) %in% c("Sequence", "MotifScore")] <- c("Sequence.EBF", "MotifScore.EBF")
at2$MotifScoreBin.EBF <- "absent"
at2$MotifScoreBin.EBF[!is.na(at2$MotifScore.EBF)] <- "intermediate"
at2$MotifScoreBin.EBF[at2$MotifScore.EBF>=11] <- "high"
at2$MotifScoreBin.EBF[at2$MotifScore.EBF<=9.5] <- "low"
at2$MotifScoreBin.EBF <- factor(at2$MotifScoreBin.EBF, levels=c("absent", "low", "intermediate", "high"))

# REH

reh.motifs <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_REH_ER.motif_hits.tsv", stringsAsFactors = F, check.names = F)
reh.motifs <- reh.motifs[order(reh.motifs$PositionID, reh.motifs$`Motif Name`, reh.motifs$MotifScore, -abs(reh.motifs$Offset), decreasing = T),]
reh.motifs.best <- reh.motifs[!duplicated(reh.motifs[,c("PositionID", "Motif Name")]),]

# REH RUNX

reh <- merge(reh, reh.motifs.best[reh.motifs.best$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(reh)[names(reh) %in% c("Sequence", "MotifScore")] <- c("Sequence.RUNX", "MotifScore.RUNX")
reh$MotifScoreBin.RUNX <- "absent"
reh$MotifScoreBin.RUNX[!is.na(reh$MotifScore.RUNX)] <- "intermediate"
reh$MotifScoreBin.RUNX[reh$MotifScore.RUNX>=10] <- "high"
reh$MotifScoreBin.RUNX[reh$MotifScore.RUNX<=7] <- "low"
reh$MotifScoreBin.RUNX <- factor(reh$MotifScoreBin.RUNX, levels=c("absent", "low", "intermediate", "high"))

# REH ETS

reh <- merge(reh, reh.motifs.best[reh.motifs.best$`Motif Name`=="ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(reh)[names(reh) %in% c("Sequence", "MotifScore")] <- c("Sequence.ETS", "MotifScore.ETS")
reh$MotifScoreBin.ETS <- "absent"
reh$MotifScoreBin.ETS[!is.na(reh$MotifScore.ETS)] <- "intermediate"
reh$MotifScoreBin.ETS[reh$MotifScore.ETS>=10] <- "high"
reh$MotifScoreBin.ETS[reh$MotifScore.ETS<=8] <- "low"
reh$MotifScoreBin.ETS <- factor(reh$MotifScoreBin.ETS, levels=c("absent", "low", "intermediate", "high"))

# REH EBF

reh <- merge(reh, reh.motifs.best[reh.motifs.best$`Motif Name`=="EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(reh)[names(reh) %in% c("Sequence", "MotifScore")] <- c("Sequence.EBF", "MotifScore.EBF")
reh$MotifScoreBin.EBF <- "absent"
reh$MotifScoreBin.EBF[!is.na(reh$MotifScore.EBF)] <- "intermediate"
reh$MotifScoreBin.EBF[reh$MotifScore.EBF>=11] <- "high"
reh$MotifScoreBin.EBF[reh$MotifScore.EBF<=9.5] <- "low"
reh$MotifScoreBin.EBF <- factor(reh$MotifScoreBin.EBF, levels=c("absent", "low", "intermediate", "high"))

# NALM-6 ER

nalm6.er.motifs <- read.delim("/mnt/projects/fiona/results/motifs/ChIP23_NALM6_ER.motif_hits.tsv", stringsAsFactors = F, check.names = F)
nalm6.er.motifs <- nalm6.er.motifs[order(nalm6.er.motifs$PositionID, nalm6.er.motifs$`Motif Name`, nalm6.er.motifs$MotifScore, -abs(nalm6.er.motifs$Offset), decreasing = T),]
nalm6.er.motifs.best <- nalm6.er.motifs[!duplicated(nalm6.er.motifs[,c("PositionID", "Motif Name")]),]

# NALM-6 ER RUNX

nalm6.er <- merge(nalm6.er, nalm6.er.motifs.best[nalm6.er.motifs.best$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(nalm6.er)[names(nalm6.er) %in% c("Sequence", "MotifScore")] <- c("Sequence.RUNX", "MotifScore.RUNX")
nalm6.er$MotifScoreBin.RUNX <- "absent"
nalm6.er$MotifScoreBin.RUNX[!is.na(nalm6.er$MotifScore.RUNX)] <- "intermediate"
nalm6.er$MotifScoreBin.RUNX[nalm6.er$MotifScore.RUNX>=10] <- "high"
nalm6.er$MotifScoreBin.RUNX[nalm6.er$MotifScore.RUNX<=7] <- "low"
nalm6.er$MotifScoreBin.RUNX <- factor(nalm6.er$MotifScoreBin.RUNX, levels=c("absent", "low", "intermediate", "high"))

# NALM-6 ER ETS

nalm6.er <- merge(nalm6.er, nalm6.er.motifs.best[nalm6.er.motifs.best$`Motif Name`=="ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(nalm6.er)[names(nalm6.er) %in% c("Sequence", "MotifScore")] <- c("Sequence.ETS", "MotifScore.ETS")
nalm6.er$MotifScoreBin.ETS <- "absent"
nalm6.er$MotifScoreBin.ETS[!is.na(nalm6.er$MotifScore.ETS)] <- "intermediate"
nalm6.er$MotifScoreBin.ETS[nalm6.er$MotifScore.ETS>=10] <- "high"
nalm6.er$MotifScoreBin.ETS[nalm6.er$MotifScore.ETS<=8] <- "low"
nalm6.er$MotifScoreBin.ETS <- factor(nalm6.er$MotifScoreBin.ETS, levels=c("absent", "low", "intermediate", "high"))

# NALM-6 ER EBF

nalm6.er <- merge(nalm6.er, nalm6.er.motifs.best[nalm6.er.motifs.best$`Motif Name`=="EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(nalm6.er)[names(nalm6.er) %in% c("Sequence", "MotifScore")] <- c("Sequence.EBF", "MotifScore.EBF")
nalm6.er$MotifScoreBin.EBF <- "absent"
nalm6.er$MotifScoreBin.EBF[!is.na(nalm6.er$MotifScore.EBF)] <- "intermediate"
nalm6.er$MotifScoreBin.EBF[nalm6.er$MotifScore.EBF>=11] <- "high"
nalm6.er$MotifScoreBin.EBF[nalm6.er$MotifScore.EBF<=9.5] <- "low"
nalm6.er$MotifScoreBin.EBF <- factor(nalm6.er$MotifScoreBin.EBF, levels=c("absent", "low", "intermediate", "high"))

# NALM-6 RUNX1

nalm6.runx1.motifs <- read.delim("/mnt/projects/fiona/results/motifs/ChIP22_NALM6_RUNX1.motif_hits.tsv", stringsAsFactors = F, check.names = F)
nalm6.runx1.motifs <- nalm6.runx1.motifs[order(nalm6.runx1.motifs$PositionID, nalm6.runx1.motifs$`Motif Name`, nalm6.runx1.motifs$MotifScore, -abs(nalm6.runx1.motifs$Offset), decreasing = T),]
nalm6.runx1.motifs.best <- nalm6.runx1.motifs[!duplicated(nalm6.runx1.motifs[,c("PositionID", "Motif Name")]),]

# NALM-6 RUNX1 RUNX

nalm6.runx1 <- merge(nalm6.runx1, nalm6.runx1.motifs.best[nalm6.runx1.motifs.best$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(nalm6.runx1)[names(nalm6.runx1) %in% c("Sequence", "MotifScore")] <- c("Sequence.RUNX", "MotifScore.RUNX")
nalm6.runx1$MotifScoreBin.RUNX <- "absent"
nalm6.runx1$MotifScoreBin.RUNX[!is.na(nalm6.runx1$MotifScore.RUNX)] <- "intermediate"
nalm6.runx1$MotifScoreBin.RUNX[nalm6.runx1$MotifScore.RUNX>=10] <- "high"
nalm6.runx1$MotifScoreBin.RUNX[nalm6.runx1$MotifScore.RUNX<=7] <- "low"
nalm6.runx1$MotifScoreBin.RUNX <- factor(nalm6.runx1$MotifScoreBin.RUNX, levels=c("absent", "low", "intermediate", "high"))

# NALM-6 RUNX1 ETS

nalm6.runx1 <- merge(nalm6.runx1, nalm6.runx1.motifs.best[nalm6.runx1.motifs.best$`Motif Name`=="ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(nalm6.runx1)[names(nalm6.runx1) %in% c("Sequence", "MotifScore")] <- c("Sequence.ETS", "MotifScore.ETS")
nalm6.runx1$MotifScoreBin.ETS <- "absent"
nalm6.runx1$MotifScoreBin.ETS[!is.na(nalm6.runx1$MotifScore.ETS)] <- "intermediate"
nalm6.runx1$MotifScoreBin.ETS[nalm6.runx1$MotifScore.ETS>=10] <- "high"
nalm6.runx1$MotifScoreBin.ETS[nalm6.runx1$MotifScore.ETS<=8] <- "low"
nalm6.runx1$MotifScoreBin.ETS <- factor(nalm6.runx1$MotifScoreBin.ETS, levels=c("absent", "low", "intermediate", "high"))

# NALM-6 RUNX1 EBF

nalm6.runx1 <- merge(nalm6.runx1, nalm6.runx1.motifs.best[nalm6.runx1.motifs.best$`Motif Name`=="EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(nalm6.runx1)[names(nalm6.runx1) %in% c("Sequence", "MotifScore")] <- c("Sequence.EBF", "MotifScore.EBF")
nalm6.runx1$MotifScoreBin.EBF <- "absent"
nalm6.runx1$MotifScoreBin.EBF[!is.na(nalm6.runx1$MotifScore.EBF)] <- "intermediate"
nalm6.runx1$MotifScoreBin.EBF[nalm6.runx1$MotifScore.EBF>=11] <- "high"
nalm6.runx1$MotifScoreBin.EBF[nalm6.runx1$MotifScore.EBF<=9.5] <- "low"
nalm6.runx1$MotifScoreBin.EBF <- factor(nalm6.runx1$MotifScoreBin.EBF, levels=c("absent", "low", "intermediate", "high"))

# build dataframe

d <- with(at2, data.frame("CL" = "AT2", "Motif" = "RUNX", "PeakScore" = `Peak Score`, "MotifScore" = MotifScore.RUNX, "MotifScoreBin" = MotifScoreBin.RUNX))
d <- rbind(d, with(at2, data.frame("CL" = "AT2", "Motif" = "ETS", "PeakScore" = `Peak Score`, "MotifScore" = MotifScore.ETS, "MotifScoreBin" = MotifScoreBin.ETS)))
d <- rbind(d, with(at2, data.frame("CL" = "AT2", "Motif" = "EBF", "PeakScore" = `Peak Score`, "MotifScore" = MotifScore.EBF, "MotifScoreBin" = MotifScoreBin.EBF)))
d <- rbind(d, with(reh, data.frame("CL" = "REH", "Motif" = "RUNX", "PeakScore" = `Peak Score`, "MotifScore" = MotifScore.RUNX, "MotifScoreBin" = MotifScoreBin.RUNX)))
d <- rbind(d, with(reh, data.frame("CL" = "REH", "Motif" = "ETS", "PeakScore" = `Peak Score`, "MotifScore" = MotifScore.ETS, "MotifScoreBin" = MotifScoreBin.ETS)))
d <- rbind(d, with(reh, data.frame("CL" = "REH", "Motif" = "EBF", "PeakScore" = `Peak Score`, "MotifScore" = MotifScore.EBF, "MotifScoreBin" = MotifScoreBin.EBF)))
d <- rbind(d, with(nalm6.er, data.frame("CL" = "NALM6 ER", "Motif" = "RUNX", "PeakScore" = `Peak Score`, "MotifScore" = MotifScore.RUNX, "MotifScoreBin" = MotifScoreBin.RUNX)))
d <- rbind(d, with(nalm6.er, data.frame("CL" = "NALM6 ER", "Motif" = "ETS", "PeakScore" = `Peak Score`, "MotifScore" = MotifScore.ETS, "MotifScoreBin" = MotifScoreBin.ETS)))
d <- rbind(d, with(nalm6.er, data.frame("CL" = "NALM6 ER", "Motif" = "EBF", "PeakScore" = `Peak Score`, "MotifScore" = MotifScore.EBF, "MotifScoreBin" = MotifScoreBin.EBF)))
d <- rbind(d, with(nalm6.runx1, data.frame("CL" = "NALM6 RUNX1", "Motif" = "RUNX", "PeakScore" = `Peak Score`, "MotifScore" = MotifScore.RUNX, "MotifScoreBin" = MotifScoreBin.RUNX)))
d <- rbind(d, with(nalm6.runx1, data.frame("CL" = "NALM6 RUNX1", "Motif" = "ETS", "PeakScore" = `Peak Score`, "MotifScore" = MotifScore.ETS, "MotifScoreBin" = MotifScoreBin.ETS)))
d <- rbind(d, with(nalm6.runx1, data.frame("CL" = "NALM6 RUNX1", "Motif" = "EBF", "PeakScore" = `Peak Score`, "MotifScore" = MotifScore.EBF, "MotifScoreBin" = MotifScoreBin.EBF)))
d$CL <- factor(d$CL, levels=c("AT2", "REH", "NALM6 ER", "NALM6 RUNX1"))
d$Motif <- factor(d$Motif, levels=c("RUNX", "ETS", "EBF"))

# ggplot

pdf("/mnt/projects/fiona/results/motif-sequence-analysis.pdf")

print(ggplot(data = d[!is.na(d$MotifScoreBin),], aes(x = MotifScoreBin, y = log(PeakScore, 2))) +   
  geom_violin(fill="#FFCDFF") +
  geom_boxplot(width=.15, size=0.3, fill="white", outlier.colour = NA, show.legend = F) +
  scale_fill_manual(values = rep("white", length(levels(d$CL)))) +
  geom_smooth(aes(group=1), method="gam", size=0.5, col="red") +
  facet_grid(Motif~CL) +
  labs(x="Motif score", y="Log2 peak score") + 
  ggtitle("Peak score by motif score") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=10),
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.3),
        axis.title.y = element_text(size=15, vjust=1.9)) +
  coord_fixed(ratio=0.5))

# most frequent motifs
library(Biostrings)

at2.motifs$Sequence.unstranded <- ifelse(at2.motifs$Strand=="+", as.character(reverseComplement(DNAStringSet(at2.motifs$Sequence))), at2.motifs$Sequence)
at2.freq <- sort(prop.table(table(at2.motifs$Sequence.unstranded[at2.motifs$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer"])), decreasing = T)

reh.motifs$Sequence.unstranded <- ifelse(reh.motifs$Strand=="+", as.character(reverseComplement(DNAStringSet(reh.motifs$Sequence))), reh.motifs$Sequence)
reh.freq <- sort(prop.table(table(reh.motifs$Sequence.unstranded[reh.motifs$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer"])), decreasing = T)

nalm6.er.motifs$Sequence.unstranded <- ifelse(nalm6.er.motifs$Strand=="+", as.character(reverseComplement(DNAStringSet(nalm6.er.motifs$Sequence))), nalm6.er.motifs$Sequence)
nalm6.er.freq <- sort(prop.table(table(nalm6.er.motifs$Sequence.unstranded[nalm6.er.motifs$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer"])), decreasing = T)

nalm6.runx1.motifs$Sequence.unstranded <- ifelse(nalm6.runx1.motifs$Strand=="+", as.character(reverseComplement(DNAStringSet(nalm6.runx1.motifs$Sequence))), nalm6.runx1.motifs$Sequence)
nalm6.runx1.freq <- sort(prop.table(table(nalm6.runx1.motifs$Sequence.unstranded[nalm6.runx1.motifs$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer"])), decreasing = T)

m <- merge(at2.freq, reh.freq, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)+0.003), ylim=c(0, max(m$x, m$y)+0.003), xlab = "kmer frequencyAT2", ylab = "kmer proportion REH", main="RUNX Motif Frequency AT2 vs. REH")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

m <- merge(at2.freq, nalm6.er.freq, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)+0.003), ylim=c(0, max(m$x, m$y)+0.003), xlab = "kmer proportion AT2", ylab = "kmer proportion NALM6 ER", main="RUNX Motif Frequency AT2 vs. NALM6 ER")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

m <- merge(reh.freq, nalm6.er.freq, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)+0.003), ylim=c(0, max(m$x, m$y)+0.003), xlab = "kmer proportion REH", ylab = "kmer proportion NALM6 ER", main="RUNX Motif Frequency REH vs. NALM6 ER")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

m <- merge(at2.freq, nalm6.runx1.freq, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)+0.003), ylim=c(0, max(m$x, m$y)+0.003), xlab = "kmer proportion AT2", ylab = "kmer proportion NALM6 RUNX1", main="RUNX Motif Frequency AT2 vs. NALM6 RUNX1")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

m <- merge(reh.freq, nalm6.runx1.freq, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)+0.003), ylim=c(0, max(m$x, m$y)+0.003), xlab = "kmer proportion REH", ylab = "kmer proportion NALM6 RUNX1", main="RUNX Motif Frequency REH vs. NALM6 RUNX1")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

m <- merge(nalm6.er.freq, nalm6.runx1.freq, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)+0.003), ylim=c(0, max(m$x, m$y)+0.003), xlab = "kmer proportion NALM6 ER", ylab = "kmer proportion NALM6 RUNX1", main="RUNX Motif Frequency NALM6 ER vs. NALM6 RUNX1")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

# AT2 E/R constitutive vs. de novo

at2.motifs.denovo <- at2.motifs[at2.motifs$PositionID %in% at2[at2$runx1_overlap == "de novo", 1],]
at2.motifs.denovo <- at2.motifs.denovo[order(at2.motifs.denovo$PositionID, at2.motifs.denovo$`Motif Name`, at2.motifs.denovo$MotifScore, -abs(at2.motifs.denovo$Offset), decreasing = T),]
at2.freq.denovo <- sort(prop.table(table(at2.motifs.denovo$Sequence.unstranded[at2.motifs.denovo$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer"])), decreasing = T)

at2.motifs.constitutive <- at2.motifs[at2.motifs$PositionID %in% at2[at2$runx1_overlap %in% c("constitutive_better", "constitutive_worse"), 1],]
at2.motifs.constitutive <- at2.motifs.constitutive[order(at2.motifs.constitutive$PositionID, at2.motifs.constitutive$`Motif Name`, at2.motifs.constitutive$MotifScore, -abs(at2.motifs.constitutive$Offset), decreasing = T),]
at2.freq.constitutive <- sort(prop.table(table(at2.motifs.constitutive$Sequence.unstranded[at2.motifs.constitutive$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer"])), decreasing = T)

m <- merge(at2.freq.denovo, at2.freq.constitutive, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)+0.003), ylim=c(0, max(m$x, m$y)+0.003), xlab = "kmer proportion AT2 E/R de novo", ylab = "kmer proportion AT2 E/R constitutive", main="RUNX Motif Frequency AT2 E/R de novo vs. constitutive")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

# REH E/R constitutive vs. de novo

reh.motifs.denovo <- reh.motifs[reh.motifs$PositionID %in% reh[reh$runx1_overlap == "de novo", 1],]
reh.motifs.denovo <- reh.motifs.denovo[order(reh.motifs.denovo$PositionID, reh.motifs.denovo$`Motif Name`, reh.motifs.denovo$MotifScore, -abs(reh.motifs.denovo$Offset), decreasing = T),]
reh.freq.denovo <- sort(prop.table(table(reh.motifs.denovo$Sequence.unstranded[reh.motifs.denovo$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer"])), decreasing = T)

reh.motifs.constitutive <- reh.motifs[reh.motifs$PositionID %in% reh[reh$runx1_overlap %in% c("constitutive_better", "constitutive_worse"), 1],]
reh.motifs.constitutive <- reh.motifs.constitutive[order(reh.motifs.constitutive$PositionID, reh.motifs.constitutive$`Motif Name`, reh.motifs.constitutive$MotifScore, -abs(reh.motifs.constitutive$Offset), decreasing = T),]
reh.freq.constitutive <- sort(prop.table(table(reh.motifs.constitutive$Sequence.unstranded[reh.motifs.constitutive$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer"])), decreasing = T)

m <- merge(reh.freq.denovo, reh.freq.constitutive, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)+0.003), ylim=c(0, max(m$x, m$y)+0.003), xlab = "kmer proportion REH E/R de novo", ylab = "kmer proportion REH E/R constitutive", main="RUNX Motif Frequency REH E/R de novo vs. constitutive")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

# NALM6 E/R constitutive vs. de novo

nalm6.er.motifs.denovo <- nalm6.er.motifs[nalm6.er.motifs$PositionID %in% nalm6.er[nalm6.er$runx1_overlap == "de novo", 1],]
nalm6.er.motifs.denovo <- nalm6.er.motifs.denovo[order(nalm6.er.motifs.denovo$PositionID, nalm6.er.motifs.denovo$`Motif Name`, nalm6.er.motifs.denovo$MotifScore, -abs(nalm6.er.motifs.denovo$Offset), decreasing = T),]
nalm6.er.freq.denovo <- sort(prop.table(table(nalm6.er.motifs.denovo$Sequence.unstranded[nalm6.er.motifs.denovo$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer"])), decreasing = T)

nalm6.er.motifs.constitutive <- nalm6.er.motifs[nalm6.er.motifs$PositionID %in% nalm6.er[nalm6.er$runx1_overlap %in% c("constitutive_better", "constitutive_worse"), 1],]
nalm6.er.motifs.constitutive <- nalm6.er.motifs.constitutive[order(nalm6.er.motifs.constitutive$PositionID, nalm6.er.motifs.constitutive$`Motif Name`, nalm6.er.motifs.constitutive$MotifScore, -abs(nalm6.er.motifs.constitutive$Offset), decreasing = T),]
nalm6.er.freq.constitutive <- sort(prop.table(table(nalm6.er.motifs.constitutive$Sequence.unstranded[nalm6.er.motifs.constitutive$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer"])), decreasing = T)

m <- merge(nalm6.er.freq.denovo, nalm6.er.freq.constitutive, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)+0.003), ylim=c(0, max(m$x, m$y)+0.003), xlab = "kmer proportion NALM6 E/R de novo", ylab = "kmer proportion NALM6 E/R constitutive", main="RUNX Motif Frequency NALM6 E/R de novo vs. constitutive")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

dev.off()


at2.freq <- sort(prop.table(table(at2.motifs$Sequence.unstranded[at2.motifs$`Motif Name`=="EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer"])), decreasing = T)
nalm6.runx1.freq <- sort(prop.table(table(nalm6.runx1.motifs$Sequence.unstranded[nalm6.runx1.motifs$`Motif Name`=="EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer"])), decreasing = T)

m <- merge(at2.freq, nalm6.runx1.freq, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)+0.003), ylim=c(0, max(m$x, m$y)+0.003), xlab = "kmer proportion AT2", ylab = "kmer proportion NALM6 RUNX1", main="RUNX Motif Frequency AT2 vs. NALM6 RUNX1")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)
