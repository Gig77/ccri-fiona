minscore <- 2
maxwidth <- Inf
max.dist.us <- Inf  # note: restricting analysis to promoter peaks (-5000/+2000 bp from TSS) gives very similar results
max.dist.ds <- Inf

# read data

# AT2

at2 <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
at2 <- at2[!grepl("^GL", at2$Chr),]
at2 <- at2[at2$`Peak Score` >= minscore & at2$End-at2$Start <= maxwidth & at2$`Distance to TSS` >= -max.dist.us & at2$`Distance to TSS` <= max.dist.ds,]

at2.motifs <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_AT2_ER.motif_hits.tsv", stringsAsFactors = F, check.names = F)
at2.motifs <- at2.motifs[order(at2.motifs$PositionID, at2.motifs$`Motif Name`, at2.motifs$MotifScore, -abs(at2.motifs$Offset), decreasing = T),]
at2.motifs.best <- at2.motifs[!duplicated(at2.motifs[,c("PositionID", "Motif Name")]),]

# AT2 RUNX

at2 <- merge(at2, at2.motifs.best[at2.motifs.best$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=2, by.y=1, all.x=T)
names(at2)[names(at2) %in% c("Sequence", "MotifScore")] <- c("Sequence.RUNX", "MotifScore.RUNX")
at2$MotifScoreBin.RUNX <- NA
at2$MotifScoreBin.RUNX[!is.na(at2$MotifScore.RUNX)] <- "intermediate"
at2$MotifScoreBin.RUNX[at2$MotifScore.RUNX>=10] <- "high"
at2$MotifScoreBin.RUNX[at2$MotifScore.RUNX<=7] <- "low"
at2$MotifScoreBin.RUNX <- factor(at2$MotifScoreBin.RUNX, levels=c("low", "intermediate", "high"))

# AT2 ETS

at2 <- merge(at2, at2.motifs.best[at2.motifs.best$`Motif Name`=="ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(at2)[names(at2) %in% c("Sequence", "MotifScore")] <- c("Sequence.ETS", "MotifScore.ETS")
at2$MotifScoreBin.ETS <- NA
at2$MotifScoreBin.ETS[!is.na(at2$MotifScore.ETS)] <- "intermediate"
at2$MotifScoreBin.ETS[at2$MotifScore.ETS>=10] <- "high"
at2$MotifScoreBin.ETS[at2$MotifScore.ETS<=8] <- "low"
at2$MotifScoreBin.ETS <- factor(at2$MotifScoreBin.ETS, levels=c("low", "intermediate", "high"))

# AT2 EBF

at2 <- merge(at2, at2.motifs.best[at2.motifs.best$`Motif Name`=="EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(at2)[names(at2) %in% c("Sequence", "MotifScore")] <- c("Sequence.EBF", "MotifScore.EBF")
at2$MotifScoreBin.EBF <- NA
at2$MotifScoreBin.EBF[!is.na(at2$MotifScore.EBF)] <- "intermediate"
at2$MotifScoreBin.EBF[at2$MotifScore.EBF>=11] <- "high"
at2$MotifScoreBin.EBF[at2$MotifScore.EBF<=9.5] <- "low"
at2$MotifScoreBin.EBF <- factor(at2$MotifScoreBin.EBF, levels=c("low", "intermediate", "high"))

# REH

reh <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_REH_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
reh <- reh[!grepl("^GL", reh$Chr),]
reh <- reh[reh$`Peak Score` >= minscore & reh$End-reh$Start <= maxwidth & reh$`Distance to TSS` >= -max.dist.us & reh$`Distance to TSS` <= max.dist.ds,]

reh.motifs <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_REH_ER.motif_hits.tsv", stringsAsFactors = F, check.names = F)
reh.motifs <- reh.motifs[order(reh.motifs$PositionID, reh.motifs$`Motif Name`, reh.motifs$MotifScore, -abs(reh.motifs$Offset), decreasing = T),]
reh.motifs.best <- reh.motifs[!duplicated(reh.motifs[,c("PositionID", "Motif Name")]),]

# REH RUNX

reh <- merge(reh, reh.motifs.best[reh.motifs.best$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=2, by.y=1, all.x=T)
names(reh)[names(reh) %in% c("Sequence", "MotifScore")] <- c("Sequence.RUNX", "MotifScore.RUNX")
reh$MotifScoreBin.RUNX <- NA
reh$MotifScoreBin.RUNX[!is.na(reh$MotifScore.RUNX)] <- "intermediate"
reh$MotifScoreBin.RUNX[reh$MotifScore.RUNX>=10] <- "high"
reh$MotifScoreBin.RUNX[reh$MotifScore.RUNX<=7] <- "low"
reh$MotifScoreBin.RUNX <- factor(reh$MotifScoreBin.RUNX, levels=c("low", "intermediate", "high"))

# REH ETS

reh <- merge(reh, reh.motifs.best[reh.motifs.best$`Motif Name`=="ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(reh)[names(reh) %in% c("Sequence", "MotifScore")] <- c("Sequence.ETS", "MotifScore.ETS")
reh$MotifScoreBin.ETS <- NA
reh$MotifScoreBin.ETS[!is.na(reh$MotifScore.ETS)] <- "intermediate"
reh$MotifScoreBin.ETS[reh$MotifScore.ETS>=10] <- "high"
reh$MotifScoreBin.ETS[reh$MotifScore.ETS<=8] <- "low"
reh$MotifScoreBin.ETS <- factor(reh$MotifScoreBin.ETS, levels=c("low", "intermediate", "high"))

# REH EBF

reh <- merge(reh, reh.motifs.best[reh.motifs.best$`Motif Name`=="EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(reh)[names(reh) %in% c("Sequence", "MotifScore")] <- c("Sequence.EBF", "MotifScore.EBF")
reh$MotifScoreBin.EBF <- NA
reh$MotifScoreBin.EBF[!is.na(reh$MotifScore.EBF)] <- "intermediate"
reh$MotifScoreBin.EBF[reh$MotifScore.EBF>=11] <- "high"
reh$MotifScoreBin.EBF[reh$MotifScore.EBF<=9.5] <- "low"
reh$MotifScoreBin.EBF <- factor(reh$MotifScoreBin.EBF, levels=c("low", "intermediate", "high"))

# NALM-6 ER

nalm6.er <- read.delim("/mnt/projects/fiona/results/homer/ChIP23_NALM6_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.er <- nalm6.er[!grepl("^GL", nalm6.er$Chr),]
nalm6.er <- nalm6.er[nalm6.er$`Peak Score` >= minscore & nalm6.er$End-nalm6.er$Start <= maxwidth & nalm6.er$`Distance to TSS` >= -max.dist.us & nalm6.er$`Distance to TSS` <= max.dist.ds,]

nalm6.er.motifs <- read.delim("/mnt/projects/fiona/results/motifs/ChIP23_NALM6_ER.motif_hits.tsv", stringsAsFactors = F, check.names = F)
nalm6.er.motifs <- nalm6.er.motifs[order(nalm6.er.motifs$PositionID, nalm6.er.motifs$`Motif Name`, nalm6.er.motifs$MotifScore, -abs(nalm6.er.motifs$Offset), decreasing = T),]
nalm6.er.motifs.best <- nalm6.er.motifs[!duplicated(nalm6.er.motifs[,c("PositionID", "Motif Name")]),]

# NALM-6 ER RUNX

nalm6.er <- merge(nalm6.er, nalm6.er.motifs.best[nalm6.er.motifs.best$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=2, by.y=1, all.x=T)
names(nalm6.er)[names(nalm6.er) %in% c("Sequence", "MotifScore")] <- c("Sequence.RUNX", "MotifScore.RUNX")
nalm6.er$MotifScoreBin.RUNX <- NA
nalm6.er$MotifScoreBin.RUNX[!is.na(nalm6.er$MotifScore.RUNX)] <- "intermediate"
nalm6.er$MotifScoreBin.RUNX[nalm6.er$MotifScore.RUNX>=10] <- "high"
nalm6.er$MotifScoreBin.RUNX[nalm6.er$MotifScore.RUNX<=7] <- "low"
nalm6.er$MotifScoreBin.RUNX <- factor(nalm6.er$MotifScoreBin.RUNX, levels=c("low", "intermediate", "high"))

# NALM-6 ER ETS

nalm6.er <- merge(nalm6.er, nalm6.er.motifs.best[nalm6.er.motifs.best$`Motif Name`=="ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(nalm6.er)[names(nalm6.er) %in% c("Sequence", "MotifScore")] <- c("Sequence.ETS", "MotifScore.ETS")
nalm6.er$MotifScoreBin.ETS <- NA
nalm6.er$MotifScoreBin.ETS[!is.na(nalm6.er$MotifScore.ETS)] <- "intermediate"
nalm6.er$MotifScoreBin.ETS[nalm6.er$MotifScore.ETS>=10] <- "high"
nalm6.er$MotifScoreBin.ETS[nalm6.er$MotifScore.ETS<=8] <- "low"
nalm6.er$MotifScoreBin.ETS <- factor(nalm6.er$MotifScoreBin.ETS, levels=c("low", "intermediate", "high"))

# NALM-6 ER EBF

nalm6.er <- merge(nalm6.er, nalm6.er.motifs.best[nalm6.er.motifs.best$`Motif Name`=="EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(nalm6.er)[names(nalm6.er) %in% c("Sequence", "MotifScore")] <- c("Sequence.EBF", "MotifScore.EBF")
nalm6.er$MotifScoreBin.EBF <- NA
nalm6.er$MotifScoreBin.EBF[!is.na(nalm6.er$MotifScore.EBF)] <- "intermediate"
nalm6.er$MotifScoreBin.EBF[nalm6.er$MotifScore.EBF>=11] <- "high"
nalm6.er$MotifScoreBin.EBF[nalm6.er$MotifScore.EBF<=9.5] <- "low"
nalm6.er$MotifScoreBin.EBF <- factor(nalm6.er$MotifScoreBin.EBF, levels=c("low", "intermediate", "high"))

# NALM-6 RUNX1

nalm6.runx1 <- read.delim("/mnt/projects/fiona/results/homer/ChIP22_NALM6_RUNX1_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.runx1 <- nalm6.runx1[!grepl("^GL", nalm6.runx1$Chr),]
nalm6.runx1 <- nalm6.runx1[nalm6.runx1$`Peak Score` >= minscore & nalm6.runx1$End-nalm6.runx1$Start <= maxwidth & nalm6.runx1$`Distance to TSS` >= -max.dist.us & nalm6.runx1$`Distance to TSS` <= max.dist.ds,]

nalm6.runx1.motifs <- read.delim("/mnt/projects/fiona/results/motifs/ChIP22_NALM6_RUNX1.motif_hits.tsv", stringsAsFactors = F, check.names = F)
nalm6.runx1.motifs <- nalm6.runx1.motifs[order(nalm6.runx1.motifs$PositionID, nalm6.runx1.motifs$`Motif Name`, nalm6.runx1.motifs$MotifScore, -abs(nalm6.runx1.motifs$Offset), decreasing = T),]
nalm6.runx1.motifs.best <- nalm6.runx1.motifs[!duplicated(nalm6.runx1.motifs[,c("PositionID", "Motif Name")]),]

# NALM-6 RUNX1 RUNX

nalm6.runx1 <- merge(nalm6.runx1, nalm6.runx1.motifs.best[nalm6.runx1.motifs.best$`Motif Name`=="RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=2, by.y=1, all.x=T)
names(nalm6.runx1)[names(nalm6.runx1) %in% c("Sequence", "MotifScore")] <- c("Sequence.RUNX", "MotifScore.RUNX")
nalm6.runx1$MotifScoreBin.RUNX <- NA
nalm6.runx1$MotifScoreBin.RUNX[!is.na(nalm6.runx1$MotifScore.RUNX)] <- "intermediate"
nalm6.runx1$MotifScoreBin.RUNX[nalm6.runx1$MotifScore.RUNX>=10] <- "high"
nalm6.runx1$MotifScoreBin.RUNX[nalm6.runx1$MotifScore.RUNX<=7] <- "low"
nalm6.runx1$MotifScoreBin.RUNX <- factor(nalm6.runx1$MotifScoreBin.RUNX, levels=c("low", "intermediate", "high"))

# NALM-6 RUNX1 ETS

nalm6.runx1 <- merge(nalm6.runx1, nalm6.runx1.motifs.best[nalm6.runx1.motifs.best$`Motif Name`=="ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(nalm6.runx1)[names(nalm6.runx1) %in% c("Sequence", "MotifScore")] <- c("Sequence.ETS", "MotifScore.ETS")
nalm6.runx1$MotifScoreBin.ETS <- NA
nalm6.runx1$MotifScoreBin.ETS[!is.na(nalm6.runx1$MotifScore.ETS)] <- "intermediate"
nalm6.runx1$MotifScoreBin.ETS[nalm6.runx1$MotifScore.ETS>=10] <- "high"
nalm6.runx1$MotifScoreBin.ETS[nalm6.runx1$MotifScore.ETS<=8] <- "low"
nalm6.runx1$MotifScoreBin.ETS <- factor(nalm6.runx1$MotifScoreBin.ETS, levels=c("low", "intermediate", "high"))

# NALM-6 RUNX1 EBF

nalm6.runx1 <- merge(nalm6.runx1, nalm6.runx1.motifs.best[nalm6.runx1.motifs.best$`Motif Name`=="EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer", c("PositionID", "Sequence", "MotifScore")], by.x=1, by.y=1, all.x=T)
names(nalm6.runx1)[names(nalm6.runx1) %in% c("Sequence", "MotifScore")] <- c("Sequence.EBF", "MotifScore.EBF")
nalm6.runx1$MotifScoreBin.EBF <- NA
nalm6.runx1$MotifScoreBin.EBF[!is.na(nalm6.runx1$MotifScore.EBF)] <- "intermediate"
nalm6.runx1$MotifScoreBin.EBF[nalm6.runx1$MotifScore.EBF>=11] <- "high"
nalm6.runx1$MotifScoreBin.EBF[nalm6.runx1$MotifScore.EBF<=9.5] <- "low"
nalm6.runx1$MotifScoreBin.EBF <- factor(nalm6.runx1$MotifScoreBin.EBF, levels=c("low", "intermediate", "high"))

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

# ggplot

pdf("/mnt/projects/fiona/results/motif-sequence-analysis.pdf")

ggplot(data = d[!is.na(d$MotifScoreBin),], aes(x = MotifScoreBin, y = log(PeakScore, 2))) +   
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
  coord_fixed(ratio=0.5)

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

pdf("/mnt/projects/fiona/results/motif-sequence-analysis.pdf")

m <- merge(at2.freq, reh.freq, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)), ylim=c(0, max(m$x, m$y)), xlab = "kmer frequencyAT2", ylab = "kmer proportion REH", main="RUNX Motif Frequency AT2 vs. REH")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

m <- merge(at2.freq, nalm6.er.freq, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)), ylim=c(0, max(m$x, m$y)), xlab = "kmer proportion AT2", ylab = "kmer proportion NALM6 ER", main="RUNX Motif Frequency AT2 vs. NALM6 ER")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

m <- merge(reh.freq, nalm6.er.freq, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)), ylim=c(0, max(m$x, m$y)), xlab = "kmer proportion REH", ylab = "kmer proportion NALM6 ER", main="RUNX Motif Frequency REH vs. NALM6 ER")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

m <- merge(at2.freq, nalm6.runx1.freq, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)), ylim=c(0, max(m$x, m$y)), xlab = "kmer proportion AT2", ylab = "kmer proportion NALM6 RUNX1", main="RUNX Motif Frequency AT2 vs. NALM6 RUNX1")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

m <- merge(reh.freq, nalm6.runx1.freq, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)), ylim=c(0, max(m$x, m$y)), xlab = "kmer proportion REH", ylab = "kmer proportion NALM6 RUNX1", main="RUNX Motif Frequency REH vs. NALM6 RUNX1")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

m <- merge(nalm6.er.freq, nalm6.runx1.freq, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)), ylim=c(0, max(m$x, m$y)), xlab = "kmer proportion NALM6 ER", ylab = "kmer proportion NALM6 RUNX1", main="RUNX Motif Frequency NALM6 ER vs. NALM6 RUNX1")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)

dev.off()


at2.freq <- sort(prop.table(table(at2.motifs$Sequence.unstranded[at2.motifs$`Motif Name`=="EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer"])), decreasing = T)
nalm6.runx1.freq <- sort(prop.table(table(nalm6.runx1.motifs$Sequence.unstranded[nalm6.runx1.motifs$`Motif Name`=="EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer"])), decreasing = T)

m <- merge(at2.freq, nalm6.runx1.freq, by="row.names", all = T)
m[is.na(m)] <- 0
plot(m$x, m$y, cex=0.3, xlim=c(0, max(m$x, m$y)), ylim=c(0, max(m$x, m$y)), xlab = "kmer proportion AT2", ylab = "kmer proportion NALM6 RUNX1", main="RUNX Motif Frequency AT2 vs. NALM6 RUNX1")
text(m$x, m$y-0.0007, m$Row.names, cex=0.6)
abline(0, 1, lty=2)
