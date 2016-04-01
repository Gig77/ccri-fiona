library(Vennerable)

minscore <- 2
minoverlap <- 100
maxwidth <- Inf
max.dist.us <- Inf  # note: restricting analysis to promoter peaks (-5000/+2000 bp from TSS) gives very similar results
max.dist.ds <- Inf

# read data

at2 <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
at2 <- at2[!grepl("^GL", at2$Chr),]
at2 <- at2[at2$`Peak Score` >= minscore & at2$End-at2$Start <= maxwidth & at2$`Distance to TSS` >= -max.dist.us & at2$`Distance to TSS` <= max.dist.ds,]

reh <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_REH_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
reh <- reh[!grepl("^GL", reh$Chr),]
reh <- reh[reh$`Peak Score` >= minscore & reh$End-reh$Start <= maxwidth & reh$`Distance to TSS` >= -max.dist.us & reh$`Distance to TSS` <= max.dist.ds,]

nalm6.er <- read.delim("/mnt/projects/fiona/results/homer/ChIP23_NALM6_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.er <- nalm6.er[!grepl("^GL", nalm6.er$Chr),]
nalm6.er <- nalm6.er[nalm6.er$`Peak Score` >= minscore & nalm6.er$End-nalm6.er$Start <= maxwidth & nalm6.er$`Distance to TSS` >= -max.dist.us & nalm6.er$`Distance to TSS` <= max.dist.ds,]

nalm6.runx1 <- read.delim("/mnt/projects/fiona/results/homer/ChIP22_NALM6_RUNX1_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.runx1 <- nalm6.runx1[!grepl("^GL", nalm6.runx1$Chr),]
nalm6.runx1 <- nalm6.runx1[nalm6.runx1$`Peak Score` >= minscore & nalm6.runx1$End-nalm6.runx1$Start <= maxwidth & nalm6.runx1$`Distance to TSS` >= -max.dist.us & nalm6.runx1$`Distance to TSS` <= max.dist.ds,]

# peak overlaps

at2.gr <- GRanges(seqnames = at2$Chr, ranges=IRanges(at2$Start, at2$End), mcols=data.frame(score=at2$`Peak Score`))
reh.gr <- GRanges(seqnames = reh$Chr, ranges=IRanges(reh$Start, reh$End), mcols=data.frame(score=reh$`Peak Score`))
nalm6.runx1.gr <- GRanges(seqnames = nalm6.runx1$Chr, ranges=IRanges(nalm6.runx1$Start, nalm6.runx1$End), mcols=data.frame(score=nalm6.runx1$`Peak Score`))
nalm6.er.gr <- GRanges(seqnames = nalm6.er$Chr, ranges=IRanges(nalm6.er$Start, nalm6.er$End), mcols=data.frame(score=nalm6.er$`Peak Score`))

o.at2.nalm6.runx1 <- findOverlaps(at2.gr, nalm6.runx1.gr, minoverlap=minoverlap)
o.reh.nalm6.runx1 <- findOverlaps(reh.gr, nalm6.runx1.gr, minoverlap=minoverlap)
o.nalm6.er.runx1 <- findOverlaps(nalm6.er.gr, nalm6.runx1.gr, minoverlap=minoverlap)


pdf("/mnt/projects/fiona/results/motif-analysis.pdf", paper="a4")

# overview barchart

d <- data.frame("CL" = character(0), "Motif" = character(0), "Pct" = numeric(0), stringsAsFactors = F)
d <- rbind(d, setNames(data.frame("AT2 E/R", "RUNX1", sum(!is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / nrow(at2)), names(d)))
d <- rbind(d, setNames(data.frame("AT2 E/R", "ETS", sum(!is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / nrow(at2)), names(d)))
d <- rbind(d, setNames(data.frame("AT2 E/R", "EBF", sum(!is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / nrow(at2)), names(d)))
d <- rbind(d, setNames(data.frame("REH E/R", "RUNX1", sum(!is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / nrow(reh)), names(d)))
d <- rbind(d, setNames(data.frame("REH E/R", "ETS", sum(!is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / nrow(reh)), names(d)))
d <- rbind(d, setNames(data.frame("REH E/R", "EBF", sum(!is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / nrow(reh)), names(d)))
d <- rbind(d, setNames(data.frame("NALM6 RUNX1", "RUNX1", sum(!is.na(nalm6.runx1$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / nrow(nalm6.runx1)), names(d)))
d <- rbind(d, setNames(data.frame("NALM6 RUNX1", "ETS", sum(!is.na(nalm6.runx1$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / nrow(nalm6.runx1)), names(d)))
d <- rbind(d, setNames(data.frame("NALM6 RUNX1", "EBF", sum(!is.na(nalm6.runx1$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / nrow(nalm6.runx1)), names(d)))
d <- rbind(d, setNames(data.frame("NALM6 E/R", "RUNX1", sum(!is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / nrow(nalm6.er)), names(d)))
d <- rbind(d, setNames(data.frame("NALM6 E/R", "ETS", sum(!is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / nrow(nalm6.er)), names(d)))
d <- rbind(d, setNames(data.frame("NALM6 E/R", "EBF", sum(!is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / nrow(nalm6.er)), names(d)))

ggplot(data = d, aes(x = CL, y = Pct, fill = Motif)) +   
  geom_bar(position = position_dodge(), stat = "identity") +
  labs(fill="Motif", x="Cell line", y="% of peaks with motif") + 
  ggtitle("Percentage of peaks with known motif") +
  theme_bw() +
  scale_fill_manual(values = c("red", "darkgray", "blue")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.9),
        axis.title.y = element_text(size=15, vjust=1.9))

# number of motifs per peak

d <- data.frame("CL" = character(0), "Motif" = character(0), "CountMotifs" = character(0), "CountPeaks" = numeric(0), stringsAsFactors = F)
d <- rbind(d, setNames(data.frame("AT2", "RUNX1", melt(prop.table(table(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "RUNX1", melt(prop.table(table(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6 ER", "RUNX1", melt(prop.table(table(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6 RUNX1", "RUNX1", melt(prop.table(table(nalm6.runx1$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "ETS", melt(prop.table(table(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "ETS", melt(prop.table(table(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6 ER", "ETS", melt(prop.table(table(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6 RUNX1", "ETS", melt(prop.table(table(nalm6.runx1$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "EBF", melt(prop.table(table(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "EBF", melt(prop.table(table(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6 ER", "EBF", melt(prop.table(table(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6 RUNX1", "EBF", melt(prop.table(table(nalm6.runx1$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)))), names(d)))
d <- d[d$CountMotifs <= 5,]
d$CountMotifs <- as.factor(d$CountMotifs)

#pdf("/mnt/projects/fiona/results/test.pdf", paper="a4")
ggplot(data = d, aes(x = CountMotifs, y = CountPeaks, fill = CL)) +
  geom_bar(position = position_dodge(), stat = "identity", width=0.7) +
  facet_wrap(~Motif) +
  coord_fixed(ratio=6) +
  labs(fill="Cell line", x="No. of motifs per peak", y="% peaks") + 
  ggtitle("No. of motifs in peaks") +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_fill_manual(values = c("red", "darkgray", "blue", "orange")) +
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.9),
        axis.title.y = element_text(size=15, vjust=1.9))
#dev.off()

# motif distance to summit

# AT2

#pdf("/mnt/projects/fiona/results/test.pdf", paper="a4")
par(mfrow=c(2,2))

at2.runx1.dist <- at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`
at2.runx1.dist <- unlist(sapply(at2.runx1.dist[at2.runx1.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
at2.ets.dist <- at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`
at2.ets.dist <- unlist(sapply(at2.ets.dist[at2.ets.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
at2.ebf.dist <- at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`
at2.ebf.dist <- unlist(sapply(at2.ebf.dist[at2.ebf.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))

plot(density(at2.runx1.dist, adjust=1), xlim=c(-200,200), col="red", lwd=3, main="AT2 motif distance from summit", xlab="Distance from peak summit (bp)")
lines(density(at2.ets.dist, adjust=1), xlim=c(-200,200), col="darkgray", lwd=3)
lines(density(at2.ebf.dist, adjust=1), xlim=c(-200,200), col="blue", lwd=3)
abline(v=0, lty=2)
legend("topright", c("RUNX1", "ETS", "EBF"), fill=c("red", "darkgray", "blue"))

# REH

reh.runx1.dist <- reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`
reh.runx1.dist <- unlist(sapply(reh.runx1.dist[reh.runx1.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
reh.ets.dist <- reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`
reh.ets.dist <- unlist(sapply(reh.ets.dist[reh.ets.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
reh.ebf.dist <- reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`
reh.ebf.dist <- unlist(sapply(reh.ebf.dist[reh.ebf.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))

plot(density(reh.runx1.dist), xlim=c(-200,200), col="red", lwd=3, main="REH motif distance from summit", xlab="Distance from peak summit (bp)")
lines(density(reh.ets.dist), xlim=c(-200,200), col="darkgray", lwd=3)
lines(density(reh.ebf.dist, adjust=1), xlim=c(-200,200), col="blue", lwd=3)
abline(v=0, lty=2)
legend("topright", c("RUNX1", "ETS", "EBF"), fill=c("red", "darkgray", "blue"))

# NALM6 ER

nalm6.er.runx1.dist <- nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`
nalm6.er.runx1.dist <- unlist(sapply(nalm6.er.runx1.dist[nalm6.er.runx1.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
nalm6.er.ets.dist <- nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`
nalm6.er.ets.dist <- unlist(sapply(nalm6.er.ets.dist[nalm6.er.ets.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
nalm6.er.ebf.dist <- nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`
nalm6.er.ebf.dist <- unlist(sapply(nalm6.er.ebf.dist[nalm6.er.ebf.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))

plot(density(nalm6.er.runx1.dist), xlim=c(-200,200), col="red", lwd=3, main="NALM-6 ER motif distance from summit", xlab="Distance from peak summit (bp)")
lines(density(nalm6.er.ets.dist), xlim=c(-200,200), col="darkgray", lwd=3)
lines(density(nalm6.er.ebf.dist, adjust=1), xlim=c(-200,200), col="blue", lwd=3)
abline(v=0, lty=2)
legend("topright", c("RUNX1", "ETS", "EBF"), fill=c("red", "darkgray", "blue"))

# NALM6 RUNX1

nalm6.runx1.runx1.dist <- nalm6.runx1$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`
nalm6.runx1.runx1.dist <- unlist(sapply(nalm6.runx1.runx1.dist[nalm6.runx1.runx1.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
nalm6.runx1.ets.dist <- nalm6.runx1$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`
nalm6.runx1.ets.dist <- unlist(sapply(nalm6.runx1.ets.dist[nalm6.runx1.ets.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
nalm6.runx1.ebf.dist <- nalm6.runx1$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`
nalm6.runx1.ebf.dist <- unlist(sapply(nalm6.runx1.ebf.dist[nalm6.runx1.ebf.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))

plot(density(nalm6.runx1.runx1.dist), xlim=c(-200,200), col="red", lwd=3, main="NALM-6 RUNX1 motif distance from summit", xlab="Distance from peak summit (bp)")
lines(density(nalm6.runx1.ets.dist), xlim=c(-200,200), col="darkgray", lwd=3)
lines(density(nalm6.runx1.ebf.dist), xlim=c(-200,200), col="blue", lwd=3)
abline(v=0, lty=2)
legend("topright", c("RUNX1", "ETS", "EBF"), fill=c("red", "darkgray", "blue"))

par(mfrow=c(1,1))
#dev.off()

# barchart constitutive vs. de novo peaks

d <- data.frame("CL" = character(0), "Motif" = character(0), "Class" = character(0), "Pct" = numeric(0), stringsAsFactors = F)
d <- rbind(d, setNames(data.frame("AT2", "RUNX1", "constitutive", sum(1:nrow(at2) %in% o.at2.nalm6.runx1@queryHits & !is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / length(unique(o.at2.nalm6.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "ETS", "constitutive", sum(1:nrow(at2) %in% o.at2.nalm6.runx1@queryHits & !is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / length(unique(o.at2.nalm6.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "EBF", "constitutive", sum(1:nrow(at2) %in% o.at2.nalm6.runx1@queryHits & !is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / length(unique(o.at2.nalm6.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "RUNX1", "constitutive", sum(1:nrow(reh) %in% o.reh.nalm6.runx1@queryHits & !is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / length(unique(o.reh.nalm6.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "ETS", "constitutive", sum(1:nrow(reh) %in% o.reh.nalm6.runx1@queryHits & !is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / length(unique(o.reh.nalm6.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "EBF", "constitutive", sum(1:nrow(reh) %in% o.reh.nalm6.runx1@queryHits & !is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / length(unique(o.reh.nalm6.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "RUNX1", "constitutive", sum(1:nrow(nalm6.er) %in% o.nalm6.er.runx1@queryHits & !is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / length(unique(o.nalm6.er.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "ETS", "constitutive", sum(1:nrow(nalm6.er) %in% o.nalm6.er.runx1@queryHits & !is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / length(unique(o.nalm6.er.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "EBF", "constitutive", sum(1:nrow(nalm6.er) %in% o.nalm6.er.runx1@queryHits & !is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / length(unique(o.nalm6.er.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "RUNX1", "de novo", sum(!(1:nrow(at2) %in% o.at2.nalm6.runx1@queryHits) & !is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / (nrow(at2)-length(unique(o.at2.nalm6.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "ETS", "de novo", sum(!(1:nrow(at2) %in% o.at2.nalm6.runx1@queryHits) & !is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / (nrow(at2)-length(unique(o.at2.nalm6.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "EBF", "de novo", sum(!(1:nrow(at2) %in% o.at2.nalm6.runx1@queryHits) & !is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / (nrow(at2)-length(unique(o.at2.nalm6.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "RUNX1", "de novo", sum(!(1:nrow(reh) %in% o.reh.nalm6.runx1@queryHits) & !is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / (nrow(reh)-length(unique(o.reh.nalm6.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "ETS", "de novo", sum(!(1:nrow(reh) %in% o.reh.nalm6.runx1@queryHits) & !is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / (nrow(reh)-length(unique(o.reh.nalm6.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "EBF", "de novo", sum(!(1:nrow(reh) %in% o.reh.nalm6.runx1@queryHits) & !is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / (nrow(reh)-length(unique(o.reh.nalm6.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "RUNX1", "de novo", sum(!(1:nrow(nalm6.er) %in% o.nalm6.er.runx1@queryHits) & !is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / (nrow(nalm6.er)-length(unique(o.nalm6.er.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "ETS", "de novo", sum(!(1:nrow(nalm6.er) %in% o.nalm6.er.runx1@queryHits) & !is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / (nrow(nalm6.er)-length(unique(o.nalm6.er.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "EBF", "de novo", sum(!(1:nrow(nalm6.er) %in% o.nalm6.er.runx1@queryHits) & !is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / (nrow(nalm6.er)-length(unique(o.nalm6.er.runx1@queryHits)))), names(d)))

#pdf("/mnt/projects/fiona/results/test.pdf", paper="a4")
ggplot(data = d, aes(x = CL, y = Pct, fill = Class)) +   
  geom_bar(position = position_dodge(), stat = "identity", width=0.7) +
  facet_wrap(~Motif) +
  labs(fill="E/R peak", x="Cell line", y="% of E/R peaks with motif") + 
  ggtitle("Percentage of E/R peaks with known motif\n") +
  theme_bw() +
  scale_fill_manual(values = c("black", "darkgray")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  coord_fixed(ratio=8) +
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.9),
        axis.title.y = element_text(size=15, vjust=1.9))
#dev.off()

# motif overlap AT2

peak.sets.at2 <- list(
  "RUNX1" = at2[,2][!is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)],
  "ETS" = at2[,2][!is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)],
  "EBF" = at2[,2][!is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)]
)

venn <- compute.Venn(Venn(peak.sets.at2), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("AT2 Motif Overlap", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

# motif overlap REH

peak.sets.reh <- list(
  "RUNX1" = reh[,2][!is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)],
  "ETS" = reh[,2][!is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)],
  "EBF" = reh[,2][!is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)]
)

venn <- compute.Venn(Venn(peak.sets.reh), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("REH Motif Overlap", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

# motif overlap NALM-6 RUNX1

peak.sets.nalm6.runx1 <- list(
  "RUNX1" = nalm6.runx1[,2][!is.na(nalm6.runx1$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)],
  "ETS" = nalm6.runx1[,2][!is.na(nalm6.runx1$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)],
  "EBF" = nalm6.runx1[,2][!is.na(nalm6.runx1$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)]
)

venn <- compute.Venn(Venn(peak.sets.nalm6.runx1), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("NALM-6 RUNX1 Motif Overlap", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

# motif overlap NALM-6 ER

peak.sets.nalm6.er <- list(
  "RUNX1" = nalm6.er[,2][!is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)],
  "ETS" = nalm6.er[,2][!is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)],
  "EBF" = nalm6.er[,2][!is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)]
)

venn <- compute.Venn(Venn(peak.sets.nalm6.er), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("NALM-6 ER Motif Overlap", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

# motif composition constitutive RUNX1 vs. E/R de novo binding sites

#pdf("/mnt/projects/fiona/results/test.pdf", paper="a4")

# AT2

peak.sets.at2.constitutive <- list(
  "RUNX1" = at2[intersect(o.at2.nalm6.runx1@queryHits, which(!is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`))),2],
  "ETS" = at2[intersect(o.at2.nalm6.runx1@queryHits, which(!is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`))),2],
  "EBF" = at2[intersect(o.at2.nalm6.runx1@queryHits, which(!is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`))),2]
)

venn <- compute.Venn(Venn(peak.sets.at2.constitutive), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("AT2 motif overlap constitutive E/R peaks", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

peak.sets.at2.denovo <- list(
  "RUNX1" = at2[intersect((1:nrow(at2))[-o.at2.nalm6.runx1@queryHits], which(!is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`))),2],
  "ETS" = at2[intersect((1:nrow(at2))[-o.at2.nalm6.runx1@queryHits], which(!is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`))),2],
  "EBF" = at2[intersect((1:nrow(at2))[-o.at2.nalm6.runx1@queryHits], which(!is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`))),2]
)

venn <- compute.Venn(Venn(peak.sets.at2.denovo), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("AT2 motif overlap de novo E/R peaks", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

# REH

peak.sets.reh.constitutive <- list(
  "RUNX1" = reh[intersect(o.reh.nalm6.runx1@queryHits, which(!is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`))),2],
  "ETS" = reh[intersect(o.reh.nalm6.runx1@queryHits, which(!is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`))),2],
  "EBF" = reh[intersect(o.reh.nalm6.runx1@queryHits, which(!is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`))),2]
)

venn <- compute.Venn(Venn(peak.sets.reh.constitutive), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("REH motif overlap constitutive E/R peaks", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

peak.sets.reh.denovo <- list(
  "RUNX1" = reh[intersect((1:nrow(reh))[-o.reh.nalm6.runx1@queryHits], which(!is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`))),2],
  "ETS" = reh[intersect((1:nrow(reh))[-o.reh.nalm6.runx1@queryHits], which(!is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`))),2],
  "EBF" = reh[intersect((1:nrow(reh))[-o.reh.nalm6.runx1@queryHits], which(!is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`))),2]
)

venn <- compute.Venn(Venn(peak.sets.reh.denovo), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("REH motif overlap de novo E/R peaks", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

# NALM-6

peak.sets.nalm6.er.constitutive <- list(
  "RUNX1" = nalm6.er[intersect(o.nalm6.er.runx1@queryHits, which(!is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`))),2],
  "ETS" = nalm6.er[intersect(o.nalm6.er.runx1@queryHits, which(!is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`))),2],
  "EBF" = nalm6.er[intersect(o.nalm6.er.runx1@queryHits, which(!is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`))),2]
)

venn <- compute.Venn(Venn(peak.sets.nalm6.er.constitutive), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("NALM-6 motif overlap constitutive E/R peaks", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

peak.sets.nalm6.er.denovo <- list(
  "RUNX1" = nalm6.er[intersect((1:nrow(nalm6.er))[-o.nalm6.er.runx1@queryHits], which(!is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`))),2],
  "ETS" = nalm6.er[intersect((1:nrow(nalm6.er))[-o.nalm6.er.runx1@queryHits], which(!is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`))),2],
  "EBF" = nalm6.er[intersect((1:nrow(nalm6.er))[-o.nalm6.er.runx1@queryHits], which(!is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`))),2]
)

venn <- compute.Venn(Venn(peak.sets.nalm6.er.denovo), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("NALM-6 motif overlap de novo E/R peaks", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

dev.off()

