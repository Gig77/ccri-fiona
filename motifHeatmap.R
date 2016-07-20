qCutoff <- 0.1
pCutoff <- 1e-10
freqCutoff <- 0.1

at2 <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_AT2_ER/knownResults.txt", stringsAsFactors = F, check.names = F)
at2.unique <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_AT2_ER_unique/knownResults.txt", stringsAsFactors = F, check.names = F)
at2.shared <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_AT2_ER_shared/knownResults.txt", stringsAsFactors = F, check.names = F)
reh <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_REH_ER/knownResults.txt", stringsAsFactors = F, check.names = F)
reh.unique <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_REH_ER_unique/knownResults.txt", stringsAsFactors = F, check.names = F)
reh.shared <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_REH_ER_shared/knownResults.txt", stringsAsFactors = F, check.names = F)
nalm6.er <- read.delim("/mnt/projects/fiona/results/motifs/ChIP23_NALM6_ER/knownResults.txt", stringsAsFactors = F, check.names = F)
nalm6.er.unique <- read.delim("/mnt/projects/fiona/results/motifs/ChIP23_NALM6_ER_unique/knownResults.txt", stringsAsFactors = F, check.names = F)
nalm6.er.shared <- read.delim("/mnt/projects/fiona/results/motifs/ChIP23_NALM6_ER_shared/knownResults.txt", stringsAsFactors = F, check.names = F)
nalm6.runx1 <- read.delim("/mnt/projects/fiona/results/motifs/ChIP22_NALM6_RUNX1/knownResults.txt", stringsAsFactors = F, check.names = F)

names(at2)[names(at2) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.at2", "q.at2", "freq.at2")
at2$freq.at2 <- as.numeric(gsub("%", "", at2$freq.at2))

names(at2.unique)[names(at2.unique) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.at2.unique", "q.at2.unique", "freq.at2.unique")
at2.unique$freq.at2.unique <- as.numeric(gsub("%", "", at2.unique$freq.at2.unique))

names(at2.shared)[names(at2.shared) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.at2.shared", "q.at2.shared", "freq.at2.shared")
at2.shared$freq.at2.shared <- as.numeric(gsub("%", "", at2.shared$freq.at2.shared))

names(reh)[names(reh) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.reh", "q.reh", "freq.reh")
reh$freq.reh <- as.numeric(gsub("%", "", reh$freq.reh))

names(reh.unique)[names(reh.unique) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.reh.unique", "q.reh.unique", "freq.reh.unique")
reh.unique$freq.reh.unique <- as.numeric(gsub("%", "", reh.unique$freq.reh.unique))

names(reh.shared)[names(reh.shared) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.reh.shared", "q.reh.shared", "freq.reh.shared")
reh.shared$freq.reh.shared <- as.numeric(gsub("%", "", reh.shared$freq.reh.shared))

names(nalm6.er)[names(nalm6.er) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.nalm6.er", "q.nalm6.er", "freq.nalm6.er")
nalm6.er$freq.nalm6.er <- as.numeric(gsub("%", "", nalm6.er$freq.nalm6.er))

names(nalm6.er.unique)[names(nalm6.er.unique) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.nalm6.er.unique", "q.nalm6.er.unique", "freq.nalm6.er.unique")
nalm6.er.unique$freq.nalm6.er.unique <- as.numeric(gsub("%", "", nalm6.er.unique$freq.nalm6.er.unique))

names(nalm6.er.shared)[names(nalm6.er.shared) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.nalm6.er.shared", "q.nalm6.er.shared", "freq.nalm6.er.shared")
nalm6.er.shared$freq.nalm6.er.shared <- as.numeric(gsub("%", "", nalm6.er.shared$freq.nalm6.er.shared))

names(nalm6.runx1)[names(nalm6.runx1) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.nalm6.runx1", "q.nalm6.runx1", "freq.nalm6.runx1")
nalm6.runx1$freq.nalm6.runx1 <- as.numeric(gsub("%", "", nalm6.runx1$freq.nalm6.runx1))

combined <- at2[,c("motif", "p.at2", "q.at2", "freq.at2")]
combined <- merge(combined, reh[,c("motif", "p.reh", "q.reh", "freq.reh")], all=T)
combined <- merge(combined, nalm6.er[,c("motif", "p.nalm6.er", "q.nalm6.er", "freq.nalm6.er")], all=T)
combined <- merge(combined, nalm6.runx1[,c("motif", "p.nalm6.runx1", "q.nalm6.runx1", "freq.nalm6.runx1")], all=T)
combined <- merge(combined, at2.unique[,c("motif", "p.at2.unique", "q.at2.unique", "freq.at2.unique")], all=T)
combined <- merge(combined, reh.unique[,c("motif", "p.reh.unique", "q.reh.unique", "freq.reh.unique")], all=T)
combined <- merge(combined, nalm6.er.unique[,c("motif", "p.nalm6.er.unique", "q.nalm6.er.unique", "freq.nalm6.er.unique")], all=T)
combined <- merge(combined, at2.shared[,c("motif", "p.at2.shared", "q.at2.shared", "freq.at2.shared")], all=T)
combined <- merge(combined, reh.shared[,c("motif", "p.reh.shared", "q.reh.shared", "freq.reh.shared")], all=T)
combined <- merge(combined, nalm6.er.shared[,c("motif", "p.nalm6.er.shared", "q.nalm6.er.shared", "freq.nalm6.er.shared")], all=T)
rownames(combined) <- combined$motif

combined$minP <- apply(combined[,grepl("^p.", names(combined))], 1, min)
combined$minQ <- apply(combined[,grepl("^q.", names(combined))], 1, min)
combined$maxFreq <- apply(combined[,grepl("^freq.", names(combined))], 1, max)
combined <- combined[combined$minQ <= qCutoff & combined$minP <= pCutoff & combined$maxFreq >= freqCutoff*100,]

d.freq <- as.matrix(combined[,c("freq.at2", "freq.reh", "freq.nalm6.er", "freq.nalm6.runx1", "freq.at2.unique", "freq.reh.unique", "freq.nalm6.er.unique", "freq.at2.shared", "freq.reh.shared", "freq.nalm6.er.shared")])
rownames(d.freq) <- combined$motif
colnames(d.freq) <- c("AT2", "REH", "NALM6 ER", "NALM6 RUNX1", "AT2 unique", "REH unique", "NALM6 ER unique", "AT2 shared", "REH shared", "NALM6 ER shared")

d.p <- as.matrix(combined[,c("p.at2", "p.reh", "p.nalm6.er", "p.nalm6.runx1", "p.at2.unique", "p.reh.unique", "p.nalm6.er.unique", "p.at2.shared", "p.reh.shared", "p.nalm6.er.shared")])
colnames(d.p) <- c("AT2", "REH", "NALM6 ER", "NALM6 RUNX1", "AT2 unique", "REH unique", "NALM6 ER unique", "AT2 shared", "REH shared", "NALM6 ER shared")

d.q <- as.matrix(combined[,c("q.at2", "q.reh", "q.nalm6.er", "q.nalm6.runx1", "q.at2.unique", "q.reh.unique", "q.nalm6.er.unique", "q.at2.shared", "q.reh.shared", "q.nalm6.er.shared")])
colnames(d.q) <- c("AT2", "REH", "NALM6 ER", "NALM6 RUNX1", "AT2 unique", "REH unique", "NALM6 ER unique", "AT2 shared", "REH shared", "NALM6 ER shared")

cut.onestar <- 1e-100
cut.twostars <- 1e-200
#d.stars <- d.p ; d.stars[,] <- NA
#d.stars[d.p <= cut.onestar & d.p > cut.twostars] <- "*"
#d.stars[d.p <= cut.twostars] <- "**"
d.stars <- matrix(sprintf("%.1g", d.p), ncol=ncol(d.p))
d.stars[as.numeric(d.stars) > pCutoff] <- ""
# heatmap

library("RColorBrewer")
library("gplots")

pdf("/mnt/projects/fiona/results/motifHeatmap.pdf", height=12, width=10)
distfun <- function(x, ...) dist(x)
hclustfun <- function(x, ...) hclust(x, method="average", ...)
heatmap.2(
  as.matrix(d.freq), 
  distfun = distfun,
  hclustfun = hclustfun,
  Colv=F, 
  Rowv=T, 
  dendrogram="row", 
  trace="none", 
  col=rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)), 
  margin=c(7, 25), 
  cexCol=0.8, 
  cexRow=0.6, 
  notecex = 0.5,
  keysize=1, 
  colsep=seq(1:ncol(d.freq)), 
  rowsep=seq(1:nrow(d.freq)), 
  sepcolor="grey92", 
  sepwidth=c(0.005,0.005),
  cellnote=d.stars, 
  notecol='white',
  main=sprintf("Enriched motifs\n(* p <= %.2g, ** p <= %.2g)", cut.onestar, cut.twostars),
  key.title="Motif frequency (%)", 
  key.xlab="", 
  key.ylab="",
  scale="none"
)
dev.off()
