qCutoff <- 0.1
pCutoff <- 1e-10
freqCutoff <- 0.1

at2 <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_AT2_ER/knownResults.txt", stringsAsFactors = F, check.names = F)
at2.denovo <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_AT2_ER_denovo/knownResults.txt", stringsAsFactors = F, check.names = F)
at2.const <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_AT2_ER_constitutive/knownResults.txt", stringsAsFactors = F, check.names = F)
reh <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_REH_ER/knownResults.txt", stringsAsFactors = F, check.names = F)
reh.denovo <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_REH_ER_denovo/knownResults.txt", stringsAsFactors = F, check.names = F)
reh.const <- read.delim("/mnt/projects/fiona/results/motifs/ChIP24_REH_ER_constitutive/knownResults.txt", stringsAsFactors = F, check.names = F)
nalm6.er <- read.delim("/mnt/projects/fiona/results/motifs/ChIP23_NALM6_ER/knownResults.txt", stringsAsFactors = F, check.names = F)
nalm6.er.denovo <- read.delim("/mnt/projects/fiona/results/motifs/ChIP23_NALM6_ER_denovo/knownResults.txt", stringsAsFactors = F, check.names = F)
nalm6.er.const <- read.delim("/mnt/projects/fiona/results/motifs/ChIP23_NALM6_ER_constitutive/knownResults.txt", stringsAsFactors = F, check.names = F)
nalm6.runx1 <- read.delim("/mnt/projects/fiona/results/motifs/ChIP22_NALM6_RUNX1/knownResults.txt", stringsAsFactors = F, check.names = F)

names(at2)[names(at2) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.at2", "q.at2", "freq.at2")
at2$freq.at2 <- as.numeric(gsub("%", "", at2$freq.at2))

names(at2.denovo)[names(at2.denovo) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.at2.denovo", "q.at2.denovo", "freq.at2.denovo")
at2.denovo$freq.at2.denovo <- as.numeric(gsub("%", "", at2.denovo$freq.at2.denovo))

names(at2.const)[names(at2.const) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.at2.const", "q.at2.const", "freq.at2.const")
at2.const$freq.at2.const <- as.numeric(gsub("%", "", at2.const$freq.at2.const))

names(reh)[names(reh) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.reh", "q.reh", "freq.reh")
reh$freq.reh <- as.numeric(gsub("%", "", reh$freq.reh))

names(reh.denovo)[names(reh.denovo) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.reh.denovo", "q.reh.denovo", "freq.reh.denovo")
reh.denovo$freq.reh.denovo <- as.numeric(gsub("%", "", reh.denovo$freq.reh.denovo))

names(reh.const)[names(reh.const) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.reh.const", "q.reh.const", "freq.reh.const")
reh.const$freq.reh.const <- as.numeric(gsub("%", "", reh.const$freq.reh.const))

names(nalm6.er)[names(nalm6.er) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.nalm6.er", "q.nalm6.er", "freq.nalm6.er")
nalm6.er$freq.nalm6.er <- as.numeric(gsub("%", "", nalm6.er$freq.nalm6.er))

names(nalm6.er.denovo)[names(nalm6.er.denovo) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.nalm6.er.denovo", "q.nalm6.er.denovo", "freq.nalm6.er.denovo")
nalm6.er.denovo$freq.nalm6.er.denovo <- as.numeric(gsub("%", "", nalm6.er.denovo$freq.nalm6.er.denovo))

names(nalm6.er.const)[names(nalm6.er.const) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.nalm6.er.const", "q.nalm6.er.const", "freq.nalm6.er.const")
nalm6.er.const$freq.nalm6.er.const <- as.numeric(gsub("%", "", nalm6.er.const$freq.nalm6.er.const))

names(nalm6.runx1)[names(nalm6.runx1) %in% c("Motif Name", "P-value", "q-value (Benjamini)", "% of Target Sequences with Motif")] <- c("motif", "p.nalm6.runx1", "q.nalm6.runx1", "freq.nalm6.runx1")
nalm6.runx1$freq.nalm6.runx1 <- as.numeric(gsub("%", "", nalm6.runx1$freq.nalm6.runx1))

combined <- at2[,c("motif", "p.at2", "q.at2", "freq.at2")]
combined <- merge(combined, reh[,c("motif", "p.reh", "q.reh", "freq.reh")], all=T)
combined <- merge(combined, nalm6.er[,c("motif", "p.nalm6.er", "q.nalm6.er", "freq.nalm6.er")], all=T)
combined <- merge(combined, nalm6.runx1[,c("motif", "p.nalm6.runx1", "q.nalm6.runx1", "freq.nalm6.runx1")], all=T)
combined <- merge(combined, at2.denovo[,c("motif", "p.at2.denovo", "q.at2.denovo", "freq.at2.denovo")], all=T)
combined <- merge(combined, reh.denovo[,c("motif", "p.reh.denovo", "q.reh.denovo", "freq.reh.denovo")], all=T)
combined <- merge(combined, nalm6.er.denovo[,c("motif", "p.nalm6.er.denovo", "q.nalm6.er.denovo", "freq.nalm6.er.denovo")], all=T)
combined <- merge(combined, at2.const[,c("motif", "p.at2.const", "q.at2.const", "freq.at2.const")], all=T)
combined <- merge(combined, reh.const[,c("motif", "p.reh.const", "q.reh.const", "freq.reh.const")], all=T)
combined <- merge(combined, nalm6.er.const[,c("motif", "p.nalm6.er.const", "q.nalm6.er.const", "freq.nalm6.er.const")], all=T)
rownames(combined) <- combined$motif

combined$minP <- apply(combined[,grepl("^p.", names(combined))], 1, min)
combined$minQ <- apply(combined[,grepl("^q.", names(combined))], 1, min)
combined$maxFreq <- apply(combined[,grepl("^freq.", names(combined))], 1, max)
combined <- combined[combined$minQ <= qCutoff & combined$minP <= pCutoff & combined$maxFreq >= freqCutoff*100,]

d.freq <- as.matrix(combined[,c("freq.at2", "freq.reh", "freq.nalm6.er", "freq.nalm6.runx1", "freq.at2.denovo", "freq.reh.denovo", "freq.nalm6.er.denovo", "freq.at2.const", "freq.reh.const", "freq.nalm6.er.const")])
colnames(d.freq) <- c("AT2", "REH", "NALM6 ER", "NALM6 RUNX1", "AT2 denovo", "REH denovo", "NALM6 ER denovo", "AT2 const", "REH const", "NALM6 ER const")

d.p <- as.matrix(combined[,c("p.at2", "p.reh", "p.nalm6.er", "p.nalm6.runx1", "p.at2.denovo", "p.reh.denovo", "p.nalm6.er.denovo", "p.at2.const", "p.reh.const", "p.nalm6.er.const")])
colnames(d.p) <- c("AT2", "REH", "NALM6 ER", "NALM6 RUNX1", "AT2 denovo", "REH denovo", "NALM6 ER denovo", "AT2 const", "REH const", "NALM6 ER const")

d.q <- as.matrix(combined[,c("q.at2", "q.reh", "q.nalm6.er", "q.nalm6.runx1", "q.at2.denovo", "q.reh.denovo", "q.nalm6.er.denovo", "q.at2.const", "q.reh.const", "q.nalm6.er.const")])
colnames(d.q) <- c("AT2", "REH", "NALM6 ER", "NALM6 RUNX1", "AT2 denovo", "REH denovo", "NALM6 ER denovo", "AT2 const", "REH const", "NALM6 ER const")

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
  cexRow=0.8, 
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
