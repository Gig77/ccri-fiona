library(GenomicRanges)

minscore <- 2
maxwidth <- Inf
max.dist.us <- Inf  # note: restricting analysis to promoter peaks (-5000/+2000 bp from TSS) gives very similar results
max.dist.ds <- Inf
minoverlap <- 100
max.motif.distance.from.summit <- 100

hasMotifNearSummit <- function(motifDistances) {
  unlist(sapply(motifDistances, function(x) {   
    if (x != "") {
      min(abs(as.numeric(unlist(strsplit(x, ","))))) <= max.motif.distance.from.summit
    } else {
      FALSE
    } }))
}

numMotifsNearSummit <- function(motifDistances) {
  unlist(sapply(motifDistances, function(x) {   
    if (x != "") {
      sum(abs(as.numeric(unlist(strsplit(x, ",")))) <= max.motif.distance.from.summit)
    } else {
      0
    } }))
}

at2 <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
at2 <- at2[!grepl("^GL", at2$Chr),]
at2 <- with(at2, at2[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])
at2$RunxNearSummit <- hasMotifNearSummit(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
at2$EtsNearSummit <- hasMotifNearSummit(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
at2$EbfNearSummit <- hasMotifNearSummit(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)
at2$NumRunxNearSummit <- numMotifsNearSummit(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
at2$NumEtsNearSummit <- numMotifsNearSummit(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
at2$NumEbfNearSummit <- numMotifsNearSummit(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)
at2.gr <- GRanges(seqnames = at2$Chr, ranges=IRanges(at2$Start, at2$End), mcols=data.frame(score=at2$`Peak Score`, RunxNearSummit=at2$RunxNearSummit))

at2.better <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_better_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
at2.better <- at2.better[!grepl("^GL", at2.better$Chr),]
at2.better <- with(at2.better, at2.better[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])
at2.better$RunxNearSummit <- hasMotifNearSummit(at2.better$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
at2.better$EtsNearSummit <- hasMotifNearSummit(at2.better$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
at2.better$EbfNearSummit <- hasMotifNearSummit(at2.better$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)
at2.better$NumRunxNearSummit <- numMotifsNearSummit(at2.better$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
at2.better$NumEtsNearSummit <- numMotifsNearSummit(at2.better$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
at2.better$NumEbfNearSummit <- numMotifsNearSummit(at2.better$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)

at2.worse <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_worse_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
at2.worse <- at2.worse[!grepl("^GL", at2.worse$Chr),]
at2.worse <- with(at2.worse, at2.worse[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])
at2.worse$RunxNearSummit <- hasMotifNearSummit(at2.worse$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
at2.worse$EtsNearSummit <- hasMotifNearSummit(at2.worse$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
at2.worse$EbfNearSummit <- hasMotifNearSummit(at2.worse$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)
at2.worse$NumRunxNearSummit <- numMotifsNearSummit(at2.worse$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
at2.worse$NumEtsNearSummit <- numMotifsNearSummit(at2.worse$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
at2.worse$NumEbfNearSummit <- numMotifsNearSummit(at2.worse$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)

at2.shuffled <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_shuffled_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
at2.shuffled <- at2.shuffled[!grepl("^GL", at2.shuffled$Chr),]
at2.shuffled <- with(at2.shuffled, at2.shuffled[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])
at2.shuffled$RunxNearSummit <- hasMotifNearSummit(at2.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
at2.shuffled$EtsNearSummit <- hasMotifNearSummit(at2.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
at2.shuffled$EbfNearSummit <- hasMotifNearSummit(at2.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)
at2.shuffled$NumRunxNearSummit <- numMotifsNearSummit(at2.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
at2.shuffled$NumEtsNearSummit <- numMotifsNearSummit(at2.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
at2.shuffled$NumEbfNearSummit <- numMotifsNearSummit(at2.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)

reh <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_REH_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
reh <- reh[!grepl("^GL", reh$Chr),]
reh <- with(reh, reh[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])
reh$RunxNearSummit <- hasMotifNearSummit(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
reh$EtsNearSummit <- hasMotifNearSummit(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
reh$EbfNearSummit <- hasMotifNearSummit(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)
reh$NumRunxNearSummit <- numMotifsNearSummit(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
reh$NumEtsNearSummit <- numMotifsNearSummit(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
reh$NumEbfNearSummit <- numMotifsNearSummit(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)
reh.gr <- GRanges(seqnames = reh$Chr, ranges=IRanges(reh$Start, reh$End), mcols=data.frame(score=reh$`Peak Score`, RunxNearSummit=reh$RunxNearSummit))

reh.shuffled <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_REH_ER_shuffled_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
reh.shuffled <- reh.shuffled[!grepl("^GL", reh.shuffled$Chr),]
reh.shuffled <- with(reh.shuffled, reh.shuffled[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])
reh.shuffled$RunxNearSummit <- hasMotifNearSummit(reh.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
reh.shuffled$EtsNearSummit <- hasMotifNearSummit(reh.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
reh.shuffled$EbfNearSummit <- hasMotifNearSummit(reh.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)
reh.shuffled$NumRunxNearSummit <- numMotifsNearSummit(reh.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
reh.shuffled$NumEtsNearSummit <- numMotifsNearSummit(reh.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
reh.shuffled$NumEbfNearSummit <- numMotifsNearSummit(reh.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)

nalm6.runx1 <- read.delim("/mnt/projects/fiona/results/homer/ChIP22_NALM6_RUNX1_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.runx1 <- nalm6.runx1[!grepl("^GL", nalm6.runx1$Chr),]
nalm6.runx1 <- with(nalm6.runx1, nalm6.runx1[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])
nalm6.runx1$RunxNearSummit <- hasMotifNearSummit(nalm6.runx1$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
nalm6.runx1$EtsNearSummit <- hasMotifNearSummit(nalm6.runx1$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
nalm6.runx1$EbfNearSummit <- hasMotifNearSummit(nalm6.runx1$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)
nalm6.runx1$NumRunxNearSummit <- numMotifsNearSummit(nalm6.runx1$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
nalm6.runx1$NumEtsNearSummit <- numMotifsNearSummit(nalm6.runx1$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
nalm6.runx1$NumEbfNearSummit <- numMotifsNearSummit(nalm6.runx1$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)
nalm6.runx1.gr <- GRanges(seqnames = nalm6.runx1$Chr, ranges=IRanges(nalm6.runx1$Start, nalm6.runx1$End), mcols=data.frame(score=nalm6.runx1$`Peak Score`, RunxNearSummit=nalm6.runx1$RunxNearSummit))

nalm6.er <- read.delim("/mnt/projects/fiona/results/homer/ChIP23_NALM6_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.er <- nalm6.er[!grepl("^GL", nalm6.er$Chr),]
nalm6.er <- with(nalm6.er, nalm6.er[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])
nalm6.er$RunxNearSummit <- hasMotifNearSummit(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
nalm6.er$EtsNearSummit <- hasMotifNearSummit(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
nalm6.er$EbfNearSummit <- hasMotifNearSummit(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)
nalm6.er$NumRunxNearSummit <- numMotifsNearSummit(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
nalm6.er$NumEtsNearSummit <- numMotifsNearSummit(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
nalm6.er$NumEbfNearSummit <- numMotifsNearSummit(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)
nalm6.er.gr <- GRanges(seqnames = nalm6.er$Chr, ranges=IRanges(nalm6.er$Start, nalm6.er$End), mcols=data.frame(score=nalm6.er$`Peak Score`, RunxNearSummit=nalm6.er$RunxNearSummit))

nalm6.runx1.shuffled <- read.delim("/mnt/projects/fiona/results/homer/ChIP22_NALM6_RUNX1_shuffled_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.runx1.shuffled <- nalm6.runx1.shuffled[!grepl("^GL", nalm6.runx1.shuffled$Chr),]
nalm6.runx1.shuffled <- with(nalm6.runx1.shuffled, nalm6.runx1.shuffled[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])
nalm6.runx1.shuffled$RunxNearSummit <- hasMotifNearSummit(nalm6.runx1.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
nalm6.runx1.shuffled$EtsNearSummit <- hasMotifNearSummit(nalm6.runx1.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
nalm6.runx1.shuffled$EbfNearSummit <- hasMotifNearSummit(nalm6.runx1.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)
nalm6.runx1.shuffled$NumRunxNearSummit <- numMotifsNearSummit(nalm6.runx1.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
nalm6.runx1.shuffled$NumEtsNearSummit <- numMotifsNearSummit(nalm6.runx1.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
nalm6.runx1.shuffled$NumEbfNearSummit <- numMotifsNearSummit(nalm6.runx1.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)

nalm6.er.shuffled <- read.delim("/mnt/projects/fiona/results/homer/ChIP23_NALM6_ER_shuffled_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.er.shuffled <- nalm6.er.shuffled[!grepl("^GL", nalm6.er.shuffled$Chr),]
nalm6.er.shuffled <- with(nalm6.er.shuffled, nalm6.er.shuffled[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])
nalm6.er.shuffled$RunxNearSummit <- hasMotifNearSummit(nalm6.er.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
nalm6.er.shuffled$EtsNearSummit <- hasMotifNearSummit(nalm6.er.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
nalm6.er.shuffled$EbfNearSummit <- hasMotifNearSummit(nalm6.er.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)
nalm6.er.shuffled$NumRunxNearSummit <- numMotifsNearSummit(nalm6.er.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`)
nalm6.er.shuffled$NumEtsNearSummit <- numMotifsNearSummit(nalm6.er.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`)
nalm6.er.shuffled$NumEbfNearSummit <- numMotifsNearSummit(nalm6.er.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`)

# compute overlaps

o.at2.reh <- findOverlaps(at2.gr, reh.gr, minoverlap=minoverlap)
o.nalm6.runx1.er <- findOverlaps(nalm6.runx1.gr, nalm6.er.gr, minoverlap=minoverlap)
o.nalm6er.at2 <- findOverlaps(nalm6.er.gr, at2.gr, minoverlap=minoverlap)
o.nalm6er.reh <- findOverlaps(nalm6.er.gr, reh.gr, minoverlap=minoverlap)
o.at2.nalm6.runx1 <- findOverlaps(at2.gr, nalm6.runx1.gr, minoverlap=minoverlap)
o.reh.nalm6.runx1 <- findOverlaps(reh.gr, nalm6.runx1.gr, minoverlap=minoverlap)
