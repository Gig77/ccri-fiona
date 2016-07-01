library(GenomicRanges)

minscore <- 2
maxwidth <- Inf
max.dist.us <- Inf  # note: restricting analysis to promoter peaks (-5000/+2000 bp from TSS) gives very similar results
max.dist.ds <- Inf
minoverlap <- 100

at2 <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
at2 <- at2[!grepl("^GL", at2$Chr),]
at2 <- with(at2, at2[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])
at2.gr <- GRanges(seqnames = at2$Chr, ranges=IRanges(at2$Start, at2$End), mcols=data.frame(score=at2$`Peak Score`))

at2.shuffled <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_shuffled_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
at2.shuffled <- at2.shuffled[!grepl("^GL", at2.shuffled$Chr),]
at2.shuffled <- with(at2.shuffled, at2.shuffled[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])

reh <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_REH_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
reh <- reh[!grepl("^GL", reh$Chr),]
reh <- with(reh, reh[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])
reh.gr <- GRanges(seqnames = reh$Chr, ranges=IRanges(reh$Start, reh$End), mcols=data.frame(score=reh$`Peak Score`))

reh.shuffled <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_REH_ER_shuffled_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
reh.shuffled <- reh.shuffled[!grepl("^GL", reh.shuffled$Chr),]
reh.shuffled <- with(reh.shuffled, reh.shuffled[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])

nalm6.runx1 <- read.delim("/mnt/projects/fiona/results/homer/ChIP22_NALM6_RUNX1_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.runx1 <- nalm6.runx1[!grepl("^GL", nalm6.runx1$Chr),]
nalm6.runx1 <- with(nalm6.runx1, nalm6.runx1[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])
nalm6.runx1.gr <- GRanges(seqnames = nalm6.runx1$Chr, ranges=IRanges(nalm6.runx1$Start, nalm6.runx1$End), mcols=data.frame(score=nalm6.runx1$`Peak Score`))

nalm6.er <- read.delim("/mnt/projects/fiona/results/homer/ChIP23_NALM6_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.er <- nalm6.er[!grepl("^GL", nalm6.er$Chr),]
nalm6.er <- with(nalm6.er, nalm6.er[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])
nalm6.er.gr <- GRanges(seqnames = nalm6.er$Chr, ranges=IRanges(nalm6.er$Start, nalm6.er$End), mcols=data.frame(score=nalm6.er$`Peak Score`))

nalm6.runx1.shuffled <- read.delim("/mnt/projects/fiona/results/homer/ChIP22_NALM6_RUNX1_shuffled_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.runx1.shuffled <- nalm6.runx1.shuffled[!grepl("^GL", nalm6.runx1.shuffled$Chr),]
nalm6.runx1.shuffled <- with(nalm6.runx1.shuffled, nalm6.runx1.shuffled[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])

nalm6.er.shuffled <- read.delim("/mnt/projects/fiona/results/homer/ChIP23_NALM6_ER_shuffled_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.er.shuffled <- nalm6.er.shuffled[!grepl("^GL", nalm6.er.shuffled$Chr),]
nalm6.er.shuffled <- with(nalm6.er.shuffled, nalm6.er.shuffled[`Peak Score` >= minscore & End-Start <= maxwidth & `Distance to TSS` >= -max.dist.us & `Distance to TSS` <= max.dist.ds,])

# compute overlaps

o.at2.reh <- findOverlaps(at2.gr, reh.gr, minoverlap=minoverlap)
o.nalm6.runx1.er <- findOverlaps(nalm6.runx1.gr, nalm6.er.gr, minoverlap=minoverlap)
o.nalm6er.at2 <- findOverlaps(nalm6.er.gr, at2.gr, minoverlap=minoverlap)
o.nalm6er.reh <- findOverlaps(nalm6.er.gr, reh.gr, minoverlap=minoverlap)
