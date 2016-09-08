a <- read.delim("/mnt/projects/fiona/data/homer-annotations/hg19.basic.annotation", header = F, stringsAsFactors = F)

unique(gsub("(.*?)(--| ).*", "\\1", a$V1))

promoterTSS <- a[grepl("^promoter-TSS", a$V1),-7]
write.table(promoterTSS, "/mnt/projects/fiona/results/hg19.basic.promoterTSS", col.names = F, row.names = F, quote = F, sep="\t")

intron <- a[grepl("^intron", a$V1),-7]
write.table(intron, "/mnt/projects/fiona/results/hg19.basic.intron", col.names = F, row.names = F, quote = F, sep="\t")

intergenic <- a[grepl("^Intergenic", a$V1),-7]
write.table(intergenic, "/mnt/projects/fiona/results/hg19.basic.intergenic", col.names = F, row.names = F, quote = F, sep="\t")

noncoding <- a[grepl("^non-coding", a$V1),-7]
write.table(noncoding, "/mnt/projects/fiona/results/hg19.basic.noncoding", col.names = F, row.names = F, quote = F, sep="\t")

TTS <- a[grepl("^TTS", a$V1),-7]
write.table(TTS, "/mnt/projects/fiona/results/hg19.basic.TTS", col.names = F, row.names = F, quote = F, sep="\t")

exon <- a[grepl("^exon", a$V1),-7]
write.table(exon, "/mnt/projects/fiona/results/hg19.basic.exon", col.names = F, row.names = F, quote = F, sep="\t")

fiveUTR <- a[grepl("^5' UTR", a$V1),-7]
write.table(fiveUTR, "/mnt/projects/fiona/results/hg19.basic.fiveUTR", col.names = F, row.names = F, quote = F, sep="\t")

threeUTR <- a[grepl("^3' UTR", a$V1),-7]
write.table(threeUTR, "/mnt/projects/fiona/results/hg19.basic.threeUTR", col.names = F, row.names = F, quote = F, sep="\t")

