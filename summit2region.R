options(warn=1)
library(optparse)

option_list <- list(
  make_option("--summit-file", type="character", help="Peak summits (MACS output, BED format)"),
  make_option("--size", type="integer", help="Extend region 'size' bp up- and downstream of summit"),
  make_option("--out-file", type="character", help="Output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

#opt <- data.frame('summit-file' = "/mnt/projects/fiona/results/macs/ChIP24_AT2_ER_summits.bed", 'size' = 100, 'out-file'= "/mnt/projects/fiona/results/motifs/ChIP24_AT2_ER_peaks.summit-region.bed.part", stringsAsFactors=F, check.names=F)
stopifnot(!is.null(opt$'summit-file'))
stopifnot(!is.null(opt$'size'))
stopifnot(!is.null(opt$'out-file'))

summit <- read.delim(opt$'summit-file', check.names = F, stringsAsFactors = F, header = F)
names(summit) <- c("chr", "start", "end", "id", "score")

summit$peak_id <- gsub("MACS_summit", "MACS_peak", summit$id)
summit$chr <- gsub("^([0-9+|X|Y])", "chr\\1", summit$chr)
summit$chr[summit$chr=="MT"] <- "chrM"

summit$start <- summit$start - opt$size
summit$end <- summit$end + opt$size

write.table(summit[,c("chr", "start", "end", "peak_id", "score")], file=opt$'out-file', row.names = F, col.names = F, quote = F, sep = "\t")
