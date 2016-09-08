options(warn=1)
library(optparse)

option_list <- list(
  make_option("--annotated-peaks", type="character", help="Annotated ChIP-seq peaks from Homer")
)
opt <- parse_args(OptionParser(option_list=option_list))
opt <- data.frame('annotated-peaks' = "/mnt/projects/fiona/results/homer/ChIP24_REH_ER_peaks.annotated.with-expr.tsv", 
                  'peak-out' = "/mnt/projects/fiona/results/macs/ChIP24_REH_ER_runx1Platinum_peaks.bed.part",
                  'summit-out' = "/mnt/projects/fiona/results/macs/ChIP24_REH_ER_runx1Platinum_summits.bed.part",
                  stringsAsFactors=F, check.names=F)
stopifnot(!is.null(opt$'annotated-peaks'))

peaks <- read.delim(opt$'annotated-peaks', check.names = F, stringsAsFactors = F)

# FILTER: peak score
min.score <- 10
peaks <- peaks[peaks$`Peak Score` >= min.score,]

# FILTER: peak size within range
max.peak.width <- 600
peaks <- peaks[peaks$End-peaks$Start <= max.peak.width,]

# FILTER: RUNX motif near peak AND no ETS motif that is closer to summit
max.motif.dist <- 50
peaks$closestRUNX <- suppressWarnings(sapply(strsplit(peaks$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`, ","), function(x) min(abs(as.numeric(x)))))
peaks$closestETS <- suppressWarnings(sapply(strsplit(peaks$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer Distance From Summit`, ","), function(x) min(abs(as.numeric(x)))))
peaks <- peaks[peaks$closestRUNX <= max.motif.dist,]
peaks <- peaks[peaks$closestRUNX < peaks$closestETS,]

