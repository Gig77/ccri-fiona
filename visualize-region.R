library(Gviz)
options(ucscChromosomeNames=FALSE)

# setup coverage tracks, will be subsetted to genomic region of interest later

at2.er.chip <-       DataTrack(range = "/mnt/projects/fiona/results/homer-tagfiles/ChIP24_AT2_ER_chip.ucsc.bedGraph.gz",       genome = "hg19", type = c("polygon", "g"), h = 3, v = 0, col.grid = "#DDDDDD", name = "AT2 E/R")
reh.er.chip <-       DataTrack(range = "/mnt/projects/fiona/results/homer-tagfiles/ChIP24_REH_ER_chip.ucsc.bedGraph.gz",       genome = "hg19", type = c("polygon", "g"), h = 3, v = 0, col.grid = "#DDDDDD", name = "REH E/R")
nalm6.runx1.chip <-  DataTrack(range = "/mnt/projects/fiona/results/homer-tagfiles/ChIP22_NALM6_RUNX1_chip.ucsc.bedGraph.gz",  genome = "hg19", type = c("polygon", "g"), h = 3, v = 0, col.grid = "#DDDDDD", name = "NALM6 RUNX1")
nalm6.er.chip <-     DataTrack(range = "/mnt/projects/fiona/results/homer-tagfiles/ChIP23_NALM6_ER_chip.ucsc.bedGraph.gz",     genome = "hg19", type = c("polygon", "g"), h = 3, v = 0, col.grid = "#DDDDDD", name = "NALM6 E/R")
at2.er.input <-      DataTrack(range = "/mnt/projects/fiona/results/homer-tagfiles/ChIP24_AT2_ER_input.ucsc.bedGraph.gz",      genome = "hg19", type = "polygon", name = "AT2 E/R", col.mountain = "darkgreen", fill.mountain = c("#FFCCFF", "#CCFFFF"))
reh.er.input <-      DataTrack(range = "/mnt/projects/fiona/results/homer-tagfiles/ChIP24_REH_ER_input.ucsc.bedGraph.gz",      genome = "hg19", type = "polygon", name = "REH E/R", col.mountain = "darkgreen", fill.mountain = c("#FFCCFF", "#CCFFFF"))
nalm6.runx1.input <- DataTrack(range = "/mnt/projects/fiona/results/homer-tagfiles/ChIP22_NALM6_RUNX1_input.ucsc.bedGraph.gz", genome = "hg19", type = "polygon", name = "NALM6 RUNX1", col.mountain = "darkgreen", fill.mountain = c("#FFCCFF", "#CCFFFF"))
nalm6.er.input <-    DataTrack(range = "/mnt/projects/fiona/results/homer-tagfiles/ChIP23_NALM6_ER_input.ucsc.bedGraph.gz",    genome = "hg19", type = "polygon", name = "NALM6 E/R", col.mountain = "darkgreen", fill.mountain = c("#FFCCFF", "#CCFFFF"))

regions = data.frame(name=character(0), chr=character(0), start=numeric(0), end=numeric(0), collapseTranscripts=character(0), stringsAsFactors = F)
regions = rbind(regions, setNames(data.frame("RUNX1", "chr21", 36137634, 36427347, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("RUNX1 Promoter", "chr21", 36408083, 36431709, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("MDM2", "chr12", 69195097, 69208992, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("ID2", "chr2", 8808479, 8832267, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("DDIT4", "chr10", 74017308, 74050859, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("SPIB", "chr19", 50915546, 50940239, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("BCL6", "chr3", 187438180, 187470003, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("EPOR", "chr19", 11479311, 11502877, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("MAD2L1", "chr4", 120972040, 120994919, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("VEGFA", "chr6", 43727798, 43770804, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("CDKN1A", "chr6", 36639099, 36659153, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("GZMB", "chr14", 25095203, 25108164, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("DNMT3a", "chr2", 25438556, 25592146, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("ETV6", "chr12", 11750000, 12100000, "meta", stringsAsFactors = F), names(regions)))

# E/R "unique" peaks

regions = rbind(regions, setNames(data.frame("LOC100129046", "chr1", 94050000, 94072512, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("SEMA4A", "chr1", 156104807, 156148343, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("TSPEAR", "chr21", 45900577, 45949376, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("SH3BP1", "chr22", 38026127, 38041161, "meta", stringsAsFactors = F), names(regions)))
regions = rbind(regions, setNames(data.frame("ZCCHC5", "chrX", 77907787, 77921669, "meta", stringsAsFactors = F), names(regions)))

# PCR primer coordinates

primer <- data.frame(name=character(0), seqnames=character(0), start=numeric(0), end=numeric(0))
primer <- rbind(primer, setNames(data.frame("RUNX1 P1 fwd", "chr21", 36421571, 36421592, stringsAsFactors = F), names(primer)))
primer <- rbind(primer, setNames(data.frame("RUNX1 P1 rev", "chr21", 36421490, 36421513, stringsAsFactors = F), names(primer)))
primer <- rbind(primer, setNames(data.frame("RUNX1 neg fwd", "chr21", 36410740, 36410759, stringsAsFactors = F), names(primer)))
primer <- rbind(primer, setNames(data.frame("RUNX1 neg rev", "chr21", 36410622, 36410641, stringsAsFactors = F), names(primer)))
primer <- rbind(primer, setNames(data.frame("MDM2 RUNX motif fwd", "chr12", 69202705, 69202726, stringsAsFactors = F), names(primer)))
primer <- rbind(primer, setNames(data.frame("MDM2 RUNX motif rev", "chr12", 69202872, 69202895, stringsAsFactors = F), names(primer)))
primer <- rbind(primer, setNames(data.frame("MDM2 RUNX neg fwd", "chr12", 69207256, 69207277, stringsAsFactors = F), names(primer)))
primer <- rbind(primer, setNames(data.frame("MDM2 RUNX neg rev", "chr12", 69207357, 69207382, stringsAsFactors = F), names(primer)))
primer <- rbind(primer, setNames(data.frame("ID-2 fwd #3", "chr2", 8822297, 8822316, stringsAsFactors = F), names(primer)))
primer <- rbind(primer, setNames(data.frame("ID-2 rev #3", "chr2", 8822358, 8822377, stringsAsFactors = F), names(primer)))
primer.track <- AnnotationTrack(makeGRangesFromDataFrame(primer), name="Primer", stacking="dense", min.height=1, rot.title=0, 
                            col="black", lwd.border=1, fontsize=10, col.title="white")

pdf("/mnt/projects/fiona/results/region-plots/regions.chipseq.coverage.pdf", height=8, width=12)

for (i in 1:nrow(regions)) {

  name <- regions$name[i]
  chr <- regions$chr[i]
  start <- regions$start[i]
  end <- regions$end[i]
  collapseTranscripts <- regions$collapseTranscripts[i]
  
  ideogram.track <- IdeogramTrack(genome = "hg19", chromosome = chr, showId=F, showTitle=T, 
                                  name=chr, fontcolor.title="black", background.title="white", rotation.title=0, fontsize=15)
  axis.track <- GenomeAxisTrack(labelPos="above", name=chr, lwd=0.5, cex=0.8, col="black", fontcolor="black", fontface=2)
  
  # weird behaviour of gviz: it won't give us any genes with ucsc chr identifiers, but if we use no 'chr' prefix we have to add it afterwards again
  gene.track <- BiomartGeneRegionTrack(genome = "hg19", chromosome = gsub("chr", "", chr), start = start, end = end, 
                                       name = "", showId=TRUE, stacking="squish", collapseTranscripts = collapseTranscripts,
                                       fontsize=16, fontcolor.group="black", background.title="white")
  seqlevels(ranges(gene.track)) <- sprintf("chr%s", seqlevels(ranges(gene.track)))
  chromosome(gene.track)=chr
  
  chromosome(primer.track) <- chr ; primer.track.region <- subset(primer.track, start, end)
  
  chromosome(at2.er.chip) <- chr ; at2.er.chip.region <- subset(at2.er.chip, start, end)
  chromosome(at2.er.input) <- chr ; at2.er.input.region <- subset(at2.er.input, start, end)
  at2.er <- OverlayTrack(list(at2.er.chip.region, at2.er.input.region))
  
  chromosome(reh.er.chip) <- chr ; reh.er.chip.region <- subset(reh.er.chip, start, end)
  chromosome(reh.er.input) <- chr ; reh.er.input.region <- subset(reh.er.input, start, end)
  reh.er <- OverlayTrack(list(reh.er.chip.region, reh.er.input.region))
  
  chromosome(nalm6.runx1.chip) <- chr ; nalm6.runx1.chip.region <- subset(nalm6.runx1.chip, start, end)
  chromosome(nalm6.runx1.input) <- chr ; nalm6.runx1.input.region <- subset(nalm6.runx1.input, start, end)
  nalm6.runx1 <- OverlayTrack(list(nalm6.runx1.chip.region, nalm6.runx1.input.region))
  
  chromosome(nalm6.er.chip) <- chr ; nalm6.er.chip.region <- subset(nalm6.er.chip, start, end)
  chromosome(nalm6.er.input) <- chr ; nalm6.er.input.region <- subset(nalm6.er.input, start, end)
  nalm6.er <- OverlayTrack(list(nalm6.er.chip.region, nalm6.er.input.region))
  
  ylims <- extendrange(range(c(values(at2.er.chip.region), values(at2.er.input.region),
                               values(reh.er.chip.region), values(reh.er.input.region),
                               values(nalm6.er.chip.region), values(nalm6.er.input.region), 
                               values(nalm6.runx1.chip.region), values(nalm6.runx1.input.region))))
  
  if (length(primer.track.region) > 0) {
    plotTracks(list(ideogram.track, axis.track, gene.track, primer.track.region, at2.er, reh.er, nalm6.er, nalm6.runx1), 
               from=start, to=end, title.width=1, main=name,
               cex.title = 0.9, cex.axis=0.5, fontcolor.title="black", col.axis="black", ylim = ylims)
  } else {
    plotTracks(list(ideogram.track, axis.track, gene.track, at2.er, reh.er, nalm6.er, nalm6.runx1), 
               from=start, to=end, title.width=1, main=name,
               cex.title = 0.9, cex.axis=0.5, fontcolor.title="black", col.axis="black", ylim = ylims)
  }
}
dev.off()
