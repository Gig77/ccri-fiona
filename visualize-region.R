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

chr <- "chr21" ; start <- 36114222 ; end <- 36433562  # RUNX1
chr <- "chr12" ; start <- 69195097 ; end <- 69208992  # MDM2
chr <- "chr2"  ; start <- 8808479 ; end <- 8832267    # ID2
chr <- "chr10" ; start <- 74017308 ; end <- 74050859  # DDIT4
chr <- "chr19" ; start <- 50915546 ; end <- 50940239  # SPIB

ideogram.track <- IdeogramTrack(genome = "hg19", chromosome = chr, showId=F, showTitle=T, 
                                name=chr, fontcolor.title="black", background.title="white", rotation.title=0, fontsize=15)
axis.track <- GenomeAxisTrack(labelPos="above", name=chr, lwd=0.5, cex=0.8, col="black", fontcolor="black", fontface=2)
gene.track <- BiomartGeneRegionTrack(genome = "hg19", chromosome = chr, start = start, end = end, 
                                     name = "", showId=TRUE, stacking="squish", collapseTranscripts = "longest",
                                     fontsize=16, fontcolor.group="black", background.title="white")


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

pdf("/mnt/projects/fiona/results/region-plots/test.pdf", height=5, width=12)
plotTracks(list(ideogram.track, axis.track, gene.track, at2.er, reh.er, nalm6.runx1, nalm6.er), 
           from=start, to=end, title.width=1,
           cex.title = 0.7, cex.axis=0.5, ylim = ylims)
dev.off()
