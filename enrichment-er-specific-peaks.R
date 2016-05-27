library(GO.db)
library(clusterProfiler)
#source("https://bioconductor.org/biocLite.R") ; biocLite("RDAVIDWebService")

at2 <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
reh <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_REH_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.runx1 <- read.delim("/mnt/projects/fiona/results/homer/ChIP22_NALM6_RUNX1_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.er <- read.delim("/mnt/projects/fiona/results/homer/ChIP23_NALM6_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
er.denovo <- read.delim("/mnt/projects/fiona/results/er-consensus-peaks-not-runx1.xls", check.names = F)
er.better <- read.delim("/mnt/projects/fiona/results/er-consensus-peaks-better-than-runx1.xls", check.names = F)
er.diffbind <- read.delim("/mnt/projects/fiona/results/homer/ER.diffbind_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)

# MSigDB 5.0 gene sets
gmt <- read.delim("/mnt/projects/generic/data/msigdb5.0/msigdb.v5.0.symbols.gmt", stringsAsFactors = F, header = FALSE)
term2gene <- do.call("rbind", apply(gmt[gmt$V3 != "",], 1, function(x) {
  data.frame(term=x[1], gene=unique(x[3:length(x)][x[3:length(x)]!=""]), row.names = NULL)
}))

enrichAll <- function(gene.entrez=NULL, gene.hgnc=NULL, universe.entrez=NULL, universe.hgnc=NULL, pvalueCutoff = 0.05, qvalueCutoff = 0.2) {
  
  result <- data.frame(Category=character(0), ID=character(0), Description=character(0), GeneRatio=character(0), BgRatio=character(0), pvalue=numeric(0), p.adjust=numeric(0), qvalue=numeric(0), Count=integer(0))
  
  universe.mf <- unique(org.Hs.egGO2ALLEGS[["GO:0003674"]])
  ego.mf <- enrichGO(
    gene = gene.entrez,
    universe = if(is.null(universe.entrez)) universe.mf else unique(intersect(universe.entrez, universe.mf)),
    organism = "human",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    readable = TRUE
  )
  if (nrow(ego.mf@result) > 0) result <- rbind(result, cbind(data.frame(Category="GO molecular function"), ego.mf@result))

  universe.bp <- unique(org.Hs.egGO2ALLEGS[["GO:0008150"]])
  ego.bp <- enrichGO(
    gene = gene.entrez,
    universe = if(is.null(universe.entrez)) universe.bp else unique(intersect(universe.entrez, universe.bp)),
    organism = "human",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    readable = TRUE
  )
  if (nrow(ego.bp@result) > 0) result <- rbind(result, cbind(data.frame(Category="GO biological process"), ego.bp@result))
  
  universe.cc <- unique(org.Hs.egGO2ALLEGS[["GO:0005575"]])
  ego.cc <- enrichGO(
    gene = gene.entrez,
    universe = if(is.null(universe.entrez)) universe.cc else unique(intersect(universe.entrez, universe.cc)),
    organism = "human",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    readable = TRUE
  )
  if (nrow(ego.cc@result) > 0) result <- rbind(result, cbind(data.frame(Category="GO cellular component"), ego.cc@result))
  
  kegg <- enrichKEGG(
    gene = gene.entrez,
    universe = universe.entrez,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    readable = TRUE
  )
  if (nrow(kegg@result) > 0) result <- rbind(result, cbind(data.frame(Category="KEGG"), kegg@result))
  
  msigdb <- enricher(  
    gene          = gene.hgnc, 
    pAdjustMethod = "BH", 
    universe      = universe.hgnc, 
    minGSSize     = 5, 
    pvalueCutoff  = pvalueCutoff,
    qvalueCutoff  = qvalueCutoff,
    TERM2GENE     = term2gene, 
    TERM2NAME     = gmt[,c(1,2)]
  )
  if (nrow(msigdb@result) > 0) result <- rbind(result, cbind(data.frame(Category="MSigDB5.0"), msigdb@result))
  
  result[order(result$pvalue),]
}

# enrichment E/R de novo peaks
enr.denovo <- enrichAll(
  gene.entrez=as.character(unique(er.denovo$`Entrez ID`)),
  gene.hgnc=unique(er.denovo$`Gene Name`),
  universe.entrez=as.character(unique(c(at2$`Entrez ID`, reh$`Entrez ID`, nalm6.runx1$`Entrez ID`, nalm6.er$`Entrez ID`))),
  universe.hgnc=unique(c(at2$`Gene Name`, reh$`Gene Name`, nalm6.runx1$`Gene Name`, nalm6.er$`Gene Name`)),
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# enrichment E/R peaks that are better than RUNX1 peaks
enr.better <- enrichAll(
  gene.entrez=as.character(unique(er.better$`Entrez ID`)),
  gene.hgnc=as.character(unique(er.better$`Gene Name`)),
  universe.entrez=as.character(unique(c(at2$`Entrez ID`, reh$`Entrez ID`, nalm6.runx1$`Entrez ID`, nalm6.er$`Entrez ID`))),
  universe.hgnc=as.character(unique(c(at2$`Gene Name`, reh$`Gene Name`, nalm6.runx1$`Gene Name`, nalm6.er$`Gene Name`))),
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

write.table(enr.better, "/mnt/projects/fiona/results/genesets-enriched-in-er-peaks-better-than-runx1-peaks.tsv", col.names = T, row.names = F, quote = F, sep = "\t")

# enrichment differentially bound E/R peaks (DiffBind)
enr.diffbind <- enrichAll(
  gene.entrez=as.character(unique(er.diffbind$`Entrez ID`)),
  gene.hgnc=unique(er.diffbind$`Gene Name`),
  universe.entrez=as.character(unique(c(at2$`Entrez ID`, reh$`Entrez ID`, nalm6.runx1$`Entrez ID`, nalm6.er$`Entrez ID`))),
  universe.hgnc=unique(c(at2$`Gene Name`, reh$`Gene Name`, nalm6.runx1$`Gene Name`, nalm6.er$`Gene Name`)),
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

write.table(enr.diffbind, "/mnt/projects/fiona/results/genesets-enriched-in-differentially-bound-peaks-ER-vs-runx1.diffbind.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
