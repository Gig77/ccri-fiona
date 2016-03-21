options(warn=1)

# read annotation data sets

fuka.d20plus <- read.delim("/mnt/projects/chrisi/results/fuka/matAnn.telamlKD.REHandAT2.esetnsF.REH.AT2.balanced.annot.tsv")
fuka.d20plus <- fuka.d20plus[,c("syms", "Padj", "logFC")]
names(fuka.d20plus) <- c("Gene", "fuka.kdER.d20plus.padj", "fuka.kdER.d20plus.logfc")

fuka.d20plus.REH <- read.delim("/mnt/projects/chrisi/results/fuka/telamlKD.REHandAT2.esetnsF.onlyG_late_REH.annot.tsv")
fuka.d20plus.REH <- fuka.d20plus.REH[,c("syms", "Padj.onlyG_late_REH", "logFC.onlyG_late_REH")]
names(fuka.d20plus.REH) <- c("Gene", "fuka.kdER.d20plus.REH.padj", "fuka.kdER.d20plus.REH.logfc")

fuka.d20plus.AT2 <- read.delim("/mnt/projects/chrisi/results/fuka/telamlKD.REHandAT2.esetnsF.onlyG_late_AT2.annot.tsv")
fuka.d20plus.AT2 <- fuka.d20plus.AT2[,c("syms", "Padj.onlyG_late_AT2", "logFC.onlyG_late_AT2")]
names(fuka.d20plus.AT2) <- c("Gene", "fuka.kdER.d20plus.AT2.padj", "fuka.kdER.d20plus.AT2.logfc")

fuka.d13 <- read.delim("/mnt/projects/veronika/data/fuka/telamlKD.esetnsF.onlyG_13.annot.xls")
fuka.d13 <- fuka.d13[,c("syms", "Padj", "logFC")]
names(fuka.d13) <- c("Gene", "fuka.kdER.d13.padj", "fuka.kdER.d13.logfc")

boer.TA.vs.noTall <- read.delim("/mnt/projects/chrisi/data/RossBoer/NordischALL.esetnsF.annot.txt", check.names=F)
boer.TA.vs.noTall <- boer.TA.vs.noTall[,c("syms", "adjPval.TAvs.mean.noTall", "TAvs.mean.noTall")]
names(boer.TA.vs.noTall) <- c("Gene", "boer.TA.vs.noTall.padj", "boer.TA.vs.noTall.logfc")

boer.TA.vs.rest <- read.delim("/mnt/projects/chrisi/data/RossBoer/matAnn.GSE13351_BOER.eset_zfilt_th3_nsF.tsv", check.names=F)
boer.TA.vs.rest <- boer.TA.vs.rest[,c("syms", "adjP.TA_vs_rest", "TA_vs_rest")]
names(boer.TA.vs.rest) <- c("Gene", "boer.TA.vs.rest.padj", "boer.TA.vs.rest.logfc")

ross <- read.delim("/mnt/projects/chrisi/data/RossBoer/ROSS2.2003.esetnsF.annot.txt", check.names=F)
ross <- ross[order(ross$adjPval.TAvs.mean_noTALL),]
ross <- ross[!duplicated(ross$syms),]
ross <- ross[,c("syms", "adjPval.TAvs.mean_noTALL", "TAvs.mean_noTALL")]
names(ross) <- c("Gene", "ross.TA.vs.noTall.padj", "ross.TA.vs.noTall.logfc")

chrisi <- read.delim("/mnt/projects/chrisi/results/deseq/zuber+strobl-expressed-vs-notexpressed.deseq2.chipseq-annotated.tsv", check.names=F)
chrisi <- chrisi[order(chrisi$padj),]
chrisi <- chrisi[!duplicated(chrisi$hgnc_symbol),]
chrisi <- chrisi[,c("hgnc_symbol", "padj", "log2FoldChange")]
names(chrisi) <- c("Gene", "chrisi.oeER.padj", "chrisi.oeER.logfc")

veronika.E1 <- read.delim("/mnt/projects/helena_veronika/results/anduril/execute/deseqAnnotated_shG1vsNT/table.csv")
veronika.E1 <- veronika.E1[order(veronika.E1$q),]
veronika.E1 <- veronika.E1[!duplicated(veronika.E1$Gene),]
veronika.E1 <- veronika.E1[,c("Gene", "q", "fc")]
names(veronika.E1) <- c("Gene", "veronika.E1.kdER.vs.empty.padj", "veronika.E1.kdER.vs.empty.logfc")

veronika.E2.d3 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseqAnnotated_ERvsNTd3/table.csv")
veronika.E2.d3 <- veronika.E2.d3[order(veronika.E2.d3$q),]
veronika.E2.d3 <- veronika.E2.d3[!duplicated(veronika.E2.d3$Gene),]
veronika.E2.d3 <- veronika.E2.d3[,c("Gene", "q", "fc")]
names(veronika.E2.d3) <- c("Gene", "veronika.E2.kdER.vs.empty.D3.padj", "veronika.E2.kdER.vs.empty.D3.logfc")

veronika.E2.d8 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseqAnnotated_ERvsNTd8/table.csv")
veronika.E2.d8 <- veronika.E2.d8[order(veronika.E2.d8$q),]
veronika.E2.d8 <- veronika.E2.d8[!duplicated(veronika.E2.d8$Gene),]
veronika.E2.d8 <- veronika.E2.d8[,c("Gene", "q", "fc")]
names(veronika.E2.d8) <- c("Gene", "veronika.E2.kdER.vs.empty.D8.padj", "veronika.E2.kdER.vs.empty.D8.logfc")

veronika.E2.d15 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseqAnnotated_ERvsNTd15/table.csv")
veronika.E2.d15 <- veronika.E2.d15[order(veronika.E2.d15$q),]
veronika.E2.d15 <- veronika.E2.d15[!duplicated(veronika.E2.d15$Gene),]
veronika.E2.d15 <- veronika.E2.d15[,c("Gene", "q", "fc")]
names(veronika.E2.d15) <- c("Gene", "veronika.E2.kdER.vs.empty.D15.padj", "veronika.E2.kdER.vs.empty.D15.logfc")

helena <- read.delim("/mnt/projects/helena_veronika/results/anduril/execute/deseqAnnotated_oeERvsEmpty/table.csv")
helena <- helena[order(helena$q),]
helena <- helena[!duplicated(helena$Gene),]
helena <- helena[,c("Gene", "q", "fc")]
names(helena) <- c("Gene", "helena.oeER.vs.empty.padj", "helena.oeER.vs.empty.logfc")

fionaB1.oeERvsEmpty <- read.delim("/mnt/projects/fiona/results/anduril/execute/deseqAnnotated_oeERvsEmptyB1/table.csv")
fionaB1.oeERvsEmpty <- fionaB1.oeERvsEmpty[order(fionaB1.oeERvsEmpty$q),]
fionaB1.oeERvsEmpty <- fionaB1.oeERvsEmpty[!duplicated(fionaB1.oeERvsEmpty$Gene),]
fionaB1.oeERvsEmpty <- fionaB1.oeERvsEmpty[,c("Gene", "q", "fc")]
names(fionaB1.oeERvsEmpty) <- c("Gene", "fionaB1.oeER.vs.empty.padj", "fionaB1.oeER.vs.empty.logfc")

fiona.chipseq.runx1 <- read.delim("/mnt/projects/fiona/results/homer/runx1_peaks.annotated.with-expr.tsv", check.names = F, stringsAsFactors = F)
#fiona.chipseq.runx1 <- fiona.chipseq.runx1[fiona.chipseq.runx1$`Distance to TSS` > -2000 & fiona.chipseq.runx1$`Distance to TSS` < 1000,]
fiona.chipseq.runx1.tssdist <- aggregate(`Distance to TSS` ~ `Gene Name`, paste, collapse=",", data=fiona.chipseq.runx1)
names(fiona.chipseq.runx1.tssdist) <- c("Gene", "fiona.chipseq.runx1.TSS.distance")

# ChIP-seq Tijssen et al. 2011 (http://www.ncbi.nlm.nih.gov/pubmed/21571218)
tijssen <- read.csv("/mnt/projects/chrisi/results/chipseq/Tijssen_all.genes.txt",  stringsAsFactors=F, sep="\t", header=T, fill=T)
tijssen <- data.frame(Gene=tijssen$Runx1_alone, tijssen2011.chipseq.runx1=tijssen$Runx1_alone)

# ChIP-seq Wilson et al. 2010 (http://www.ncbi.nlm.nih.gov/pubmed/20887958)
wilson <- read.csv("/mnt/projects/chrisi/results/chipseq/Wilson_Gottgens_ChIPseq.txt",  stringsAsFactors=F, sep="\t", header=T, fill=T)
wilson <- wilson[!is.na(wilson$Runx1) & wilson$Runx1 != "",]
library(biomaRt)
human = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") # GRCh37, v75
mouse = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", dataset="mmusculus_gene_ensembl") # GRCm38, v75
humOrt <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values=wilson$Runx1, mart = mouse, attributesL = c("hgnc_symbol", "entrezgene"), martL = human)
wilson <- humOrt[!is.na(humOrt$EntrezGene.ID), c("HGNC.symbol", "MGI.symbol")]
wilson <- wilson[!duplicated(wilson),]
colnames(wilson) <- c("Gene", "wilson2010.chipseq.runx1.mouse")
wilson <- aggregate(wilson2010.chipseq.runx1.mouse~Gene, paste, collapse="|", data=wilson)

# ChIP-seq Niebuhr et. al 2013 (http://www.ncbi.nlm.nih.gov/pubmed/23704093)
niebuhr <-  read.csv("/mnt/projects/chrisi/results/chipseq/Niebuhr_TableS3_Runx1 Peaks Called in ProB-Cells.txt", stringsAsFactors=F, sep="\t", header=T, fill=T)
#niebuhr <- niebuhr[niebuhr$dist_tss > -5000 & niebuhr$dist_tss < 1000,]
#niebuhr <- niebuhr[niebuhr$score >= 100,]
humOrt <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values=niebuhr$nearest.gene, mart = mouse, attributesL = c("hgnc_symbol", "entrezgene"), martL = human)
niebuhr <- humOrt[!is.na(humOrt$EntrezGene.ID), c("HGNC.symbol", "MGI.symbol")]
niebuhr <- niebuhr[!duplicated(niebuhr),]
colnames(niebuhr) <- c("Gene", "niebuhr2013.chipseq.runx1.mouse")
niebuhr <- aggregate(niebuhr2013.chipseq.runx1.mouse~Gene, paste, collapse="|", data=niebuhr)

# merge

fiona.oeERvsEmptyB2 <- read.delim("/mnt/projects/fiona/results/anduril/execute/deseqAnnotated_oeERvsEmptyB2/table.csv", check.names = F, stringsAsFactors = F)

fiona.oeERvsEmptyB2 <- fiona.oeERvsEmptyB2[!is.na(fiona.oeERvsEmptyB2$Gene),]
fiona.oeERvsEmptyB2 <- fiona.oeERvsEmptyB2[order(fiona.oeERvsEmptyB2$q),]
fiona.oeERvsEmptyB2 <- fiona.oeERvsEmptyB2[!duplicated(fiona.oeERvsEmptyB2$Gene),]
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, fuka.d13, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, fuka.d20plus, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, fuka.d20plus.REH, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, fuka.d20plus.AT2, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, boer.TA.vs.noTall, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, boer.TA.vs.rest, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, ross, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, chrisi, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, veronika.E1, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, veronika.E2.d3, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, veronika.E2.d8, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, veronika.E2.d15, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, helena, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, fionaB1.oeERvsEmpty, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, fiona.chipseq.runx1.tssdist, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, tijssen, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, wilson, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, niebuhr, all.x=T)

# compute logFC averages over all data sets
fcs <- with(fiona.oeERvsEmptyB2, data.frame(
                  fc, 
                  -fuka.kdER.d13.logfc, 
                  -fuka.kdER.d20plus.logfc,
                  boer.TA.vs.noTall.logfc,
                  boer.TA.vs.rest.logfc,
                  ross.TA.vs.noTall.logfc,
                  chrisi.oeER.logfc,
                  -veronika.E1.kdER.vs.empty.logfc,
                  -veronika.E2.kdER.vs.empty.D3.logfc,
                  -veronika.E2.kdER.vs.empty.D8.logfc,
                  -veronika.E2.kdER.vs.empty.D15.logfc,
                  helena.oeER.vs.empty.logfc,
                  fionaB1.oeER.vs.empty.logfc,
                  check.names = F, 
                  stringsAsFactors = F))
fcs[is.na(fcs)] <- 0
fiona.oeERvsEmptyB2$meanFC <- rowMeans(fcs)

write.table(fiona.oeERvsEmptyB2, "/mnt/projects/fiona/results/oeERvsEmpty.combined.tsv", sep = "\t", quote = F, col.names = T, row.names = F, na="")
