options(warn=1)

# read annotation data sets

fuka.d20plus <- read.delim("/mnt/projects/chrisi/results/fuka/matAnn.telamlKD.REHandAT2.esetnsF.REH.AT2.balanced.annot.tsv")
fuka.d20plus <- fuka.d20plus[,c("syms", "Padj", "logFC")]
names(fuka.d20plus) <- c("Gene", "fuka.kdER.d20plus.padj", "fuka.kdER.d20plus.logfc")

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
veronika.E1 <- veronika.E1[order(veronika.E1$qValue),]
veronika.E1 <- veronika.E1[!duplicated(veronika.E1$Gene),]
veronika.E1 <- veronika.E1[,c("Gene", "qValue", "fc")]
names(veronika.E1) <- c("Gene", "veronika.E1.kdER.vs.empty.padj", "veronika.E1.kdER.vs.empty.logfc")

veronika.E2.d3 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseqAnnotated_ERvsNTd3/table.csv")
veronika.E2.d3 <- veronika.E2.d3[order(veronika.E2.d3$qValue),]
veronika.E2.d3 <- veronika.E2.d3[!duplicated(veronika.E2.d3$Gene),]
veronika.E2.d3 <- veronika.E2.d3[,c("Gene", "qValue", "fc")]
names(veronika.E2.d3) <- c("Gene", "veronika.E2.kdER.vs.empty.D3.padj", "veronika.E2.kdER.vs.empty.D3.logfc")

veronika.E2.d8 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseqAnnotated_ERvsNTd8/table.csv")
veronika.E2.d8 <- veronika.E2.d8[order(veronika.E2.d8$qValue),]
veronika.E2.d8 <- veronika.E2.d8[!duplicated(veronika.E2.d8$Gene),]
veronika.E2.d8 <- veronika.E2.d8[,c("Gene", "qValue", "fc")]
names(veronika.E2.d8) <- c("Gene", "veronika.E2.kdER.vs.empty.D8.padj", "veronika.E2.kdER.vs.empty.D8.logfc")

veronika.E2.d15 <- read.delim("/mnt/projects/veronika/results/anduril/execute/deseqAnnotated_ERvsNTd15/table.csv")
veronika.E2.d15 <- veronika.E2.d15[order(veronika.E2.d15$qValue),]
veronika.E2.d15 <- veronika.E2.d15[!duplicated(veronika.E2.d15$Gene),]
veronika.E2.d15 <- veronika.E2.d15[,c("Gene", "qValue", "fc")]
names(veronika.E2.d15) <- c("Gene", "veronika.E2.kdER.vs.empty.D15.padj", "veronika.E2.kdER.vs.empty.D15.logfc")

helena <- read.delim("/mnt/projects/helena_veronika/results/anduril/execute/deseqAnnotated_oeERvsEmpty/table.csv")
helena <- helena[order(helena$qValue),]
helena <- helena[!duplicated(helena$Gene),]
helena <- helena[,c("Gene", "qValue", "fc")]
names(helena) <- c("Gene", "helena.oeER.vs.empty.padj", "helena.oeER.vs.empty.logfc")

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

fiona.oeERvsEmpty <- read.delim("/mnt/projects/fiona/results/anduril/execute/deseqAnnotated_oeERvsEmpty/table.csv", check.names = F, stringsAsFactors = F)

fiona.oeERvsEmpty <- fiona.oeERvsEmpty[!is.na(fiona.oeERvsEmpty$Gene),]
fiona.oeERvsEmpty <- fiona.oeERvsEmpty[order(fiona.oeERvsEmpty$qValue),]
fiona.oeERvsEmpty <- fiona.oeERvsEmpty[!duplicated(fiona.oeERvsEmpty$Gene),]
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, fuka.d13, all.x=T)
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, fuka.d20plus, all.x=T)
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, boer.TA.vs.noTall, all.x=T)
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, boer.TA.vs.rest, all.x=T)
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, ross, all.x=T)
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, chrisi, all.x=T)
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, veronika.E1, all.x=T)
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, veronika.E2.d3, all.x=T)
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, veronika.E2.d8, all.x=T)
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, veronika.E2.d15, all.x=T)
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, helena, all.x=T)
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, fiona.chipseq.runx1.tssdist, all.x=T)
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, tijssen, all.x=T)
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, wilson, all.x=T)
fiona.oeERvsEmpty <- merge(fiona.oeERvsEmpty, niebuhr, all.x=T)

# compute logFC averages over all data sets
fcs <- data.frame(fiona.oeERvsEmpty$fc, 
                  -fiona.oeERvsEmpty$fuka.kdER.d13.logfc, 
                  -fiona.oeERvsEmpty$fuka.kdER.d20plus.logfc,
                  fiona.oeERvsEmpty$boer.TA.vs.noTall.logfc,
                  fiona.oeERvsEmpty$boer.TA.vs.rest.logfc,
                  fiona.oeERvsEmpty$ross.TA.vs.noTall.logfc,
                  fiona.oeERvsEmpty$chrisi.oeER.logfc,
                  -fiona.oeERvsEmpty$veronika.E1.kdER.vs.empty.logfc,
                  -fiona.oeERvsEmpty$veronika.E2.kdER.vs.empty.D3.logfc,
                  -fiona.oeERvsEmpty$veronika.E2.kdER.vs.empty.D8.logfc,
                  -fiona.oeERvsEmpty$veronika.E2.kdER.vs.empty.D15.logfc,
                  fiona.oeERvsEmpty$helena.oeER.vs.empty.logfc,
                  check.names = F, 
                  stringsAsFactors = F)
fcs[is.na(fcs)] <- 0
fiona.oeERvsEmpty$meanFC <- rowMeans(fcs)

write.table(fiona.oeERvsEmpty, "/mnt/projects/fiona/results/oeERvsEmpty.combined.tsv", sep = "\t", quote = F, col.names = T, row.names = F)
