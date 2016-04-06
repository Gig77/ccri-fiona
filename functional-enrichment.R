library(GO.db)
library(clusterProfiler)

# get background gene ids
library(biomaRt)
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") # GRCh37, v75
bm <- getBM(attributes=c("entrezgene"), mart=mart)
genes.bg <- as.character(unique(bm$entrezgene))

at2 <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
at2 <- at2[!grepl("^GL", at2$Chr),]
at2.entrez <- as.character(unique(at2$`Entrez ID`))

reh <- read.delim("/mnt/projects/fiona/results/homer/ChIP24_REH_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
reh <- reh[!grepl("^GL", reh$Chr),]
reh.entrez <- as.character(unique(reh$`Entrez ID`))

nalm6.er <- read.delim("/mnt/projects/fiona/results/homer/ChIP23_NALM6_ER_peaks.annotated.with-expr.tsv", stringsAsFactors = F, check.names = F)
nalm6.er <- nalm6.er[!grepl("^GL", nalm6.er$Chr),]
nalm6.er.entrez <- as.character(unique(nalm6.er$`Entrez ID`))

readHomerEnrichments <- function(dir) {
  df <- cbind(Category="GOTERM_BP_FAT", read.delim(paste0(dir, "biological_process.txt"), stringsAsFactor = FALSE, check.names = F))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "cellular_component.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_MF_FAT", read.delim(paste0(dir, "molecular_function.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="MSigDB", read.delim(paste0(dir, "msigdb.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="KEGG", read.delim(paste0(dir, "kegg.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "biocyc.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "chromosome.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "cosmic.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "gene3d.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "gwas.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "interactions.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "interpro.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "kegg.txt"), stringsAsFactor = FALSE, check.names = F)))
  #df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "lipidmaps.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "pathwayInteractionDB.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "pfam.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "prints.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "prosite.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "reactome.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "smart.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "smpdb.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- rbind(df, cbind(Category="GOTERM_CC_FAT", read.delim(paste0(dir, "wikipathways.txt"), stringsAsFactor = FALSE, check.names = F)))
  df <- df[order(df$logP),]
  names(df)[names(df)=="Genes in Term"] <- "Pop.Hits"
  names(df)[names(df)=="Gene Symbols"] <- "Genes"
  df
}

at2.gs <- readHomerEnrichments("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER/")
reh.gs <- readHomerEnrichments("/mnt/projects/fiona/results/homer/ChIP24_REH_ER/")
nalm6.runx1.gs <- readHomerEnrichments("/mnt/projects/fiona/results/homer/ChIP22_NALM6_RUNX1/")
nalm6.er.gs <- readHomerEnrichments("/mnt/projects/fiona/results/homer/ChIP23_NALM6_ER/")
at2.gs.wRunx1Motif <- readHomerEnrichments("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_runx1Motif/")
at2.gs.woRunx1Motif <- readHomerEnrichments("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_norunx1Motif/")
at2.gs.denovo <- readHomerEnrichments("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_denovo/")
at2.gs.constitutive <- readHomerEnrichments("/mnt/projects/fiona/results/homer/ChIP24_AT2_ER_constitutive/")
reh.gs.denovo <- readHomerEnrichments("/mnt/projects/fiona/results/homer/ChIP24_REH_ER_denovo/")
reh.gs.constitutive <- readHomerEnrichments("/mnt/projects/fiona/results/homer/ChIP24_REH_ER_constitutive/")
nalm6.er.gs.denovo <- readHomerEnrichments("/mnt/projects/fiona/results/homer/ChIP23_NALM6_ER_denovo/")
nalm6.er.gs.constitutive <- readHomerEnrichments("/mnt/projects/fiona/results/homer/ChIP23_NALM6_ER_constitutive/")

# scatter plot AT2 vs REH

at2.vs.reh <- merge(at2.gs[,c("Category", "TermID", "Term", "logP")], reh.gs[,c("Category", "TermID", "Term", "logP")], by=c("Category", "TermID", "Term"), suffixes = c(".at2", ".reh"), all=T)
sign <- at2.vs.reh$logP.at2 <= -1 | at2.vs.reh$logP.reh <= -1
pdf("/mnt/projects/fiona/results/enrichment/at2-vs-reh.pdf")
plot(-at2.vs.reh$logP.at2[sign], -at2.vs.reh$logP.reh[sign], cex=.3, xlab="-logP Enrichment AT2", ylab="-logP Enrichment REH")
dev.off()

# scatter plot NALM-6 RUNX1 vs. ER

nalm6.runx1.vs.er <- merge(nalm6.runx1.gs[,c("Category", "TermID", "Term", "logP", "Pop.Hits")], nalm6.er.gs[,c("Category", "TermID", "Term", "logP", "Pop.Hits")], by=c("Category", "TermID", "Term", "Pop.Hits"), suffixes = c(".runx1", ".er"), all=T)
sign <- nalm6.runx1.vs.er$Pop.Hits <= 5000 & (nalm6.runx1.vs.er$logP.runx1 <= -1 | nalm6.runx1.vs.er$logP.er <= -1)
pdf("/mnt/projects/fiona/results/enrichment/nalm6-runx1-vs-er.pdf")
plot(-nalm6.runx1.vs.er$logP.runx1[sign], -nalm6.runx1.vs.er$logP.er[sign], cex=.3, xlab="-logP Enrichment NALM-6 RUNX1", ylab="-logP Enrichment NALM-6 ER")
dev.off()

# scatter plot AT2 ER vs NALM-6 ER

at2.vs.nalm6.er <- merge(at2.gs[,c("Category", "TermID", "Term", "logP", "Pop.Hits")], nalm6.er.gs[,c("Category", "TermID", "Term", "logP", "Pop.Hits")], by=c("Category", "TermID", "Term", "Pop.Hits"), suffixes = c(".at2", ".nalm6.er"), all=T)
sign <- at2.vs.nalm6.er$Pop.Hits <= 5000 & (at2.vs.nalm6.er$logP.at2 <= -1 | at2.vs.nalm6.er$logP.nalm6.er <= -1)
pdf("/mnt/projects/fiona/results/enrichment/at2-vs-nalm6-er.pdf")
plot(-at2.vs.nalm6.er$logP.at2[sign], -at2.vs.nalm6.er$logP.nalm6.er[sign], cex=.3, xlab="-logP Enrichment AT2 ER", ylab="-logP Enrichment NALM-6 ER")
abline(0, 1)
#do.label <- order(abs(at2.vs.nalm6.er$logP.at2[sign]-at2.vs.nalm6.er$logP.nalm6.er[sign]), decreasing = T)[1:50]
#with(at2.vs.nalm6.er[sign,], text(-logP.at2[do.label]+2, -logP.nalm6.er[do.label]-2, Term[do.label], cex=0.3, adj=0))
dev.off()
at2.vs.nalm6.er[sign,][order(at2.vs.nalm6.er$logP.at2[sign]),][1:20,]
at2.vs.nalm6.er[sign,][order(at2.vs.nalm6.er$logP.nalm6.er[sign]),][1:20,]

# scatter plot AT2 w/ RUNX1 motif vs AT2 wo/ RUNX1 motif

at2.wRunx1Motif.vs.woRunx1Motif <- merge(at2.gs.wRunx1Motif[,c("Category", "TermID", "Term", "logP")], at2.gs.woRunx1Motif[,c("Category", "TermID", "Term", "logP")], by=c("Category", "TermID", "Term"), suffixes = c(".wRunx1Motif", ".woRunx1Motif"), all=T)
sign <- at2.wRunx1Motif.vs.woRunx1Motif$logP.wRunx1Motif <= -1 | at2.wRunx1Motif.vs.woRunx1Motif$logP.woRunx1Motif <= -1
pdf("/mnt/projects/fiona/results/enrichment/at2-withRunx1Motif-vs-withoutRunx1Motif.pdf")
plot(-at2.wRunx1Motif.vs.woRunx1Motif$logP.wRunx1Motif[sign], -at2.wRunx1Motif.vs.woRunx1Motif$logP.woRunx1Motif[sign], cex=.3, xlab="-logP Enrichment AT2 w/ RUNX1 motif", ylab="-logP Enrichment AT2 wo/ RUNX1 motif")
dev.off()
at2.wRunx1Motif.vs.woRunx1Motif[sign,][order(at2.wRunx1Motif.vs.woRunx1Motif$logP.wRunx1Motif[sign]),][1:10,]
at2.wRunx1Motif.vs.woRunx1Motif[sign,][order(at2.wRunx1Motif.vs.woRunx1Motif$logP.woRunx1Motif[sign]),][1:10,]

# scatter plot AT2 denovo vs. constitutive

at2.denovo.vs.constitutive <- merge(at2.gs.denovo[,c("Category", "TermID", "Term", "logP", "Fraction of Targets in Term")], at2.gs.constitutive[,c("Category", "TermID", "Term", "logP", "Fraction of Targets in Term")], by=c("Category", "TermID", "Term"), suffixes = c(".denovo", ".constitutive"), all=T)
sign <- (at2.denovo.vs.constitutive$logP.denovo <= -4 | at2.denovo.vs.constitutive$logP.constitutive <= -4) & (at2.denovo.vs.constitutive$`Fraction of Targets in Term.denovo` > 0.05 | at2.denovo.vs.constitutive$`Fraction of Targets in Term.constitutive` > 0.05)
pdf("/mnt/projects/fiona/results/enrichment/at2-denovo-vs-constitutive.pdf")
plot(at2.denovo.vs.constitutive$`Fraction of Targets in Term.denovo`[sign], 
     at2.denovo.vs.constitutive$`Fraction of Targets in Term.constitutive`[sign], cex=.3, 
     xlab="AT2 de novo", 
     ylab="AT2 constitutive",
     main="AT2 fraction of genes in significant gene set (logP <= -4)")
label <- at2.denovo.vs.constitutive$`Fraction of Targets in Term.denovo`[sign] > 0.15 & at2.denovo.vs.constitutive$`Fraction of Targets in Term.constitutive`[sign] < 0.15
topdiff <- at2.denovo.vs.constitutive[sign,][order(at2.denovo.vs.constitutive$`Fraction of Targets in Term.denovo`[sign]-at2.denovo.vs.constitutive$`Fraction of Targets in Term.constitutive`[sign], decreasing = TRUE),]
text(topdiff$`Fraction of Targets in Term.denovo`[1:10]+0.01, topdiff$`Fraction of Targets in Term.constitutive`[1:10], topdiff$Term[1:10], cex=.6, adj=0)
dev.off()

# scatter plot AT2 denovo vs. constitutive

reh.denovo.vs.constitutive <- merge(reh.gs.denovo[,c("Category", "TermID", "Term", "logP", "Fraction of Targets in Term")], reh.gs.constitutive[,c("Category", "TermID", "Term", "logP", "Fraction of Targets in Term")], by=c("Category", "TermID", "Term"), suffixes = c(".denovo", ".constitutive"), all=T)
sign <- (reh.denovo.vs.constitutive$logP.denovo <= -4 | reh.denovo.vs.constitutive$logP.constitutive <= -4) & (reh.denovo.vs.constitutive$`Fraction of Targets in Term.denovo` > 0.05 | reh.denovo.vs.constitutive$`Fraction of Targets in Term.constitutive` > 0.05)
pdf("/mnt/projects/fiona/results/enrichment/reh-denovo-vs-constitutive.pdf")
plot(reh.denovo.vs.constitutive$`Fraction of Targets in Term.denovo`[sign], 
     reh.denovo.vs.constitutive$`Fraction of Targets in Term.constitutive`[sign], cex=.3, 
     xlab="REH de novo", 
     ylab="REH constitutive",
     main="REH fraction of genes in significant gene set (logP <= -4)")
label <- reh.denovo.vs.constitutive$`Fraction of Targets in Term.denovo`[sign] > 0.15 & reh.denovo.vs.constitutive$`Fraction of Targets in Term.constitutive`[sign] < 0.15
topdiff <- reh.denovo.vs.constitutive[sign,][order(reh.denovo.vs.constitutive$`Fraction of Targets in Term.denovo`[sign]-reh.denovo.vs.constitutive$`Fraction of Targets in Term.constitutive`[sign], decreasing = TRUE),]
text(topdiff$`Fraction of Targets in Term.denovo`[1:15]+0.01, topdiff$`Fraction of Targets in Term.constitutive`[1:15], topdiff$Term[1:15], cex=.6, adj=0)
dev.off()

# scatter plot NALM-6 ER denovo vs. constitutive

nalm6.er.denovo.vs.constitutive <- merge(nalm6.er.gs.denovo[,c("Category", "TermID", "Term", "logP", "Fraction of Targets in Term")], nalm6.er.gs.constitutive[,c("Category", "TermID", "Term", "logP", "Fraction of Targets in Term")], by=c("Category", "TermID", "Term"), suffixes = c(".denovo", ".constitutive"), all=T)
sign <- (nalm6.er.denovo.vs.constitutive$logP.denovo <= -4 | nalm6.er.denovo.vs.constitutive$logP.constitutive <= -4) & (nalm6.er.denovo.vs.constitutive$`Fraction of Targets in Term.denovo` > 0.05 | nalm6.er.denovo.vs.constitutive$`Fraction of Targets in Term.constitutive` > 0.05)
pdf("/mnt/projects/fiona/results/enrichment/nalm6.er-denovo-vs-constitutive.pdf")
plot(nalm6.er.denovo.vs.constitutive$`Fraction of Targets in Term.denovo`[sign], 
     nalm6.er.denovo.vs.constitutive$`Fraction of Targets in Term.constitutive`[sign], cex=.3, 
     xlab="NALM6 ER de novo", 
     ylab="NALM6 ER constitutive",
     main="NALM6 ER fraction of genes in significant gene set (logP <= -4)")
label <- nalm6.er.denovo.vs.constitutive$`Fraction of Targets in Term.denovo`[sign] > 0.15 & nalm6.er.denovo.vs.constitutive$`Fraction of Targets in Term.constitutive`[sign] < 0.15
topdiff <- nalm6.er.denovo.vs.constitutive[sign,][order(nalm6.er.denovo.vs.constitutive$`Fraction of Targets in Term.denovo`[sign]-nalm6.er.denovo.vs.constitutive$`Fraction of Targets in Term.constitutive`[sign], decreasing = TRUE),]
text(topdiff$`Fraction of Targets in Term.denovo`[1:15]+0.01, topdiff$`Fraction of Targets in Term.constitutive`[1:15], topdiff$Term[1:15], cex=.6, adj=0)
dev.off()


# enrichment maps with top-100 enriched gene sets (script from max)

write.table(at2.gs[at2.gs$Category %in% c("GOTERM_BP_FAT", "GOTERM_CC_FAT", "GOTERM_MF_FAT"),], "/mnt/projects/fiona/results/enrichment/at2.enrichedGeneSets.homer.GO.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(at2.gs[at2.gs$Category %in% c("MSigDB") & !grepl("^MODULE", at2.gs$TermID),], "/mnt/projects/fiona/results/enrichment/at2.enrichedGeneSets.homer.MSigDB.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(at2.gs[at2.gs$Category %in% c("KEGG") & grepl("^hsa", at2.gs$TermID),], "/mnt/projects/fiona/results/enrichment/at2.enrichedGeneSets.homer.KEGG.tsv", col.names = T, row.names = F, sep = "\t", quote = F)

options(error=recover)
source("/mnt/projects/fiona/scripts/pGSEAfunctions.R")
source("/mnt/projects/fiona/scripts/clusterFctCatDAVID.v3.R")

clusterFctCat(
  Files=c("/mnt/projects/fiona/results/enrichment/at2.enrichedGeneSets.homer.GO.tsv"),
  whichcols=c(4),
  whichPcols=c(4),
  PvalornumRow="numRow",
  Pval=0.05,
  howmany=100,
  largeCat=1000,
  pdftif="pdf",
  cut.col.scheme=0.4,
  col.overRide="vhighsign",
  cex.labels=0.8)

clusterFctCat(
  Files=c("/mnt/projects/fiona/results/enrichment/at2.enrichedGeneSets.homer.MSigDB.tsv"),
  whichcols=c(4),
  whichPcols=c(4),
  PvalornumRow="numRow",
  Pval=0.05,
  howmany=100,
  largeCat=1000,
  pdftif="pdf",
  cut.col.scheme=0.4,
  col.overRide="vhighsign",
  cex.labels=0.8)

clusterFctCat(
  Files=c("/mnt/projects/fiona/results/enrichment/at2.enrichedGeneSets.homer.KEGG.tsv"),
  whichcols=c(4),
  whichPcols=c(4),
  PvalornumRow="numRow",
  Pval=0.05,
  howmany=100,
  largeCat=1000,
  pdftif="pdf",
  cut.col.scheme=0.4,
  col.overRide="highsign",
  cex.labels=0.8)


# GO enrichment

ego.at2 <- enrichGO(
  gene = at2.entrez,
  universe = unique(org.Hs.egGO2ALLEGS[["GO:0008150"]]),
  organism = "human",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)

pdf("/mnt/projects/fiona/results/enrichment-maps.pdf")
enrichMap(ego.at2, n=50, vertex.label.font = 1, vertex.label.cex = 0.7)
dev.off()

# compare enrichments

ck <- compareCluster(
  geneClusters = list(AT2=at2.entrez, REH=reh.entrez, "NALM6 ER"=nalm6.er.entrez), 
  fun = "enrichGO", 
  ont="BP", 
  universe=unique(org.Hs.egGO2ALLEGS[["GO:0008150"]]),
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)
plot(ck)

