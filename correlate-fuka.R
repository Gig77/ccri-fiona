p <- 0.01

fuka.d20plus <- read.delim("/mnt/projects/chrisi/results/fuka/matAnn.telamlKD.REHandAT2.esetnsF.REH.AT2.balanced.annot.tsv")
fuka.d20plus <- fuka.d20plus[,c("syms", "Padj", "logFC")]
names(fuka.d20plus) <- c("Gene", "fuka.kdER.d20plus.padj", "fuka.kdER.d20plus.logfc")

fuka.d20plus.REH <- read.delim("/mnt/projects/chrisi/results/fuka/telamlKD.REHandAT2.esetnsF.onlyG_late_REH.annot.tsv")
fuka.d20plus.REH <- fuka.d20plus.REH[,c("syms", "Padj.onlyG_late_REH", "logFC.onlyG_late_REH")]
names(fuka.d20plus.REH) <- c("Gene", "fuka.kdER.d20plus.REH.padj", "fuka.kdER.d20plus.REH.logfc")

fuka.d20plus.AT2 <- read.delim("/mnt/projects/chrisi/results/fuka/telamlKD.REHandAT2.esetnsF.onlyG_late_AT2.annot.tsv")
fuka.d20plus.AT2 <- fuka.d20plus.AT2[,c("syms", "Padj.onlyG_late_AT2", "logFC.onlyG_late_AT2")]
names(fuka.d20plus.AT2) <- c("Gene", "fuka.kdER.d20plus.AT2.padj", "fuka.kdER.d20plus.AT2.logfc")

fiona.oeERvsEmptyB1 <- read.delim("/mnt/projects/fiona/results/anduril/execute/deseqAnnotated_oeERvsEmptyB1/table.csv", check.names = F, stringsAsFactors = F)
fiona.oeERvsEmptyB1 <- fiona.oeERvsEmptyB1[!is.na(fiona.oeERvsEmptyB1$Gene) & !is.na(fiona.oeERvsEmptyB1$q) & fiona.oeERvsEmptyB1$p <= p,]
fiona.oeERvsEmptyB1 <- fiona.oeERvsEmptyB1[order(fiona.oeERvsEmptyB1$q),]
fiona.oeERvsEmptyB1 <- fiona.oeERvsEmptyB1[!duplicated(fiona.oeERvsEmptyB1$Gene),]

fiona.oeERvsEmptyB2 <- read.delim("/mnt/projects/fiona/results/anduril/execute/deseqAnnotated_oeERvsEmptyB2/table.csv", check.names = F, stringsAsFactors = F)
fiona.oeERvsEmptyB2 <- fiona.oeERvsEmptyB2[!is.na(fiona.oeERvsEmptyB2$Gene) & !is.na(fiona.oeERvsEmptyB2$q) & fiona.oeERvsEmptyB2$p <= p,]
fiona.oeERvsEmptyB2 <- fiona.oeERvsEmptyB2[order(fiona.oeERvsEmptyB2$q),]
fiona.oeERvsEmptyB2 <- fiona.oeERvsEmptyB2[!duplicated(fiona.oeERvsEmptyB2$Gene),]

pdf("/mnt/projects/fiona/results/correlation-fuka.pdf", width=11, height=11)
par(mfrow=c(3,3), oma=c(1,1,1,1))

#-------------------
# batch 1 vs. batch2
#-------------------

plot.new()

fiona.oeERvsEmptyBoth <- merge(fiona.oeERvsEmptyB1[,c("ids", "fc")], fiona.oeERvsEmptyB2[,c("ids", "fc")], by="ids", suffixes=c(".B1", ".B2"))
notna <- !is.na(fiona.oeERvsEmptyBoth$fc.B1) & !is.na(fiona.oeERvsEmptyBoth$fc.B2)
x <- fiona.oeERvsEmptyBoth[notna, "fc.B1"]
y <- fiona.oeERvsEmptyBoth[notna, "fc.B2"]
plot(x, y, xlab="log2FC Fiona ERvsEmpty B1", ylab="log2FC Fiona ERvsEmpty B2", cex=0.4, pch=19, main="Fiona Batch 1 vs. Batch 2")
fit <- lm(y~x)
abline(fit, col="red", lwd=2)
test <- cor.test(x, y, method="spearman")
text((min(x)+max(x))/2, max(y), sprintf("R=%.2f/p=%.3g", test$estimate, test$p.value), adj=0.5)

plot.new()

#-------------------
# batch 1 vs. fuka
#-------------------

fiona.oeERvsEmptyB1 <- merge(fiona.oeERvsEmptyB1, fuka.d20plus, all.x=T)
fiona.oeERvsEmptyB1 <- merge(fiona.oeERvsEmptyB1, fuka.d20plus.REH, all.x=T)
fiona.oeERvsEmptyB1 <- merge(fiona.oeERvsEmptyB1, fuka.d20plus.AT2, all.x=T)

notna <- !is.na(fiona.oeERvsEmptyB1$fuka.kdER.d20plus.logfc) & !is.na(fiona.oeERvsEmptyB1$fc)
x <- fiona.oeERvsEmptyB1[notna, "fc"]
y <- fiona.oeERvsEmptyB1[notna, "fuka.kdER.d20plus.logfc"]
plot(x, y, xlab="log2FC Fiona ERvsEmpty B1", ylab="log2FC Fuka REH+AT2 D20+", cex=0.4, pch=19, main="Fiona Batch 1 vs. Fuka REH+AT2 D20+")
fit <- lm(y~x)
abline(fit, col="red", lwd=2)
test <- cor.test(x, y, method="spearman")
text((min(x)+max(x))/2, max(y), sprintf("R=%.2f/p=%.3g", test$estimate, test$p.value), adj=0.5)

# REH D20+
notna <- !is.na(fiona.oeERvsEmptyB1$fuka.kdER.d20plus.REH.logfc) & !is.na(fiona.oeERvsEmptyB1$fc)
x <- fiona.oeERvsEmptyB1[notna, "fc"]
y <- fiona.oeERvsEmptyB1[notna, "fuka.kdER.d20plus.REH.logfc"]
plot(x, y, xlab="log2FC Fiona ERvsEmpty B1", ylab="log2FC Fuka REH D20+", cex=0.4, pch=19, main="Fiona Batch 1 vs. Fuka REH D20+")
fit <- lm(y~x)
abline(fit, col="red", lwd=2)
test <- cor.test(x, y, method="spearman")
text((min(x)+max(x))/2, max(y), sprintf("R=%.2f/p=%.3g", test$estimate, test$p.value), adj=0.5)

# AT2 D20+
notna <- !is.na(fiona.oeERvsEmptyB1$fuka.kdER.d20plus.AT2.logfc) & !is.na(fiona.oeERvsEmptyB1$fc)
x <- fiona.oeERvsEmptyB1[notna, "fc"]
y <- fiona.oeERvsEmptyB1[notna, "fuka.kdER.d20plus.AT2.logfc"]
plot(x, y, xlab="log2FC Fiona ERvsEmpty B1", ylab="log2FC Fuka AT2 D20+", cex=0.4, pch=19, main="Fiona Batch 1 vs. Fuka AT2 D20+")
fit <- lm(y~x)
abline(fit, col="red", lwd=2)
test <- cor.test(x, y, method="spearman")
text((min(x)+max(x))/2, max(y), sprintf("R=%.2f/p=%.3g", test$estimate, test$p.value), adj=0.5)

#-------------------
# batch 2 vs. fuka
#-------------------

fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, fuka.d20plus, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, fuka.d20plus.REH, all.x=T)
fiona.oeERvsEmptyB2 <- merge(fiona.oeERvsEmptyB2, fuka.d20plus.AT2, all.x=T)

# REH+AT2 D20+
notna <- !is.na(fiona.oeERvsEmptyB2$fuka.kdER.d20plus.logfc) & !is.na(fiona.oeERvsEmptyB2$fc)
x <- fiona.oeERvsEmptyB2[notna, "fc"]
y <- fiona.oeERvsEmptyB2[notna, "fuka.kdER.d20plus.logfc"]
plot(x, y, xlab="log2FC Fiona ERvsEmpty B2", ylab="log2FC Fuka REH+AT2 D20+", cex=0.4, pch=19, main="Fiona Batch 2 vs. Fuka REH+AT2 D20+")
fit <- lm(y~x)
abline(fit, col="red", lwd=2)
test <- cor.test(x, y, method="spearman")
text((min(x)+max(x))/2, max(y), sprintf("R=%.2f/p=%.3g", test$estimate, test$p.value), adj=0.5)

# REH D20+
notna <- !is.na(fiona.oeERvsEmptyB2$fuka.kdER.d20plus.REH.logfc) & !is.na(fiona.oeERvsEmptyB2$fc)
x <- fiona.oeERvsEmptyB2[notna, "fc"]
y <- fiona.oeERvsEmptyB2[notna, "fuka.kdER.d20plus.REH.logfc"]
plot(x, y, xlab="log2FC Fiona ERvsEmpty B2", ylab="log2FC Fuka REH D20+", cex=0.4, pch=19, main="Fiona Batch 2 vs. Fuka REH D20+")
fit <- lm(y~x)
abline(fit, col="red", lwd=2)
test <- cor.test(x, y, method="spearman")
text((min(x)+max(x))/2, max(y), sprintf("R=%.2f/p=%.3g", test$estimate, test$p.value), adj=0.5)

# AT2 D20+
notna <- !is.na(fiona.oeERvsEmptyB2$fuka.kdER.d20plus.AT2.logfc) & !is.na(fiona.oeERvsEmptyB2$fc)
x <- fiona.oeERvsEmptyB2[notna, "fc"]
y <- fiona.oeERvsEmptyB2[notna, "fuka.kdER.d20plus.AT2.logfc"]
plot(x, y, xlab="log2FC Fiona ERvsEmpty B2", ylab="log2FC Fuka AT2 D20+", cex=0.4, pch=19, main="Fiona Batch 2 vs. Fuka AT2 D20+")
fit <- lm(y~x)
abline(fit, col="red", lwd=2)
test <- cor.test(x, y, method="spearman")
text((min(x)+max(x))/2, max(y), sprintf("R=%.2f/p=%.3g", test$estimate, test$p.value), adj=0.5)

dev.off()