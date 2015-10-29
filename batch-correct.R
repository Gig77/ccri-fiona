library(DESeq2)
library(edgeR)
library(sva)

samples <- read.delim("/mnt/projects/fiona/results/anduril/execute/_qcReport_counts_array1/array/_index")
rownames(samples) <- samples$Key

samples$type <- NA
samples$type[grepl("ER", rownames(samples))] <- "ER"
samples$type[grepl("Em", rownames(samples))] <- "Empty"
samples$type[grepl("RHD", rownames(samples))] <- "RHD"
samples$type <- as.factor(samples$type)

samples$batch <- NA
samples$batch[grepl("_1", rownames(samples))] <- "harvest_early"
samples$batch[grepl("_2", rownames(samples))] <- "harvest_late"
samples$batch <- as.factor(samples$batch)

cds <- DESeqDataSetFromHTSeqCount(sampleTable=samples, directory="/", design=~1)
expressed <- rowSums(counts(cds)) >= 10
dge <- DGEList(counts=counts(cds[expressed,]))
dge.norm <- calcNormFactors(dge, method="TMM")
y <- voom(dge.norm)

# correct batch
model.full <- model.matrix(~type, data=samples) 
y.corrected <- ComBat(dat=y$E, batch=samples$batch, mod=model.full)

# plot PCA before and after batch effect correction
pdf("/mnt/projects/fiona/results/batch-effect-correction.pdf", height=10, width=6)
par(mar=c(4,4,4,7), mfrow=c(2, 1))
par(xpd=TRUE)

pca <- prcomp(t(y$E), center = TRUE, scale. = TRUE)
xspan <- max(pca$x[,1])-min(pca$x[,1])
plot(pca$x[,1], pca$x[,2], col="white", xlab="PC1", ylab="PC2", xlim=c(min(pca$x[,1])-xspan/10, max(pca$x[,1])+xspan/10), main="Before batch correction")
text(pca$x[,1], pca$x[,2], samples$Key, col=as.numeric(samples$type), cex=0.8)
legend(x=par("usr")[2], y=max(par("usr")[3:4]), levels(samples$type), col=1:length(levels(samples$type)), pch=19)

pca <- prcomp(t(y.corrected), center = TRUE, scale. = TRUE)
xspan <- max(pca$x[,1])-min(pca$x[,1])
plot(pca$x[,1], pca$x[,2], col="white", xlab="PC1", ylab="PC2", xlim=c(min(pca$x[,1])-xspan/10, max(pca$x[,1])+xspan/10), main="After batch correction")
text(pca$x[,1], pca$x[,2], samples$Key, col=as.numeric(samples$type), cex=0.8)
legend(x=par("usr")[2], y=max(par("usr")[3:4]), levels(samples$type), col=1:length(levels(samples$type)), pch=19)
dev.off()