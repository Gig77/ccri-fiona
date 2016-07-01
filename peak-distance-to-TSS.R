source("/mnt/projects/fiona/scripts/read-peaks.R")

# E/R distance across cell lines

pdf ("/mnt/projects/fiona/results/peak-distance-to-TSS.pdf")

plot(NA, xlim=c(-50,50), ylim=c(0,0.05), xlab="Distance to nearest TSS (kb)", main="E/R peak distance to nearest TSS", ylab="Density", cex.lab=1.5, cex.axis=1.3)
polygon(density(c(at2.shuffled$`Distance to TSS`/1000, reh.shuffled$`Distance to TSS`/1000, nalm6.runx1.shuffled$`Distance to TSS`/1000, nalm6.er.shuffled$`Distance to TSS`/1000), n=5000, bw=3), col="lightgray", border="lightgray")
lines(density(at2$`Distance to TSS`/1000, n=5000, bw=3), xlim=c(-50,50), col="blue", lwd=4)
lines(density(reh$`Distance to TSS`/1000, n=5000, bw=3), xlim=c(-50,50), col="red", lwd=4)
lines(density(nalm6.er$`Distance to TSS`/1000, n=5000, bw=3), xlim=c(-50,50), col="black", lwd=4)
legend("topright", c("NALM6", "REH", "AT2", "random"), fill=c("black", "red", "blue", "lightgray"), cex=1.6)

# E/R unique vs. shared peaks comparing AT2 and REH

plot(NA, xlim=c(-50,50), ylim=c(0,0.05), xlab="Distance to nearest TSS (kb)", main="AT2 E/R peak distance to nearest TSS", ylab="Density", cex.lab=1.5, cex.axis=1.3)
polygon(density(at2.shuffled$`Distance to TSS`/1000, n=5000, bw=3), col="lightgray", border="lightgray")
lines(density(at2$`Distance to TSS`[unique(o.at2.reh@queryHits)]/1000, n=5000, bw=3), xlim=c(-50,50), col="blue", lwd=4)
lines(density(at2$`Distance to TSS`[-unique(o.at2.reh@queryHits)]/1000, n=5000, bw=3), xlim=c(-50,50), col="red", lwd=4)
legend("topright", c("shared with REH", "only AT2"), fill=c("blue", "red"), cex=1.6)

plot(NA, xlim=c(-50,50), ylim=c(0,0.05), xlab="Distance to nearest TSS (kb)", main="REH E/R peak distance to nearest TSS", ylab="Density", cex.lab=1.5, cex.axis=1.3)
polygon(density(reh.shuffled$`Distance to TSS`/1000, n=5000, bw=3), xlim=c(-50,50), col="lightgray", border="lightgray")
lines(density(reh$`Distance to TSS`[unique(o.at2.reh@subjectHits)]/1000, n=5000, bw=3), xlim=c(-50,50), col="blue", lwd=4)
lines(density(reh$`Distance to TSS`[-unique(o.at2.reh@subjectHits)]/1000, n=5000, bw=3), xlim=c(-50,50), col="red", lwd=4)
legend("topright", c("shared with AT2", "only REH"), fill=c("blue", "red"), cex=1.6)

# high vs. low scoring peaks

score.cutoff <- 30

plot(NA, xlim=c(-50,50), ylim=c(0,0.05), xlab="Distance to nearest TSS (kb)", main="AT E/R peak distance to nearest TSS", ylab="Density", cex.lab=1.5, cex.axis=1.3)
polygon(density(at2.shuffled$`Distance to TSS`/1000, n=5000, bw=3), col="lightgray", border="lightgray")
lines(density(at2$`Distance to TSS`[at2$`Peak Score` <= score.cutoff]/1000, n=5000, bw=3), col="blue", lwd=4)
lines(density(at2$`Distance to TSS`[at2$`Peak Score` > score.cutoff]/1000, n=5000, bw=3), col="red", lwd=4)
legend("topright", c(paste0("Peak score <= ", score.cutoff), paste0("Peak score > ", score.cutoff)), fill=c("blue", "red"), cex=1.6)

plot(NA, xlim=c(-50,50), ylim=c(0,0.05), xlab="Distance to nearest TSS (kb)", main="AT E/R peak distance to nearest TSS", ylab="Density", cex.lab=1.5, cex.axis=1.3)
polygon(density(at2.shuffled$`Distance to TSS`/1000, n=5000, bw=3), col="lightgray", border="lightgray")
lines(density(at2$`Distance to TSS`[at2$`Peak Score` <= score.cutoff]/1000, n=5000, bw=3), col="blue", lwd=4)
lines(density(at2$`Distance to TSS`[at2$`Peak Score` > score.cutoff]/1000, n=5000, bw=3), col="red", lwd=4)
legend("topright", c(paste0("Peak score <= ", score.cutoff), paste0("Peak score > ", score.cutoff)), fill=c("blue", "red"), cex=1.6)

plot(NA, xlim=c(-50,50), ylim=c(0,0.05), xlab="Distance to nearest TSS (kb)", main="NALM6 E/R peak distance to nearest TSS", ylab="Density", cex.lab=1.5, cex.axis=1.3)
polygon(density(nalm6.er.shuffled$`Distance to TSS`/1000, n=5000, bw=3), col="lightgray", border="lightgray")
lines(density(nalm6.er$`Distance to TSS`[nalm6.er$`Peak Score` <= score.cutoff]/1000, n=5000, bw=3), col="blue", lwd=4)
lines(density(nalm6.er$`Distance to TSS`[nalm6.er$`Peak Score` > score.cutoff]/1000, n=5000, bw=3), col="red", lwd=4)
legend("topright", c(paste0("Peak score <= ", score.cutoff), paste0("Peak score > ", score.cutoff)), fill=c("blue", "red"), cex=1.6)

plot(NA, xlim=c(-50,50), ylim=c(0,0.05), xlab="Distance to nearest TSS (kb)", main="NALM6 RUNX1 peak distance to nearest TSS", ylab="Density", cex.lab=1.5, cex.axis=1.3)
polygon(density(nalm6.runx1.shuffled$`Distance to TSS`/1000, n=5000, bw=3), col="lightgray", border="lightgray")
lines(density(nalm6.runx1$`Distance to TSS`[nalm6.runx1$`Peak Score` <= score.cutoff]/1000, n=5000, bw=3), col="blue", lwd=4)
lines(density(nalm6.runx1$`Distance to TSS`[nalm6.runx1$`Peak Score` > score.cutoff]/1000, n=5000, bw=3), col="red", lwd=4)
legend("topright", c(paste0("Peak score <= ", score.cutoff), paste0("Peak score > ", score.cutoff)), fill=c("blue", "red"), cex=1.6)

# NALM6 E/R vs. RUNX1

plot(NA, xlim=c(-50,50), ylim=c(0,0.05), xlab="Distance to nearest TSS (kb)", main="NALM6 peak distance to nearest TSS", ylab="Density", cex.lab=1.5, cex.axis=1.3)
#hist(nalm6.er$`Distance to TSS`/1000, xlim=c(-50,50), breaks=1000, freq = F, col="#0000FF50", add=T)
#hist(nalm6.runx1$`Distance to TSS`/1000, xlim=c(-50,50), breaks=1000, freq = F, col="#FF000050", add=T)
#hist(nalm6.shuffled$`Distance to TSS`/1000, xlim=c(-50,50), breaks=1000, freq = F, col="#00000030", add=T)
polygon(density(c(nalm6.er.shuffled$`Distance to TSS`/1000, nalm6.runx1.shuffled$`Distance to TSS`/1000), n=5000, bw=3), col="lightgray", border="lightgray")
lines(density(nalm6.er$`Distance to TSS`/1000, n=5000, bw=3), col="blue", lwd=4)
lines(density(nalm6.runx1$`Distance to TSS`/1000, n=5000, bw=3), col="red", lwd=4)
legend("topright", c("E/R", "RUNX1", "random"), fill=c("blue", "red", "lightgray"), cex=1.6)

# NALM6 E/R constitutive vs. de novo

plot(NA, xlim=c(-50,50), ylim=c(0,0.05), xlab="Distance to nearest TSS (kb)", main="NALM6 E/R peak distance to nearest TSS", ylab="Density", cex.lab=1.5, cex.axis=1.3)
polygon(density(nalm6.er.shuffled$`Distance to TSS`/1000, n=5000, bw=3), col="lightgray", border="lightgray")
lines(density(nalm6.er$`Distance to TSS`[unique(o.nalm6.runx1.er@subjectHits)]/1000, n=5000, bw=3), col="blue", lwd=4)
lines(density(nalm6.er$`Distance to TSS`[-unique(o.nalm6.runx1.er@subjectHits)]/1000, n=5000, bw=3), col="red", lwd=4)
legend("topright", c("constitutive", "de novo", "random"), fill=c("blue", "red", "lightgray"), cex=1.6)

dev.off()