source("/mnt/projects/fiona/scripts/read-peaks.R")

pdf("/mnt/projects/fiona/results/peak-width-distribution.pdf")

plot(density(reh$End-reh$Start, bw=40), col="blue", lwd=3, xlim=c(0,2000), xlab="Peak width (bp)", main="Peak width distribution")
lines(density(at2$End-at2$Start, bw=40), col="red", lwd=3)
lines(density(nalm6.er$End-nalm6.er$Start, bw=40), col="gray", lwd=3)
lines(density(nalm6.runx1$End-nalm6.runx1$Start, bw=40), col="black", lwd=3)
legend("topright", c("REH", "AT2", "NALM-6 E/R", "NALM-6 RUNX1"), fill=c("blue", "red", "gray", "black"))

dev.off()
