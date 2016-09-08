source("/mnt/projects/fiona/scripts/read-peaks.R")

#npeaks <- 1000
minpos <- -100
maxpos <- 100
missingpos <- 10000

runxCoords <- t(unlist(apply(at2, 1, function(x) {
    runx <- as.numeric(unlist(strsplit(x["RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit"], ",")))
    runx <- runx[runx >= minpos & runx <= maxpos]
    runx <- runx[order(abs(runx))][1:2]
    runx[is.na(runx)] <- mean(runx, na.rm = TRUE)
    ets <- as.numeric(unlist(strsplit(x["ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit"], ",")))
    ets <- ets[ets >= minpos & ets <= maxpos]
    ets <- ets[order(abs(ets))][1:2]
    ets[is.na(ets)] <- mean(ets, na.rm = TRUE)
    #pos <- pos[order(abs(pos))]
    #pos <- pos[order(pos)]
    pos <- c(runx, ets)
    pos[is.na(pos)] <- mean(pos, na.rm = TRUE)
    pos
})))
rownames(runxCoords) <- at2[,1]
runxCoords <- runxCoords[complete.cases(runxCoords),]
#runxCoords <- runxCoords[apply(runxCoords, 1, function(x) sum(x == missingpos) != 4),]
#peakOrder <- hclust(dist(runxCoords[,1:2]), method="average")$order
peakOrder <- order(runxCoords[,1], runxCoords[,2])
peakOrder <- order(rowMeans(runxCoords, na.rm=T))

plot(0, yaxt='n', bty='n', pch='', ylab='', xlab='Motif position', xlim=c(minpos,maxpos), ylim=c(0,nrow(runxCoords)))
line <- nrow(runxCoords)
for (peakid in rownames(runxCoords)[peakOrder]) {
  runx <- at2[which(at2[,1] == peakid), "RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit"]
  runx <- as.numeric(unlist(strsplit(runx, ",")))
  runx <- runx[runx >= minpos & runx <= maxpos]
  runx <- runx[order(abs(runx))][1:2]
  runx <- runx[!is.na(runx)]

  ets <- at2[which(at2[,1] == peakid), "ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit"]
  ets <- as.numeric(unlist(strsplit(ets, ",")))
  ets <- ets[ets >= minpos & ets <= maxpos]
  ets <- ets[order(abs(ets))][1:2]
  ets <- ets[!is.na(ets)]
  
  if (length(runx) > 0) {
    for (mpos in runx) {
      segments(mpos-1, line, mpos+1, line, col="red")
    }
  }
  if (length(ets) > 0) {
    for (mpos in ets) {
      segments(mpos-1, line, mpos+1, line, col="blue")
    }
  }
  
  line <- line - 1
}
abline(v=0, lty=2, col="gray")

