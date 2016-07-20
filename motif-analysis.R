library(reshape2)

# read data
source("/mnt/projects/fiona/scripts/read-peaks.R")

pdf("/mnt/projects/fiona/results/motif-analysis.pdf", paper="a4")

# overview barchart

d <- data.frame("CL" = character(0), "Motif" = character(0), "Pct" = numeric(0), stringsAsFactors = F)
d <- rbind(d, setNames(data.frame("AT2 E/R", "RUNX1", sum(!is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / nrow(at2)), names(d)))
d <- rbind(d, setNames(data.frame("AT2 E/R", "ETS", sum(!is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / nrow(at2)), names(d)))
d <- rbind(d, setNames(data.frame("AT2 E/R", "EBF", sum(!is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / nrow(at2)), names(d)))
d <- rbind(d, setNames(data.frame("REH E/R", "RUNX1", sum(!is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / nrow(reh)), names(d)))
d <- rbind(d, setNames(data.frame("REH E/R", "ETS", sum(!is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / nrow(reh)), names(d)))
d <- rbind(d, setNames(data.frame("REH E/R", "EBF", sum(!is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / nrow(reh)), names(d)))
d <- rbind(d, setNames(data.frame("NALM6 RUNX1", "RUNX1", sum(!is.na(nalm6.runx1$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / nrow(nalm6.runx1)), names(d)))
d <- rbind(d, setNames(data.frame("NALM6 RUNX1", "ETS", sum(!is.na(nalm6.runx1$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / nrow(nalm6.runx1)), names(d)))
d <- rbind(d, setNames(data.frame("NALM6 RUNX1", "EBF", sum(!is.na(nalm6.runx1$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / nrow(nalm6.runx1)), names(d)))
d <- rbind(d, setNames(data.frame("NALM6 E/R", "RUNX1", sum(!is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / nrow(nalm6.er)), names(d)))
d <- rbind(d, setNames(data.frame("NALM6 E/R", "ETS", sum(!is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / nrow(nalm6.er)), names(d)))
d <- rbind(d, setNames(data.frame("NALM6 E/R", "EBF", sum(!is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / nrow(nalm6.er)), names(d)))
d$CL <- factor(d$CL, levels=c("AT2 E/R", "REH E/R", "NALM6 E/R", "NALM6 RUNX1"))
d$Motif <- factor(d$Motif, levels=c("RUNX1", "ETS", "EBF"))

d.shuffled <- data.frame("CL" = character(0), "Motif" = character(0), "Pct" = numeric(0), stringsAsFactors = F)
d.shuffled <- rbind(d.shuffled, setNames(data.frame("AT2 E/R", "RUNX1", sum(!is.na(at2.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / nrow(at2.shuffled)), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("AT2 E/R", "ETS", sum(!is.na(at2.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / nrow(at2.shuffled)), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("AT2 E/R", "EBF", sum(!is.na(at2.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / nrow(at2.shuffled)), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("REH E/R", "RUNX1", sum(!is.na(reh.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / nrow(reh.shuffled)), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("REH E/R", "ETS", sum(!is.na(reh.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / nrow(reh.shuffled)), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("REH E/R", "EBF", sum(!is.na(reh.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / nrow(reh.shuffled)), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("NALM6 RUNX1", "RUNX1", sum(!is.na(nalm6.runx1.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / nrow(nalm6.runx1.shuffled)), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("NALM6 RUNX1", "ETS", sum(!is.na(nalm6.runx1.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / nrow(nalm6.runx1.shuffled)), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("NALM6 RUNX1", "EBF", sum(!is.na(nalm6.runx1.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / nrow(nalm6.runx1.shuffled)), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("NALM6 E/R", "RUNX1", sum(!is.na(nalm6.er.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / nrow(nalm6.er.shuffled)), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("NALM6 E/R", "ETS", sum(!is.na(nalm6.er.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / nrow(nalm6.er.shuffled)), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("NALM6 E/R", "EBF", sum(!is.na(nalm6.er.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / nrow(nalm6.er.shuffled)), names(d.shuffled)))
d.shuffled$CL <- factor(d.shuffled$CL, levels=c("AT2 E/R", "REH E/R", "NALM6 E/R", "NALM6 RUNX1"))
d.shuffled$Motif <- factor(d.shuffled$Motif, levels=c("RUNX1", "ETS", "EBF"))

#pdf("/mnt/projects/fiona/results/temp.pdf", paper="a4")
ggplot(data = d, aes(x = CL, y = Pct, fill = Motif)) +   
  geom_bar(position = position_dodge(), stat = "identity") +
#  geom_bar(data = d.shuffled, aes(group=Motif), fill="black", alpha=0.2, position=position_dodge(), stat = "identity") +
  geom_text(data = d.shuffled, label="---", color = "white", alpha = 1, size=5, position=position_dodge(width=0.9), stat = "identity") +
  labs(fill="Motif", x="Cell line", y="% of peaks with motif") + 
  ggtitle("Percentage of peaks with known motif") +
  theme_bw() +
  scale_fill_manual(values = c("red", "darkgray", "blue")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.9),
        axis.title.y = element_text(size=15, vjust=1.9))
#dev.off()

# number of motifs per peak

d <- data.frame("CL" = character(0), "Motif" = character(0), "CountMotifs" = character(0), "CountPeaks" = numeric(0), stringsAsFactors = F)
d <- rbind(d, setNames(data.frame("AT2", "RUNX1", melt(prop.table(table(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "RUNX1", melt(prop.table(table(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6 ER", "RUNX1", melt(prop.table(table(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6 RUNX1", "RUNX1", melt(prop.table(table(nalm6.runx1$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "ETS", melt(prop.table(table(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "ETS", melt(prop.table(table(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6 ER", "ETS", melt(prop.table(table(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6 RUNX1", "ETS", melt(prop.table(table(nalm6.runx1$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "EBF", melt(prop.table(table(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "EBF", melt(prop.table(table(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6 ER", "EBF", melt(prop.table(table(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6 RUNX1", "EBF", melt(prop.table(table(nalm6.runx1$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)))), names(d)))
#d <- rbind(d, setNames(data.frame("AT2", "ETS:RUNX", melt(prop.table(table(at2$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d)))
#d <- rbind(d, setNames(data.frame("REH", "ETS:RUNX", melt(prop.table(table(reh$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d)))
#d <- rbind(d, setNames(data.frame("NALM-6 ER", "ETS:RUNX", melt(prop.table(table(nalm6.er$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d)))
#d <- rbind(d, setNames(data.frame("NALM-6 RUNX1", "ETS:RUNX", melt(prop.table(table(nalm6.runx1$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d)))
d <- d[d$CountMotifs <= 5,]
d$CountMotifs <- as.factor(d$CountMotifs)
d$CL <- factor(d$CL, levels=c("AT2", "REH", "NALM-6 ER", "NALM-6 RUNX1"))
d$Motif <- factor(d$Motif, levels=c("RUNX1", "ETS", "EBF"))

d.shuffled <- data.frame("CL" = character(0), "Motif" = character(0), "CountMotifs" = character(0), "CountPeaks" = numeric(0), stringsAsFactors = F)
d.shuffled <- rbind(d.shuffled, setNames(data.frame("AT2", "RUNX1", melt(prop.table(table(at2.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)))), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("REH", "RUNX1", melt(prop.table(table(reh.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)))), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("NALM-6 ER", "RUNX1", melt(prop.table(table(nalm6.er.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)))), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("NALM-6 RUNX1", "RUNX1", melt(prop.table(table(nalm6.runx1.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)))), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("AT2", "ETS", melt(prop.table(table(at2.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("REH", "ETS", melt(prop.table(table(reh.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("NALM-6 ER", "ETS", melt(prop.table(table(nalm6.er.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("NALM-6 RUNX1", "ETS", melt(prop.table(table(nalm6.runx1.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("AT2", "EBF", melt(prop.table(table(at2.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)))), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("REH", "EBF", melt(prop.table(table(reh.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)))), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("NALM-6 ER", "EBF", melt(prop.table(table(nalm6.er.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)))), names(d.shuffled)))
d.shuffled <- rbind(d.shuffled, setNames(data.frame("NALM-6 RUNX1", "EBF", melt(prop.table(table(nalm6.runx1.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)))), names(d.shuffled)))
#d.shuffled <- rbind(d.shuffled, setNames(data.frame("AT2", "ETS:RUNX", melt(prop.table(table(at2.shuffled$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d.shuffled)))
#d.shuffled <- rbind(d.shuffled, setNames(data.frame("REH", "ETS:RUNX", melt(prop.table(table(reh.shuffled$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d.shuffled)))
#d.shuffled <- rbind(d.shuffled, setNames(data.frame("NALM-6 ER", "ETS:RUNX", melt(prop.table(table(nalm6.er.shuffled$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d.shuffled)))
#d.shuffled <- rbind(d.shuffled, setNames(data.frame("NALM-6 RUNX1", "ETS:RUNX", melt(prop.table(table(nalm6.runx1.shuffled$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)))), names(d.shuffled)))
d.shuffled <- d.shuffled[d.shuffled$CountMotifs <= 5,]
d.shuffled$CountMotifs <- as.factor(d.shuffled$CountMotifs)
d.shuffled$CL <- factor(d.shuffled$CL, levels=c("AT2", "REH", "NALM-6 ER", "NALM-6 RUNX1"))
d.shuffled$Motif <- factor(d.shuffled$Motif, levels=c("RUNX1", "ETS", "EBF"))

#pdf("/mnt/projects/fiona/results/temp.pdf", paper="a4")
ggplot(data = d, aes(x = CountMotifs, y = CountPeaks, fill = CL)) +
  geom_bar(position = position_dodge(), stat = "identity", width=0.7) +
  facet_wrap(~Motif, ncol=2) +
  coord_fixed(ratio=6) +
  labs(fill="Cell line", x="No. of motifs per peak", y="% peaks") + 
  ggtitle("No. of motifs in peaks") +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_fill_manual(values = c("red", "darkgray", "blue", "orange")) +
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.9),
        axis.title.y = element_text(size=15, vjust=1.9))

ggplot(data = d.shuffled, aes(x = CountMotifs, y = CountPeaks, fill = CL)) +
  geom_bar(position = position_dodge(), stat = "identity", width=0.7) +
  facet_wrap(~Motif, ncol=3) +
  coord_fixed(ratio=6) +
  labs(fill="Cell line", x="No. of motifs per peak", y="% peaks") + 
  ggtitle("No. of motifs in random (shuffled) peaks") +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_fill_manual(values = c("red", "darkgray", "blue", "orange")) +
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.9),
        axis.title.y = element_text(size=15, vjust=1.9))

#dev.off()

# motif distance to summit

# AT2

#pdf("/mnt/projects/fiona/results/test.pdf", paper="a4")
par(mfrow=c(2,2))

at2.runx1.dist <- at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`
at2.runx1.dist <- unlist(sapply(at2.runx1.dist[at2.runx1.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
at2.ets.dist <- at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`
at2.ets.dist <- unlist(sapply(at2.ets.dist[at2.ets.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
at2.ebf.dist <- at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`
at2.ebf.dist <- unlist(sapply(at2.ebf.dist[at2.ebf.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))

plot(density(at2.runx1.dist, adjust=1), xlim=c(-200,200), col="red", lwd=3, main="AT2 motif distance from summit", xlab="Distance from peak summit (bp)")
lines(density(at2.ets.dist, adjust=1), xlim=c(-200,200), col="darkgray", lwd=3)
lines(density(at2.ebf.dist, adjust=1), xlim=c(-200,200), col="blue", lwd=3)
abline(v=0, lty=2)
legend("topright", c("RUNX1", "ETS", "EBF"), fill=c("red", "darkgray", "blue"))

# distance b/w RUNX and ETS motif

at2.runx1.ets.dist <- at2[,c("RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit", "ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit")]
at2.runx1.ets.dist <- at2.runx1.ets.dist[at2.runx1.ets.dist[,1] != "" & at2.runx1.ets.dist[,2] != "",]
at2.runx1.ets.dist <- unlist(apply(at2.runx1.ets.dist, 1, function(x) { 
    m1 <- as.numeric(unlist(strsplit(x[1], ","))) 
    m1 <- m1[abs(m1)==min(abs(m1))] 
    m2 <- as.numeric(unlist(strsplit(x[2], ","))) 
    m2 <- m2[abs(m2)==min(abs(m2))] 
    m1-m2
}))

at2.better.runx1.ets.dist <- at2.better[,c("RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit", "ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit")]
at2.better.runx1.ets.dist <- at2.better.runx1.ets.dist[at2.better.runx1.ets.dist[,1] != "" & at2.better.runx1.ets.dist[,2] != "",]
at2.better.runx1.ets.dist <- unlist(apply(at2.better.runx1.ets.dist, 1, function(x) { 
  m1 <- as.numeric(unlist(strsplit(x[1], ","))) 
  m1 <- m1[abs(m1)==min(abs(m1))] 
  m2 <- as.numeric(unlist(strsplit(x[2], ","))) 
  m2 <- m2[abs(m2)==min(abs(m2))] 
  m1-m2
}))

plot(density(at2.runx1.ets.dist, adjust=1), xlim=c(-200,200), col="red", lwd=3, main="AT2 distance between RUNX and ETS motif", xlab="Distance between motifs (bp)")
lines(density(at2.better.runx1.ets.dist, adjust=1), xlim=c(-200,200), col="blue", lwd=3)
abline(v=0, lty=2)
legend("topright", c("all peaks", "better than RUNX1"), fill=c("red", "blue"))

# REH

reh.runx1.dist <- reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`
reh.runx1.dist <- unlist(sapply(reh.runx1.dist[reh.runx1.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
reh.ets.dist <- reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`
reh.ets.dist <- unlist(sapply(reh.ets.dist[reh.ets.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
reh.ebf.dist <- reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`
reh.ebf.dist <- unlist(sapply(reh.ebf.dist[reh.ebf.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))

plot(density(reh.runx1.dist), xlim=c(-200,200), col="red", lwd=3, main="REH motif distance from summit", xlab="Distance from peak summit (bp)")
lines(density(reh.ets.dist), xlim=c(-200,200), col="darkgray", lwd=3)
lines(density(reh.ebf.dist, adjust=1), xlim=c(-200,200), col="blue", lwd=3)
abline(v=0, lty=2)
legend("topright", c("RUNX1", "ETS", "EBF"), fill=c("red", "darkgray", "blue"))

# NALM6 ER

nalm6.er.runx1.dist <- nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`
nalm6.er.runx1.dist <- unlist(sapply(nalm6.er.runx1.dist[nalm6.er.runx1.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
nalm6.er.ets.dist <- nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`
nalm6.er.ets.dist <- unlist(sapply(nalm6.er.ets.dist[nalm6.er.ets.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
nalm6.er.ebf.dist <- nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`
nalm6.er.ebf.dist <- unlist(sapply(nalm6.er.ebf.dist[nalm6.er.ebf.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))

plot(density(nalm6.er.runx1.dist), xlim=c(-200,200), col="red", lwd=3, main="NALM-6 ER motif distance from summit", xlab="Distance from peak summit (bp)")
lines(density(nalm6.er.ets.dist), xlim=c(-200,200), col="darkgray", lwd=3)
lines(density(nalm6.er.ebf.dist, adjust=1), xlim=c(-200,200), col="blue", lwd=3)
abline(v=0, lty=2)
legend("topright", c("RUNX1", "ETS", "EBF"), fill=c("red", "darkgray", "blue"))

# NALM6 RUNX1

nalm6.runx1.runx1.dist <- nalm6.runx1$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`
nalm6.runx1.runx1.dist <- unlist(sapply(nalm6.runx1.runx1.dist[nalm6.runx1.runx1.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
nalm6.runx1.ets.dist <- nalm6.runx1$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`
nalm6.runx1.ets.dist <- unlist(sapply(nalm6.runx1.ets.dist[nalm6.runx1.ets.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))
nalm6.runx1.ebf.dist <- nalm6.runx1$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`
nalm6.runx1.ebf.dist <- unlist(sapply(nalm6.runx1.ebf.dist[nalm6.runx1.ebf.dist != ""], function(x) { d <- as.numeric(unlist(strsplit(x, ","))) ; d[abs(d)==min(abs(d))] }))

plot(density(nalm6.runx1.runx1.dist), xlim=c(-200,200), col="red", lwd=3, main="NALM-6 RUNX1 motif distance from summit", xlab="Distance from peak summit (bp)")
lines(density(nalm6.runx1.ets.dist), xlim=c(-200,200), col="darkgray", lwd=3)
lines(density(nalm6.runx1.ebf.dist), xlim=c(-200,200), col="blue", lwd=3)
abline(v=0, lty=2)
legend("topright", c("RUNX1", "ETS", "EBF"), fill=c("red", "darkgray", "blue"))

par(mfrow=c(1,1))
#dev.off()

# barchart shared vs. unique peaks

d <- data.frame("CL" = character(0), "Motif" = character(0), "Class" = character(0), "Pct" = numeric(0), stringsAsFactors = F)
d <- rbind(d, setNames(data.frame("AT2", "RUNX1", "shared", sum(1:nrow(at2) %in% o.at2.nalm6.runx1@queryHits & !is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / length(unique(o.at2.nalm6.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "ETS", "shared", sum(1:nrow(at2) %in% o.at2.nalm6.runx1@queryHits & !is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / length(unique(o.at2.nalm6.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "EBF", "shared", sum(1:nrow(at2) %in% o.at2.nalm6.runx1@queryHits & !is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / length(unique(o.at2.nalm6.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "ETS:RUNX", "shared", sum(1:nrow(at2) %in% o.at2.nalm6.runx1@queryHits & !is.na(at2$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / length(unique(o.at2.nalm6.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "RUNX1", "shared", sum(1:nrow(reh) %in% o.reh.nalm6.runx1@queryHits & !is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / length(unique(o.reh.nalm6.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "ETS", "shared", sum(1:nrow(reh) %in% o.reh.nalm6.runx1@queryHits & !is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / length(unique(o.reh.nalm6.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "EBF", "shared", sum(1:nrow(reh) %in% o.reh.nalm6.runx1@queryHits & !is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / length(unique(o.reh.nalm6.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "ETS:RUNX", "shared", sum(1:nrow(reh) %in% o.reh.nalm6.runx1@queryHits & !is.na(reh$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / length(unique(o.reh.nalm6.runx1@queryHits))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "RUNX1", "shared", sum(1:nrow(nalm6.er) %in% o.nalm6.runx1.er@subjectHits & !is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / length(unique(o.nalm6.runx1.er@subjectHits))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "ETS", "shared", sum(1:nrow(nalm6.er) %in% o.nalm6.runx1.er@subjectHits & !is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / length(unique(o.nalm6.runx1.er@subjectHits))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "EBF", "shared", sum(1:nrow(nalm6.er) %in% o.nalm6.runx1.er@subjectHits & !is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / length(unique(o.nalm6.runx1.er@subjectHits))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "ETS:RUNX", "shared", sum(1:nrow(nalm6.er) %in% o.nalm6.runx1.er@subjectHits & !is.na(nalm6.er$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / length(unique(o.nalm6.runx1.er@subjectHits))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "RUNX1", "unique", sum(!(1:nrow(at2) %in% o.at2.nalm6.runx1@queryHits) & !is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / (nrow(at2)-length(unique(o.at2.nalm6.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "ETS", "unique", sum(!(1:nrow(at2) %in% o.at2.nalm6.runx1@queryHits) & !is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / (nrow(at2)-length(unique(o.at2.nalm6.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "EBF", "unique", sum(!(1:nrow(at2) %in% o.at2.nalm6.runx1@queryHits) & !is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / (nrow(at2)-length(unique(o.at2.nalm6.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "ETS:RUNX", "unique", sum(!(1:nrow(at2) %in% o.at2.nalm6.runx1@queryHits) & !is.na(at2$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / (nrow(at2)-length(unique(o.at2.nalm6.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "RUNX1", "unique", sum(!(1:nrow(reh) %in% o.reh.nalm6.runx1@queryHits) & !is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / (nrow(reh)-length(unique(o.reh.nalm6.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "ETS", "unique", sum(!(1:nrow(reh) %in% o.reh.nalm6.runx1@queryHits) & !is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / (nrow(reh)-length(unique(o.reh.nalm6.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "EBF", "unique", sum(!(1:nrow(reh) %in% o.reh.nalm6.runx1@queryHits) & !is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / (nrow(reh)-length(unique(o.reh.nalm6.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "ETS:RUNX", "unique", sum(!(1:nrow(reh) %in% o.reh.nalm6.runx1@queryHits) & !is.na(reh$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / (nrow(reh)-length(unique(o.reh.nalm6.runx1@queryHits)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "RUNX1", "unique", sum(!(1:nrow(nalm6.er) %in% o.nalm6.runx1.er@subjectHits) & !is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / (nrow(nalm6.er)-length(unique(o.nalm6.runx1.er@subjectHits)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "ETS", "unique", sum(!(1:nrow(nalm6.er) %in% o.nalm6.runx1.er@subjectHits) & !is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / (nrow(nalm6.er)-length(unique(o.nalm6.runx1.er@subjectHits)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "EBF", "unique", sum(!(1:nrow(nalm6.er) %in% o.nalm6.runx1.er@subjectHits) & !is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / (nrow(nalm6.er)-length(unique(o.nalm6.runx1.er@subjectHits)))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "ETS:RUNX", "unique", sum(!(1:nrow(nalm6.er) %in% o.nalm6.runx1.er@subjectHits) & !is.na(nalm6.er$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / (nrow(nalm6.er)-length(unique(o.nalm6.runx1.er@subjectHits)))), names(d)))

#pdf("/mnt/projects/fiona/results/test.pdf", paper="a4")
ggplot(data = d, aes(x = CL, y = Pct, fill = Class)) +   
  geom_bar(position = position_dodge(), stat = "identity", width=0.7) +
  facet_wrap(~Motif, ncol=2) +
  labs(fill="E/R peak", x="Cell line", y="% of E/R peaks with motif") + 
  ggtitle("Percentage of E/R peaks with known motif\n") +
  theme_bw() +
  scale_fill_manual(values = c("black", "darkgray")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  coord_fixed(ratio=8) +
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.9),
        axis.title.y = element_text(size=15, vjust=1.9))
#dev.off()

# barchart E/R better vs. worse

d <- data.frame("CL" = character(0), "Motif" = character(0), "Class" = character(0), "Pct" = numeric(0), stringsAsFactors = F)
d <- rbind(d, setNames(data.frame("AT2", "RUNX1", "unique or better", sum(at2$runx1_overlap %in% c("shared_better","unique") & !is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / sum(at2$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "ETS", "unique or better", sum(at2$runx1_overlap %in% c("shared_better","unique") & !is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / sum(at2$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "EBF", "unique or better", sum(at2$runx1_overlap %in% c("shared_better","unique") & !is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / sum(at2$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "ETS:RUNX", "unique or better", sum(at2$runx1_overlap %in% c("shared_better","unique") & !is.na(at2$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / sum(at2$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "RUNX1", "unique or better", sum(reh$runx1_overlap %in% c("shared_better","unique") & !is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / sum(reh$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "ETS", "unique or better", sum(reh$runx1_overlap %in% c("shared_better","unique") & !is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / sum(reh$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "EBF", "unique or better", sum(reh$runx1_overlap %in% c("shared_better","unique") & !is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / sum(reh$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "ETS:RUNX", "unique or better", sum(reh$runx1_overlap %in% c("shared_better","unique") & !is.na(reh$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / sum(reh$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "RUNX1", "unique or better", sum(nalm6.er$runx1_overlap %in% c("shared_better","unique") & !is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / sum(nalm6.er$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "ETS", "unique or better", sum(nalm6.er$runx1_overlap %in% c("shared_better","unique") & !is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / sum(nalm6.er$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "EBF", "unique or better", sum(nalm6.er$runx1_overlap %in% c("shared_better","unique") & !is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / sum(nalm6.er$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "ETS:RUNX", "unique or better", sum(nalm6.er$runx1_overlap %in% c("shared_better","unique") & !is.na(nalm6.er$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / sum(nalm6.er$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "RUNX1", "worse", sum(!at2$runx1_overlap %in% c("shared_better","unique") & !is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / sum(!at2$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "ETS", "worse", sum(!at2$runx1_overlap %in% c("shared_better","unique") & !is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / sum(!at2$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "EBF", "worse", sum(!at2$runx1_overlap %in% c("shared_better","unique") & !is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / sum(!at2$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("AT2", "ETS:RUNX", "worse", sum(!at2$runx1_overlap %in% c("shared_better","unique") & !is.na(at2$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / sum(!at2$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "RUNX1", "worse", sum(!reh$runx1_overlap %in% c("shared_better","unique") & !is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / sum(!reh$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "ETS", "worse", sum(!reh$runx1_overlap %in% c("shared_better","unique") & !is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / sum(!reh$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "EBF", "worse", sum(!reh$runx1_overlap %in% c("shared_better","unique") & !is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / sum(!reh$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("REH", "ETS:RUNX", "worse", sum(!reh$runx1_overlap %in% c("shared_better","unique") & !is.na(reh$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / sum(!reh$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "RUNX1", "worse", sum(!nalm6.er$runx1_overlap %in% c("shared_better","unique") & !is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)) / sum(!nalm6.er$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "ETS", "worse", sum(!nalm6.er$runx1_overlap %in% c("shared_better","unique") & !is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / sum(!nalm6.er$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "EBF", "worse", sum(!nalm6.er$runx1_overlap %in% c("shared_better","unique") & !is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)) / sum(!nalm6.er$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6", "ETS:RUNX", "worse", sum(!nalm6.er$runx1_overlap %in% c("shared_better","unique") & !is.na(nalm6.er$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)) / sum(!nalm6.er$runx1_overlap %in% c("shared_better","unique"))), names(d)))
d$Class <- factor(d$Class, levels=c("worse", "unique or better"))

#pdf("/mnt/projects/fiona/results/test.pdf", paper="a4")
ggplot(data = d, aes(x = CL, y = Pct, fill = Class)) +   
  geom_bar(position = position_dodge(), stat = "identity", width=0.7) +
  facet_wrap(~Motif, ncol=2) +
  labs(fill="E/R peak", x="Cell line", y="% of E/R peaks with motif") + 
  ggtitle("Percentage of E/R peaks with known motif\n") +
  theme_bw() +
  scale_fill_manual(values = c("black", "darkgray")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  coord_fixed(ratio=8) +
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.9),
        axis.title.y = element_text(size=15, vjust=1.9))
#dev.off()

# motif overlap AT2

library(Vennerable)

idcol <- which(grepl("PeakID", names(at2)))
peak.sets.at2 <- list(
  "RUNX1" = at2[,idcol][!is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)],
  "ETS" = at2[,idcol][!is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)],
  "EBF" = at2[,idcol][!is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)]
)

idcol <- which(grepl("PeakID", names(at2.shuffled)))
peak.sets.at2.shuffled <- list(
  "RUNX1" = at2.shuffled[,idcol][!is.na(at2.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)],
  "ETS" = at2.shuffled[,idcol][!is.na(at2.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)],
  "EBF" = at2.shuffled[,idcol][!is.na(at2.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)]
)

#pdf("/mnt/projects/fiona/results/temp.pdf", paper="a4")

venn <- compute.Venn(Venn(peak.sets.at2), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("AT2 Motif Overlap", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

venn <- compute.Venn(Venn(peak.sets.at2.shuffled), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("AT2 Motif Overlap (shuffled peaks)", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

#dev.off()

# motif overlap REH

idcol <- which(grepl("PeakID", names(reh)))
peak.sets.reh <- list(
  "RUNX1" = reh[,idcol][!is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)],
  "ETS" = reh[,idcol][!is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)],
  "EBF" = reh[,idcol][!is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)]
)

idcol <- which(grepl("PeakID", names(reh.shuffled)))
peak.sets.reh.shuffled <- list(
  "RUNX1" = reh.shuffled[,idcol][!is.na(reh.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)],
  "ETS" = reh.shuffled[,idcol][!is.na(reh.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)],
  "EBF" = reh.shuffled[,idcol][!is.na(reh.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)]
)

#pdf("/mnt/projects/fiona/results/temp.pdf", paper="a4")

venn <- compute.Venn(Venn(peak.sets.reh), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("REH Motif Overlap", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

venn <- compute.Venn(Venn(peak.sets.reh.shuffled), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("REH Motif Overlap (shuffled peaks)", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

#dev.off()

# motif overlap NALM-6 RUNX1

idcol <- which(grepl("PeakID", names(nalm6.runx1)))
peak.sets.nalm6.runx1 <- list(
  "RUNX1" = nalm6.runx1[,idcol][!is.na(nalm6.runx1$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)],
  "ETS" = nalm6.runx1[,idcol][!is.na(nalm6.runx1$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)],
  "EBF" = nalm6.runx1[,idcol][!is.na(nalm6.runx1$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)]
)

idcol <- which(grepl("PeakID", names(nalm6.runx1.shuffled)))
peak.sets.nalm6.runx1.shuffled <- list(
  "RUNX1" = nalm6.runx1.shuffled[,idcol][!is.na(nalm6.runx1.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)],
  "ETS" = nalm6.runx1.shuffled[,idcol][!is.na(nalm6.runx1.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)],
  "EBF" = nalm6.runx1.shuffled[,idcol][!is.na(nalm6.runx1.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)]
)

#pdf("/mnt/projects/fiona/results/temp.pdf", paper="a4")

venn <- compute.Venn(Venn(peak.sets.nalm6.runx1), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("NALM-6 RUNX1 Motif Overlap", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

venn <- compute.Venn(Venn(peak.sets.nalm6.runx1.shuffled), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("NALM-6 RUNX1 Motif Overlap (shuffled peaks)", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

#dev.off()

# motif overlap NALM-6 ER

idcol <- which(grepl("PeakID", names(nalm6.er)))
peak.sets.nalm6.er <- list(
  "RUNX1" = nalm6.er[,idcol][!is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)],
  "ETS" = nalm6.er[,idcol][!is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)],
  "EBF" = nalm6.er[,idcol][!is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)]
)

idcol <- which(grepl("PeakID", names(nalm6.er.shuffled)))
peak.sets.nalm6.er.shuffled <- list(
  "RUNX1" = nalm6.er.shuffled[,idcol][!is.na(nalm6.er.shuffled$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)],
  "ETS" = nalm6.er.shuffled[,idcol][!is.na(nalm6.er.shuffled$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)],
  "EBF" = nalm6.er.shuffled[,idcol][!is.na(nalm6.er.shuffled$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)]
)

#pdf("/mnt/projects/fiona/results/temp.pdf", paper="a4")

venn <- compute.Venn(Venn(peak.sets.nalm6.er), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("NALM-6 ER Motif Overlap", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

venn <- compute.Venn(Venn(peak.sets.nalm6.er.shuffled), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("NALM-6 ER Motif Overlap (shuffled peaks)", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

#dev.off()

# motif composition shared RUNX1 vs. E/R unique binding sites

#pdf("/mnt/projects/fiona/results/test.pdf", paper="a4")

# AT2

idcol <- which(grepl("PeakID", names(at2)))

print(sprintf("Number of shared peaks in AT2: %d", length(unique(o.at2.nalm6.runx1@queryHits))))

peak.sets.at2.shared <- list(
  "RUNX1" = at2[intersect(o.at2.nalm6.runx1@queryHits, which(!is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`))),idcol],
  "ETS" = at2[intersect(o.at2.nalm6.runx1@queryHits, which(!is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`))),idcol],
  "EBF" = at2[intersect(o.at2.nalm6.runx1@queryHits, which(!is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`))),idcol]
)

venn <- compute.Venn(Venn(peak.sets.at2.shared), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("AT2 motif overlap shared E/R peaks", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))


print(sprintf("Number of unique peaks in AT2: %d", nrow(at2)-length(unique(o.at2.nalm6.runx1@queryHits))))

peak.sets.at2.unique <- list(
  "RUNX1" = at2[intersect((1:nrow(at2))[-o.at2.nalm6.runx1@queryHits], which(!is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`))),idcol],
  "ETS" = at2[intersect((1:nrow(at2))[-o.at2.nalm6.runx1@queryHits], which(!is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`))),idcol],
  "EBF" = at2[intersect((1:nrow(at2))[-o.at2.nalm6.runx1@queryHits], which(!is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`))),idcol]
)

venn <- compute.Venn(Venn(peak.sets.at2.unique), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("AT2 motif overlap unique E/R peaks", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

# REH

idcol <- which(grepl("PeakID", names(reh)))

print(sprintf("Number of shared peaks in REH: %d", length(unique(o.reh.nalm6.runx1@queryHits))))

peak.sets.reh.shared <- list(
  "RUNX1" = reh[intersect(o.reh.nalm6.runx1@queryHits, which(!is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`))),idcol],
  "ETS" = reh[intersect(o.reh.nalm6.runx1@queryHits, which(!is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`))),idcol],
  "EBF" = reh[intersect(o.reh.nalm6.runx1@queryHits, which(!is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`))),idcol]
)

venn <- compute.Venn(Venn(peak.sets.reh.shared), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("REH motif overlap shared E/R peaks", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

print(sprintf("Number of unique peaks in REH: %d", nrow(reh)-length(unique(o.reh.nalm6.runx1@queryHits))))

peak.sets.reh.unique <- list(
  "RUNX1" = reh[intersect((1:nrow(reh))[-o.reh.nalm6.runx1@queryHits], which(!is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`))),idcol],
  "ETS" = reh[intersect((1:nrow(reh))[-o.reh.nalm6.runx1@queryHits], which(!is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`))),idcol],
  "EBF" = reh[intersect((1:nrow(reh))[-o.reh.nalm6.runx1@queryHits], which(!is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`))),idcol]
)

venn <- compute.Venn(Venn(peak.sets.reh.unique), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("REH motif overlap unique E/R peaks", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

# NALM-6

idcol <- which(grepl("PeakID", names(nalm6.er)))

print(sprintf("Number of shared peaks in NALM-6 E/R: %d", length(unique(o.nalm6.runx1.er@subjectHits))))

peak.sets.nalm6.er.shared <- list(
  "RUNX1" = nalm6.er[intersect(o.nalm6.runx1.er@subjectHits, which(!is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`))),idcol],
  "ETS" = nalm6.er[intersect(o.nalm6.runx1.er@subjectHits, which(!is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`))),idcol],
  "EBF" = nalm6.er[intersect(o.nalm6.runx1.er@subjectHits, which(!is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`))),idcol]
)

venn <- compute.Venn(Venn(peak.sets.nalm6.er.shared), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("NALM-6 motif overlap shared E/R peaks", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

print(sprintf("Number of unique peaks in NALM-6 E/R: %d", nrow(nalm6.er)-length(unique(o.nalm6.runx1.er@subjectHits))))

peak.sets.nalm6.er.unique <- list(
  "RUNX1" = nalm6.er[intersect((1:nrow(nalm6.er))[-o.nalm6.runx1.er@subjectHits], which(!is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`))),idcol],
  "ETS" = nalm6.er[intersect((1:nrow(nalm6.er))[-o.nalm6.runx1.er@subjectHits], which(!is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`))),idcol],
  "EBF" = nalm6.er[intersect((1:nrow(nalm6.er))[-o.nalm6.runx1.er@subjectHits], which(!is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`))),idcol]
)

venn <- compute.Venn(Venn(peak.sets.nalm6.er.unique), doWeights = TRUE, type = "circles")
gp <- VennThemes(venn, colourAlgorithm = "signature")
grid.newpage()
plot(venn, gpList = gp)
grid.text("NALM-6 motif overlap unique E/R peaks", vp = viewport(x = 0.5, y = 0.98, w=unit(1, "npc"), h=unit(1, "npc")))

# motif frequency proximal vs. distal

#pdf("/mnt/projects/fiona/results/test.pdf", paper="a4")

library(vcd)

# AT2 RUNX1

in.promoter <- at2$`Distance to TSS` >= -5000 & at2$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & at2$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(at2$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("RUNX+", "RUNX-")))
)
mosaic(t, pop=F, main=sprintf("AT2 p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# AT2 ETS

in.promoter <- at2$`Distance to TSS` >= -5000 & at2$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & at2$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(at2$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("ETS+", "ETS-")))
)
mosaic(t, pop=F, main=sprintf("AT2 p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# AT2 EBF

in.promoter <- at2$`Distance to TSS` >= -5000 & at2$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & at2$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(at2$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("EBF+", "EBF-")))
)
mosaic(t, pop=F, main=sprintf("AT2 p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# AT2 ETS:RUNX

in.promoter <- at2$`Distance to TSS` >= -5000 & at2$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & at2$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(at2$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("ETS:RUNX+", "ETS:RUNX-")))
)
mosaic(t, pop=F, main=sprintf("AT2 p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# REH RUNX1

in.promoter <- reh$`Distance to TSS` >= -5000 & reh$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & reh$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(reh$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("RUNX+", "RUNX-")))
)
mosaic(t, pop=F, main=sprintf("REH p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# REH ETS

in.promoter <- reh$`Distance to TSS` >= -5000 & reh$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & reh$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(reh$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("ETS+", "ETS-")))
)
mosaic(t, pop=F, main=sprintf("REH p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# REH EBF

in.promoter <- reh$`Distance to TSS` >= -5000 & reh$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & reh$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(reh$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("EBF+", "EBF-")))
)
mosaic(t, pop=F, main=sprintf("REH p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# REH ETS:RUNX

in.promoter <- reh$`Distance to TSS` >= -5000 & reh$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & reh$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(reh$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("ETS:RUNX+", "ETS:RUNX-")))
)
mosaic(t, pop=F, main=sprintf("REH p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# NALM-6 ER RUNX1

in.promoter <- nalm6.er$`Distance to TSS` >= -5000 & nalm6.er$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & nalm6.er$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(nalm6.er$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("RUNX+", "RUNX-")))
)
mosaic(t, pop=F, main=sprintf("NALM-6 ER p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# NALM-6 ER ETS

in.promoter <- nalm6.er$`Distance to TSS` >= -5000 & nalm6.er$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & nalm6.er$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(nalm6.er$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("ETS+", "ETS-")))
)
mosaic(t, pop=F, main=sprintf("NALM-6 ER p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# NALM-6 ER EBF

in.promoter <- nalm6.er$`Distance to TSS` >= -5000 & nalm6.er$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & nalm6.er$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(nalm6.er$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("EBF+", "EBF-")))
)
mosaic(t, pop=F, main=sprintf("NALM-6 ER p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# NALM-6 ER ETS:RUNX

in.promoter <- nalm6.er$`Distance to TSS` >= -5000 & nalm6.er$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & nalm6.er$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(nalm6.er$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("ETS:RUNX+", "ETS:RUNX-")))
)
mosaic(t, pop=F, main=sprintf("NALM-6 ER p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# NALM-6 RUNX1 RUNX

in.promoter <- nalm6.runx1$`Distance to TSS` >= -5000 & nalm6.runx1$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & nalm6.runx1$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(nalm6.runx1$`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("RUNX+", "RUNX-")))
)
mosaic(t, pop=F, main=sprintf("NALM-6 RUNX1 p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# NALM-6 RUNX1 ETS

in.promoter <- nalm6.runx1$`Distance to TSS` >= -5000 & nalm6.runx1$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & nalm6.runx1$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(nalm6.runx1$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("ETS+", "ETS-")))
)
mosaic(t, pop=F, main=sprintf("NALM-6 RUNX1 p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# NALM-6 RUNX1 EBF

in.promoter <- nalm6.runx1$`Distance to TSS` >= -5000 & nalm6.runx1$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & nalm6.runx1$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(nalm6.runx1$`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("EBF+", "EBF-")))
)
mosaic(t, pop=F, main=sprintf("NALM-6 RUNX1 p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# NALM-6 RUNX1 ETS:RUNX

in.promoter <- nalm6.runx1$`Distance to TSS` >= -5000 & nalm6.runx1$`Distance to TSS` <= 1000
in.enhancer <- !in.promoter & nalm6.runx1$overlaps_enhancer_in_celllines != ""
has.motif <- !is.na(nalm6.runx1$`ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq(GSE17954)/Homer No. motifs`)
t <- as.table(matrix(c(sum(in.promoter & has.motif), 
                       sum(in.enhancer & has.motif), 
                       sum(!in.promoter & !in.enhancer & has.motif), 
                       sum(in.promoter & !has.motif), 
                       sum(in.enhancer & !has.motif), 
                       sum(!in.promoter & !in.enhancer & !has.motif)),
                     nrow = 3,
                     dimnames = list('Region'=c("Promoter", "Enhancer", "Other"), 
                                     'Motif'=c("ETS:RUNX+", "ETS:RUNX-")))
)
mosaic(t, pop=F, main=sprintf("NALM-6 RUNX1 p=%.2g", fisher.test(t, workspace = 2e9)$p.value))
labeling_cells(text=t)(t)

# peak score stratified by motif occurrence

windowSize <- 100

d <- with(at2, data.frame("CL" = "AT2", "ID" = at2[,2], "Score" = `Peak Score`, "width" = End-Start, "TSSdist" = `Distance to TSS`,
                          "numRUNX" = `RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`,
                          "numETS" = `ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`,
                          "numEBF" = `EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`,
                          "numRUNXSummit" = suppressWarnings(sapply(strsplit(`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`, ","), function(x) sum(abs(as.numeric(x)) <= windowSize))),
                          "numETSSummit" = suppressWarnings(sapply(strsplit(`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`, ","), function(x) sum(abs(as.numeric(x)) <= windowSize))),
                          "numEBFSummit" = suppressWarnings(sapply(strsplit(`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`, ","), function(x) sum(abs(as.numeric(x)) <= windowSize))),
                          "distRUNX" = suppressWarnings(sapply(strsplit(`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`, ","), function(x) min(abs(as.numeric(x)), na.rm=T))),
                          "distETS" = suppressWarnings(sapply(strsplit(`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`, ","), function(x) min(abs(as.numeric(x)), na.rm=T))),
                          "distEBF" = suppressWarnings(sapply(strsplit(`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`, ","), function(x) min(abs(as.numeric(x)), na.rm=T)))
))
d <- rbind(d, with(reh, data.frame("CL" = "REH", "ID" = reh[,2], "Score" = `Peak Score`, "width" = End-Start, "TSSdist" = `Distance to TSS`, 
                          "numRUNX" = `RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`,
                          "numETS" = `ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`,
                          "numEBF" = `EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`,
                          "numRUNXSummit" = suppressWarnings(sapply(strsplit(`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`, ","), function(x) sum(abs(as.numeric(x)) <= windowSize))),
                          "numETSSummit" = suppressWarnings(sapply(strsplit(`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`, ","), function(x) sum(abs(as.numeric(x)) <= windowSize))),
                          "numEBFSummit" = suppressWarnings(sapply(strsplit(`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`, ","), function(x) sum(abs(as.numeric(x)) <= windowSize))),
                          "distRUNX" = suppressWarnings(sapply(strsplit(`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`, ","), function(x) min(abs(as.numeric(x)), na.rm=T))),
                          "distETS" = suppressWarnings(sapply(strsplit(`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`, ","), function(x) min(abs(as.numeric(x)), na.rm=T))),
                          "distEBF" = suppressWarnings(sapply(strsplit(`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`, ","), function(x) min(abs(as.numeric(x)), na.rm=T)))
)))
d <- rbind(d, with(nalm6.er, data.frame("CL" = "NALM-6 ER", "ID" = nalm6.er[,2], "Score" = `Peak Score`, "width" = End-Start, "TSSdist" = `Distance to TSS`, 
                                   "numRUNX" = `RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`,
                                   "numETS" = `ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`,
                                   "numEBF" = `EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`,
                                   "numRUNXSummit" = suppressWarnings(sapply(strsplit(`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`, ","), function(x) sum(abs(as.numeric(x)) <= windowSize))),
                                   "numETSSummit" = suppressWarnings(sapply(strsplit(`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`, ","), function(x) sum(abs(as.numeric(x)) <= windowSize))),
                                   "numEBFSummit" = suppressWarnings(sapply(strsplit(`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`, ","), function(x) sum(abs(as.numeric(x)) <= windowSize))),
                                   "distRUNX" = suppressWarnings(sapply(strsplit(`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`, ","), function(x) min(abs(as.numeric(x)), na.rm=T))),
                                   "distETS" = suppressWarnings(sapply(strsplit(`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`, ","), function(x) min(abs(as.numeric(x)), na.rm=T))),
                                   "distEBF" = suppressWarnings(sapply(strsplit(`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`, ","), function(x) min(abs(as.numeric(x)), na.rm=T)))
)))
d <- rbind(d, with(nalm6.runx1, data.frame("CL" = "NALM-6 RUNX1", "ID" = nalm6.runx1[,2], "Score" = `Peak Score`, "width" = End-Start, "TSSdist" = `Distance to TSS`,
                                        "numRUNX" = `RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer No. motifs`,
                                        "numETS" = `ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer No. motifs`,
                                        "numEBF" = `EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer No. motifs`,
                                        "numRUNXSummit" = suppressWarnings(sapply(strsplit(`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`, ","), function(x) sum(abs(as.numeric(x)) <= windowSize))),
                                        "numETSSummit" = suppressWarnings(sapply(strsplit(`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`, ","), function(x) sum(abs(as.numeric(x)) <= windowSize))),
                                        "numEBFSummit" = suppressWarnings(sapply(strsplit(`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`, ","), function(x) sum(abs(as.numeric(x)) <= windowSize))),
                                        "distRUNX" = suppressWarnings(sapply(strsplit(`RUNX1(Runt)/Jurkat-RUNX1-ChIP-Seq(GSE29180)/Homer Distance From Summit`, ","), function(x) min(abs(as.numeric(x)), na.rm=T))),
                                        "distETS" = suppressWarnings(sapply(strsplit(`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer Distance From Summit`, ","), function(x) min(abs(as.numeric(x)), na.rm=T))),
                                        "distEBF" = suppressWarnings(sapply(strsplit(`EBF(EBF)/proBcell-EBF-ChIP-Seq(GSE21978)/Homer Distance From Summit`, ","), function(x) min(abs(as.numeric(x)), na.rm=T)))
)))

d$numRUNXorETSSummit <- d$numRUNXSummit + d$numETSSummit
d$numRUNXSummit[d$numRUNXSummit >= 3] <- "3+"
d$numETSSummit[d$numETSSummit >= 3] <- "3+"
d$numEBFSummit[d$numEBFSummit >= 3] <- "3+"
d$numRUNXorETSSummit[d$numRUNXorETSSummit >= 3] <- "3+"

d$category[d$numRUNXSummit == 0 & d$numETSSummit == 0 & d$numEBFSummit == 0] <- "None"
d$category[d$numRUNXSummit > 0 & d$numETSSummit == 0 & d$numEBFSummit == 0] <- "only RUNX"
d$category[d$numRUNXSummit == 0 & d$numETSSummit > 0 & d$numEBFSummit == 0] <- "only ETS"
d$category[d$numRUNXSummit == 0 & d$numETSSummit == 0 & d$numEBFSummit > 0] <- "only EBF"
d$category[d$numRUNXSummit > 0 & d$numETSSummit > 0 & d$numEBFSummit == 0] <- "RUNX and ETS"
d$category[d$numRUNXSummit > 0 & d$numETSSummit == 0 & d$numEBFSummit > 0] <- "RUNX and EBF"
d$category[d$numRUNXSummit == 0 & d$numETSSummit > 0 & d$numEBFSummit > 0] <- "ETS and EBF"
d$category[d$numRUNXSummit > 0 & d$numETSSummit > 0 & d$numEBFSummit > 0] <- "All three"
d$category <- factor(d$category, levels=c("All three", "RUNX and ETS", "RUNX and EBF", "only RUNX", "ETS and EBF", "only ETS", "only EBF", "None"))
d$CL <- factor(d$CL, levels=c("AT2", "REH", "NALM-6 ER", "NALM-6 RUNX1"))

#pdf("/mnt/projects/fiona/results/test.pdf", paper="a4")

ggplot(data = d, aes(x = category, y = log(Score, 2))) +
  geom_violin(fill="#FFCDFF") +
  geom_boxplot(width=.15, size=0.3, aes(fill=CL), outlier.colour = NA, show.legend = F) +
  scale_fill_manual(values = rep("white", length(levels(d$CL)))) +
  facet_wrap(~CL, ncol = 2) +
  labs(x=sprintf("Motifs within %d bp of peak summit", windowSize), y="Log2 Peak Score") + 
  ggtitle("Peak score by motif composition near summit") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=90, hjust=1, size=10),
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.3),
        axis.title.y = element_text(size=15, vjust=1.9)) +
  coord_fixed(ratio=0.8)

ggplot(data = d, aes(x = factor(numRUNXSummit), y = log(Score, 2))) +   
  geom_violin(fill="#FFCDFF") +
  geom_boxplot(width=.15, size=0.3, aes(fill=CL), outlier.colour = NA, show.legend = F) +
  scale_fill_manual(values = rep("white", length(levels(d$CL)))) +
  geom_smooth(aes(group=1), method="gam", size=1.5) +
  facet_wrap(~CL, ncol = 2) +
  labs(x=sprintf("Number of motifs within %d bp of peak summit", windowSize), y="Log2 Peak Score") + 
  ggtitle("Peak score by number of RUNX motifs near summit") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=10),
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.3),
        axis.title.y = element_text(size=15, vjust=1.9)) +
  coord_fixed(ratio=0.5)

ggplot(data = d, aes(x = factor(numETSSummit), y = log(Score, 2))) +   
  geom_violin(fill="#FFCDFF") +
  geom_boxplot(width=.15, size=0.3, aes(fill=CL), outlier.colour = NA, show.legend = F) +
  scale_fill_manual(values = rep("white", length(levels(d$CL)))) +
  geom_smooth(aes(group=1), method="gam", size=1.5) +
  facet_wrap(~CL, ncol = 2) +
  labs(x=sprintf("Number of motifs within %d bp of peak summit", windowSize), y="Log2 Peak Score") + 
  ggtitle("Peak score by number of ETS motifs near summit") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=90, hjust=1, size=10),
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.3),
        axis.title.y = element_text(size=15, vjust=1.9)) +
  coord_fixed(ratio=0.5)

ggplot(data = d, aes(x = factor(numRUNXorETSSummit), y = log(Score, 2))) +   
  geom_violin(fill="#FFCDFF") +
  geom_boxplot(width=.15, size=0.3, aes(fill=CL), outlier.colour = NA, show.legend = F) +
  scale_fill_manual(values = rep("white", length(levels(d$CL)))) +
  geom_smooth(aes(group=1), method="gam", size=1.5) +
  facet_wrap(~CL, ncol = 2) +
  labs(x=sprintf("Number of motifs within %d bp of peak summit", windowSize), y="Log2 Peak Score") + 
  ggtitle("Peak score by number of RUNX or ETS motifs near summit") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=90, hjust=1, size=10),
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.3),
        axis.title.y = element_text(size=15, vjust=1.9)) +
  coord_fixed(ratio=0.5)

ggplot(data = d, aes(x = factor(numEBFSummit), y = log(Score, 2))) +   
  geom_violin(fill="#FFCDFF") +
  geom_boxplot(width=.15, size=0.3, aes(fill=CL), outlier.colour = NA, show.legend = F) +
  scale_fill_manual(values = rep("white", length(levels(d$CL)))) +
  geom_smooth(aes(group=1), method="gam", size=1.5) +
  facet_wrap(~CL, ncol = 2) +
  labs(x=sprintf("Number of motifs within %d bp of peak summit", windowSize), y="Log2 Peak Score") + 
  ggtitle("Peak score by number of EBF motifs near summit") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=90, hjust=1, size=10),
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.3),
        axis.title.y = element_text(size=15, vjust=1.9)) +
  coord_fixed(ratio=0.5)

dev.off()

# distance between co-occuring motifs

