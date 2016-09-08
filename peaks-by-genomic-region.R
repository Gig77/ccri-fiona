source("/mnt/projects/fiona/scripts/read-peaks.R")

getSimpleAnnotation <- function(df) {
  gsub("(.*?)( \\(|-[0-9]).*", "\\1", df$Annotation)
}

d <- data.frame("CL" = character(0), "Region" = character(0), "Frequency" = numeric(0), stringsAsFactors = F)
d <- rbind(d, setNames(data.frame("AT2", melt(prop.table(table(getSimpleAnnotation(at2))))), names(d)))
d$Region <- factor(d$Region, levels=d$Region[order(d$Frequency, decreasing = T)])
d <- rbind(d, setNames(data.frame("REH", melt(prop.table(table(getSimpleAnnotation(reh))))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6 E/R", melt(prop.table(table(getSimpleAnnotation(nalm6.er))))), names(d)))
d <- rbind(d, setNames(data.frame("NALM-6 RUNX1", melt(prop.table(table(getSimpleAnnotation(nalm6.runx1))))), names(d)))
d$CL <- factor(d$CL, levels=c("AT2", "REH", "NALM-6 E/R", "NALM-6 RUNX1"))

pdf("/mnt/projects/fiona/results/peaks-by-genomic-region.pdf", width=9, height=4)
print(ggplot(data = d, aes(x = Region, y = Frequency, fill = CL)) +
  geom_bar(position = position_dodge(), stat = "identity", width=0.7) +
  labs(fill="ChIP", x="Region", y="Frequency") + 
  ggtitle("Peaks by genomic region") +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_fill_manual(values = c("red", "blue", "darkgray", "black")) +
  theme(axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5), 
        axis.text.y = element_text(size=10), 
        axis.title.x = element_text(size=15, vjust=-0.5),
        axis.title.y = element_text(size=15, vjust=1.4)))
dev.off()
