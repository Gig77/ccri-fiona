export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PROJECT_HOME=/mnt/projects/fiona
DOCKER=docker run --rm --net=host -e DOCKER_UID=$$(id -u) -e DOCKER_UNAME=$$(id -un) -e DOCKER_GID=$$(id -g) -e DOCKER_GNAME=$$(id -gn) -e DOCKER_HOME=$$HOME -v /home:/home -v /data_synology:/data_synology -v /data:/data -v /data2:/data2 -v /home/cf/fiona/results/current:$(PROJECT_HOME)/results -v /home/cf/fiona/data:$(PROJECT_HOME)/data -w $$(pwd)
MACS2=$(DOCKER) biowaste:5000/ccri/macs-2.0.9 macs2
HOMER=$(DOCKER) biowaste:5000/ccri/homer-4.7
BWA=/data_synology/software/bwa-0.7.12/bwa
SAMTOOLS=/data_synology/software/samtools-1.1/samtools

SAMPLES_RNASEQ=32232_CGATGT_C80BJANXX_3_20150925B_20150925 \
			   32233_TGACCA_C80BJANXX_3_20150925B_20150925 \
			   32234_ACAGTG_C80BJANXX_3_20150925B_20150925 \
			   32235_GCCAAT_C80BJANXX_3_20150925B_20150925 \
			   32236_CAGATC_C80BJANXX_3_20150925B_20150925 \
			   32237_CTTGTA_C80BJANXX_3_20150925B_20150925

SAMPLES_CHIPSEQ=32238_CGATGT_C80K5ANXX_6_20150930B_20150930 \
				32239_TGACCA_C80K5ANXX_6_20150930B_20150930 \
				32240_ACAGTG_C80K5ANXX_6_20150930B_20150930 \
				32241_GCCAAT_C80K5ANXX_6_20150930B_20150930 \
				32242_CAGATC_C80K5ANXX_6_20150930B_20150930 \
				32243_CTTGTA_C80K5ANXX_6_20150930B_20150930

all: fastqc bwa homer

#-----------------------------------------------------------------------------------------

.PHONY: download
download:
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32232_CGATGT_C80BJANXX_3_20150925B_20150925.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32233_TGACCA_C80BJANXX_3_20150925B_20150925.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32234_ACAGTG_C80BJANXX_3_20150925B_20150925.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32235_GCCAAT_C80BJANXX_3_20150925B_20150925.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32236_CAGATC_C80BJANXX_3_20150925B_20150925.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32237_CTTGTA_C80BJANXX_3_20150925B_20150925.bam

	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32238_CGATGT_C80K5ANXX_6_20150930B_20150930.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32239_TGACCA_C80K5ANXX_6_20150930B_20150930.bam
	wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32240_ACAGTG_C80K5ANXX_6_20150930B_20150930.bam
	wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32241_GCCAAT_C80K5ANXX_6_20150930B_20150930.bam
	wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32242_CAGATC_C80K5ANXX_6_20150930B_20150930.bam
	wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32243_CTTGTA_C80K5ANXX_6_20150930B_20150930.bam

# --------------------------------------------------------------------------------
# FASTQC
# --------------------------------------------------------------------------------

.PHONY: fastqc
fastqc: $(foreach S, $(SAMPLES_CHIPSEQ), fastqc/$S_fastqc.html)
	
fastqc/%_fastqc.html: $(PROJECT_HOME)/data/bam/chipseq/%.bam
	mkdir -p fastqc/$*.part
	/data_synology/software/FastQC-0.11.2/fastqc -o fastqc/$*.part -f bam $^
	mv fastqc/$*.part/* fastqc
	rmdir fastqc/$*.part

# --------------------------------------------------------------------------------
# samtools flagstat
# --------------------------------------------------------------------------------

.PHONY: flagstat
flagstat: $(foreach S, $(SAMPLES_CHIPSEQ), qc/$S.samtools.flagstat)

qc/%.samtools.flagstat: bwa/%.bwa.sorted.bam
	mkdir -p qc
	$(SAMTOOLS) flagstat $< 2>&1 1>$@.part | $(LOG)
	mv $@.part $@

# --------------------------------------------------------------------------------
# alignment, sorting, duplication marking, filtering
# --------------------------------------------------------------------------------

.PHONY: bwa
bwa: $(foreach S, $(SAMPLES_CHIPSEQ), bwa/$S.bwa.sorted.bam.bai)

bwa/%.bwa.sorted.bam: $(PROJECT_HOME)/data/bam/chipseq/%.bam
	mkdir -p bwa
	$(BWA) aln -b -t 10 /mnt/projects/generic/data/broad/human_g1k_v37.fasta $< \
		| $(BWA) samse /mnt/projects/generic/data/broad/human_g1k_v37.fasta - $< \
		| $(SAMTOOLS) view -Shb -@ 2 - \
		| $(SAMTOOLS) sort -T bwa/$* -o $@.part -O bam -@ 2 -m 2G \
		| 2>&1 | $(LOG)		
	mv $@.part $@

bwa/%.bwa.sorted.dupmarked.bam: bwa/%.bwa.sorted.bam
	mkdir -p picard
	java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar /data_synology/software/picard-tools-1.114/MarkDuplicates.jar \
		INPUT=$< \
		OUTPUT=$@.part \
		METRICS_FILE=picard/$*.mark_duplicates_metrics \
		VALIDATION_STRINGENCY=LENIENT \
		2>&1 1>$@.part | $(LOG)
	mv $@.part $@

bwa/%.bwa.sorted.filtered.bam: bwa/%.bwa.sorted.bam
	$(SAMTOOLS) view -bhq 30 $< > $@.part
	mv $@.part $@
	 
bwa/%.bwa.sorted.bam.bai: bwa/%.bwa.sorted.bam
	rm -f $@
	$(SAMTOOLS) index $^ $@.part 2>&1 | $(LOG)
	
# --------------------------------------------------------------------------------
# peak calling (MACS)
# --------------------------------------------------------------------------------

.PHONY: macs
macs: macs/runx1_peaks.bed macs/er_peaks.bed macs/rhd_peaks.bed

macs/runx1_peaks.bed: bwa/32243_CTTGTA_C80K5ANXX_6_20150930B_20150930.bwa.sorted.filtered.bam bwa/32242_CAGATC_C80K5ANXX_6_20150930B_20150930.bwa.sorted.filtered.bam
	mkdir -p macs
	#WD=$$(pwd) && cd macs && $(MACS2) callpeak -t $$WD/$(word 1, $^) -c $$WD/$(word 2, $^) -f BAM -g hs -n runx1 -q 0.01 --broad 2>&1 | $(LOG)
	WD=$$(pwd) && cd macs && $(MACS2) callpeak -t $$WD/$(word 1, $^) -c $$WD/$(word 2, $^) -f BAM -g hs -n runx1 -q 0.01 --bw 1000 --nomodel --shiftsize=100 --broad 2>&1 | $(LOG)

macs/er_peaks.bed: bwa/32239_TGACCA_C80K5ANXX_6_20150930B_20150930.bwa.sorted.filtered.bam bwa/32238_CGATGT_C80K5ANXX_6_20150930B_20150930.bwa.sorted.filtered.bam
	mkdir -p macs
	#WD=$$(pwd) && cd macs && $(MACS2) callpeak -t $$WD/$(word 1, $^) -c $$WD/$(word 2, $^) -f BAM -g hs -n er -q 0.01 --broad 2>&1 | $(LOG)
	WD=$$(pwd) && cd macs && $(MACS2) callpeak -t $$WD/$(word 1, $^) -c $$WD/$(word 2, $^) -f BAM -g hs -n er -q 0.01 --bw 1000 --nomodel --shiftsize=100 --broad 2>&1 | $(LOG)

macs/rhd_peaks.bed: bwa/32241_GCCAAT_C80K5ANXX_6_20150930B_20150930.bwa.sorted.filtered.bam bwa/32240_ACAGTG_C80K5ANXX_6_20150930B_20150930.bwa.sorted.filtered.bam
	mkdir -p macs
	#WD=$$(pwd) && cd macs && $(MACS2) callpeak -t $$WD/$(word 1, $^) -c $$WD/$(word 2, $^) -f BAM -g hs -n rhd -q 0.01 --broad 2>&1 | $(LOG)
	WD=$$(pwd) && cd macs && $(MACS2) callpeak -t $$WD/$(word 1, $^) -c $$WD/$(word 2, $^) -f BAM -g hs -n rhd -q 0.01 --bw 1000 --nomodel --shiftsize=100 --broad 2>&1 | $(LOG)

# --------------------------------------------------------------------------------
# peak annotation (Homer)
# --------------------------------------------------------------------------------
.PHONY: homer
homer: homer/runx1_peaks.annotated.with-expr.tsv homer/er_peaks.annotated.with-expr.tsv homer/rhd_peaks.annotated.with-expr.tsv \
       motifs/runx1_motifs.homer motifs/er_motifs.homer motifs/rhd_motifs.homer
	
homer/%_peaks.ucsc.bed: macs/%_peaks.bed
	mkdir -p homer
	cat $< | perl -ne 'if (/^(\d|X|Y)/) { print "chr$$_" } else { print $$_ }' > $@.part
	mv $@.part $@

homer/%_peaks.annotated.tsv: homer/%_peaks.ucsc.bed /mnt/projects/fiona/data/runx1.motif
	mkdir -p homer/$*
	$(HOMER) annotatePeaks.pl \
		$< \
		hg19 \
		-go homer/$*/ \
		-annStats homer/$*.annStats \
		-genomeOntology homer/$*/ \
		-cons \
		-CpG \
		-m /mnt/projects/fiona/data/runx1.motif \
		   /mnt/projects/fiona/data/ets1.motif \
		   /mnt/projects/fiona/data/ebf.motif \
		-mbed homer/$*_peaks.runx1-motif.bed \
		> $@.part
	mv $@.part $@
	cat homer/$*_peaks.runx1-motif.bed | grep -v "^track" | sort -k 1,1 -k2g,2g > homer/$*-peak-motifs.sorted.bed
	
homer/%_peaks.annotated.with-expr.tsv: homer/%_peaks.annotated.tsv macs/%_summits.bed anduril/execute/deseqAnnotated_oeERvsEmpty/table.csv anduril/execute/deseqAnnotated_oeRHDvsEmpty/table.csv /mnt/projects/fiona/scripts/annotate-peaks-with-expression.R
	Rscript /mnt/projects/fiona/scripts/annotate-peaks-with-expression.R --peak-file $(word 1, $^) --summit-file $(word 2, $^) --out-file $@.part
	mv $@.part $@
	
motifs/%_motifs.homer: homer/%_peaks.ucsc.bed
	mkdir -p motifs/$*
	$(HOMER) findMotifsGenome.pl $< hg19 motifs/$* -size 200 -mask -p 15 -fdr 100 > $@.part
	mv $@.part $@

# --------------------------------------------------------------------------------
# annotate peaks with expression data
# --------------------------------------------------------------------------------
	 