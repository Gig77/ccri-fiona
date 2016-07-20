export SHELLOPTS:=errexit:pipefail
SHELL=/bin/bash  # required to make pipefail work
.SECONDARY:      # do not delete any intermediate files
LOG = perl -ne 'use POSIX qw(strftime); $$|=1; print strftime("%F %02H:%02M:%S ", localtime), $$ARGV[0], "$@: $$_";'

PROJECT_HOME=/mnt/projects/fiona
DOCKER=docker run -i --rm --net=host -e DOCKER_UID=$$(id -u) -e DOCKER_UNAME=$$(id -un) -e DOCKER_GID=$$(id -g) -e DOCKER_GNAME=$$(id -gn) -e DOCKER_HOME=$$HOME -v /home:/home -v /data_synology:/data_synology -v /data:/data -v /data2:/data2 -v /home/cf/fiona/results/current:$(PROJECT_HOME)/results -v /home/cf/fiona/data:$(PROJECT_HOME)/data -w $$(pwd)
MACS2=$(DOCKER) biowaste:5000/ccri/macs-2.0.9 macs2
HOMER=$(DOCKER) -v /home/cf/fiona/results/current/homer/wgEncodeAwgSegmentationCombinedEnhancers.ann.txt:/usr/local/homer/data/genomes/hg19/annotations/custom/Enhancer.ann.txt biowaste:5000/ccri/homer-4.7
BWA=/data_synology/software/bwa-0.7.12/bwa
SAMTOOLS=/data_synology/software/samtools-1.1/samtools
BEDTOOLS=/data_synology/software/bedtools-2.25.0/bin/bedtools
FASTQC=/data_synology/software/FastQC-0.11.2/fastqc
PICARD_MARKDUPLICATES=java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar /data_synology/software/picard-tools-1.114/MarkDuplicates.jar
TOPN=1000

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

SAMPLES_CHIPSEQ_2ND_BATCH=35115_ACAGTG_C81DHANXX_8_20160104B_20160104 \
                          35116_GCCAAT_C81DHANXX_8_20160104B_20160104 \
                          35117_CAGATC_C81DHANXX_8_20160104B_20160104 \
                          35118_CTTGTA_C81DHANXX_8_20160104B_20160104 \
                          35119_CGATGT_C8202ANXX_3_20160105B_20160105 \
                          35120_TGACCA_C8202ANXX_3_20160105B_20160105 \
                          35121_ACTTGA_C8202ANXX_3_20160105B_20160105 \
                          35122_GATCAG_C8202ANXX_3_20160105B_20160105 \
                          35123_TAGCTT_C8202ANXX_3_20160105B_20160105 \
                          35124_GGCTAC_C8202ANXX_3_20160105B_20160105

all: fastqc bwa flagstat homer motifs

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
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32240_ACAGTG_C80K5ANXX_6_20150930B_20150930.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32241_GCCAAT_C80K5ANXX_6_20150930B_20150930.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32242_CAGATC_C80K5ANXX_6_20150930B_20150930.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/32243_CTTGTA_C80K5ANXX_6_20150930B_20150930.bam

	wget -c --no-check-certificate --auth-no-challenge --limit-rate=11500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35109_CGATGT_C8918ANXX_8_20151231B_20151231.bam
	wget -c --no-check-certificate --auth-no-challenge --limit-rate=11500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35110_TGACCA_C8918ANXX_8_20151231B_20151231.bam
	wget -c --no-check-certificate --auth-no-challenge --limit-rate=11500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35111_ACAGTG_C8918ANXX_8_20151231B_20151231.bam
	wget -c --no-check-certificate --auth-no-challenge --limit-rate=11500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35112_GCCAAT_C8918ANXX_8_20151231B_20151231.bam
	wget -c --no-check-certificate --auth-no-challenge --limit-rate=11500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35113_CAGATC_C8918ANXX_8_20151231B_20151231.bam
	wget -c --no-check-certificate --auth-no-challenge --limit-rate=11500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35114_CTTGTA_C8918ANXX_8_20151231B_20151231.bam

	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35115_ACAGTG_C81DHANXX_8_20160104B_20160104.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35116_GCCAAT_C81DHANXX_8_20160104B_20160104.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35117_CAGATC_C81DHANXX_8_20160104B_20160104.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35118_CTTGTA_C81DHANXX_8_20160104B_20160104.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35119_CGATGT_C8202ANXX_3_20160105B_20160105.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35120_TGACCA_C8202ANXX_3_20160105B_20160105.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35121_ACTTGA_C8202ANXX_3_20160105B_20160105.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35122_GATCAG_C8202ANXX_3_20160105B_20160105.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35123_TAGCTT_C8202ANXX_3_20160105B_20160105.bam
	#wget -c --no-check-certificate --auth-no-challenge --limit-rate=1500k --user 'Fiona.Kraler' --password '4oLYn0a6le' http://ngs.csf.ac.at/data/35124_GGCTAC_C8202ANXX_3_20160105B_20160105.bam

# --------------------------------------------------------------------------------
# FASTQC
# --------------------------------------------------------------------------------

.PHONY: fastqc
fastqc: $(foreach S, $(SAMPLES_CHIPSEQ) $(SAMPLES_CHIPSEQ_2ND_BATCH), fastqc/$S_fastqc.html)

fastqc/%_fastqc.html: $(PROJECT_HOME)/data/bam/%.bam
	mkdir -p fastqc/$*.part
	$(FASTQC) -o fastqc/$*.part -f bam $^
	mv fastqc/$*.part/* fastqc
	rmdir fastqc/$*.part

# --------------------------------------------------------------------------------
# samtools flagstat
# --------------------------------------------------------------------------------

.PHONY: flagstat
flagstat: $(foreach S, $(SAMPLES_CHIPSEQ) $(SAMPLES_CHIPSEQ_2ND_BATCH), qc/$S.samtools.flagstat)

qc/%.samtools.flagstat: bwa/%.bwa.sorted.bam
	mkdir -p qc
	$(SAMTOOLS) flagstat $< 2>&1 1>$@.part | $(LOG)
	mv $@.part $@

# --------------------------------------------------------------------------------
# alignment, sorting, duplication marking, filtering
# --------------------------------------------------------------------------------

.PHONY: bwa
bwa: $(foreach S, $(SAMPLES_CHIPSEQ) $(SAMPLES_CHIPSEQ_2ND_BATCH), bwa/$S.bwa.sorted.bam.bai)

bwa/%.bwa.sorted.bam: $(PROJECT_HOME)/data/bam/%.bam
	mkdir -p bwa
	$(BWA) aln -b -t 10 /mnt/projects/generic/data/broad/human_g1k_v37.fasta $< \
		| $(BWA) samse /mnt/projects/generic/data/broad/human_g1k_v37.fasta - $< \
		| $(SAMTOOLS) view -Shb -@ 2 - \
		| $(SAMTOOLS) sort -T bwa/$* -o $@.part -O bam -@ 2 -m 2G \
		| 2>&1 | $(LOG)
	mv $@.part $@

bwa/%.bwa.sorted.dupmarked.bam: bwa/%.bwa.sorted.bam
	mkdir -p picard
	$(PICARD_MARKDUPLICATES) \
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
	/data_synology/software/samtools-0.1.19/samtools index $^ $@.part 2>&1 | $(LOG)

# --------------------------------------------------------------------------------
# bedtools input coverage
# --------------------------------------------------------------------------------

.PHONY: bedtools
bedtools: bedtools/ChIP24_AT2_ER.zeroInputCoverage.bed \
		  bedtools/ChIP24_REH_ER.zeroInputCoverage.bed \
		  bedtools/ChIP22_NALM6_RUNX1.zeroInputCoverage.bed \
		  bedtools/ChIP23_NALM6_ER.zeroInputCoverage.bed \
		  bedtools/ChIP23_NALM6_RHD.zeroInputCoverage.bed

hg19.5kb.windows.bed: /mnt/projects/generic/data/broad/human_g1k_v37.chromsizes.tsv $(BEDTOOLS)
	$(BEDTOOLS) makewindows -g $< -w 5000 > $@.part
	mv $@.part $@

bedtools/ChIP24_AT2_ER.zeroInputCoverage.bed: bwa/35115_ACAGTG_C81DHANXX_8_20160104B_20160104.bwa.sorted.bam hg19.5kb.windows.bed /mnt/projects/generic/data/broad/human_g1k_v37.chromsizes.tsv $(BEDTOOLS)
	$(SAMTOOLS) view -bhq 10 $< | $(BEDTOOLS) intersect -a $(word 2, $^) -b stdin -g $(word 3, $^) -v -sorted -nobuf > $@.part
	mv $@.part $@

bedtools/ChIP24_REH_ER.zeroInputCoverage.bed: bwa/35117_CAGATC_C81DHANXX_8_20160104B_20160104.bwa.sorted.bam hg19.5kb.windows.bed /mnt/projects/generic/data/broad/human_g1k_v37.chromsizes.tsv $(BEDTOOLS)
	$(SAMTOOLS) view -bhq 10 $< | $(BEDTOOLS) intersect -a $(word 2, $^) -b stdin -g $(word 3, $^) -v -sorted -nobuf > $@.part
	mv $@.part $@

bedtools/ChIP22_NALM6_RUNX1.zeroInputCoverage.bed: bwa/35119_CGATGT_C8202ANXX_3_20160105B_20160105.bwa.sorted.bam hg19.5kb.windows.bed /mnt/projects/generic/data/broad/human_g1k_v37.chromsizes.tsv $(BEDTOOLS)
	$(SAMTOOLS) view -bhq 10 $< | $(BEDTOOLS) intersect -a $(word 2, $^) -b stdin -g $(word 3, $^) -v -sorted -nobuf > $@.part
	mv $@.part $@

bedtools/ChIP23_NALM6_ER.zeroInputCoverage.bed: bwa/35121_ACTTGA_C8202ANXX_3_20160105B_20160105.bwa.sorted.bam hg19.5kb.windows.bed /mnt/projects/generic/data/broad/human_g1k_v37.chromsizes.tsv $(BEDTOOLS)
	$(SAMTOOLS) view -bhq 10 $< | $(BEDTOOLS) intersect -a $(word 2, $^) -b stdin -g $(word 3, $^) -v -sorted -nobuf > $@.part
	mv $@.part $@

bedtools/ChIP23_NALM6_RHD.zeroInputCoverage.bed: bwa/35123_TAGCTT_C8202ANXX_3_20160105B_20160105.bwa.sorted.bam hg19.5kb.windows.bed /mnt/projects/generic/data/broad/human_g1k_v37.chromsizes.tsv $(BEDTOOLS)
	$(SAMTOOLS) view -bhq 10 $< | $(BEDTOOLS) intersect -a $(word 2, $^) -b stdin -g $(word 3, $^) -v -sorted -nobuf > $@.part
	mv $@.part $@

# --------------------------------------------------------------------------------
# peak calling (MACS)
# --------------------------------------------------------------------------------

.PHONY: macs
macs: macs/runx1_peaks.bed \
	  macs/er_peaks.bed \
	  macs/rhd_peaks.bed \
	  macs/ChIP24_AT2_ER_peaks.bed macs/ChIP24_AT2_ER_model.pdf \
	  macs/ChIP24_REH_ER_peaks.bed macs/ChIP24_REH_ER_model.pdf \
	  macs/ChIP22_NALM6_RUNX1_peaks.bed macs/ChIP22_NALM6_RUNX1_model.pdf \
	  macs/ChIP23_NALM6_ER_peaks.bed macs/ChIP23_NALM6_ER_model.pdf \
	  macs/ChIP23_NALM6_RHD_peaks.bed

macs/runx1_peaks.bed: bwa/32243_CTTGTA_C80K5ANXX_6_20150930B_20150930.bwa.sorted.filtered.bam bwa/32242_CAGATC_C80K5ANXX_6_20150930B_20150930.bwa.sorted.filtered.bam
	mkdir -p macs
	WD=$$(pwd) && cd macs && $(MACS2) callpeak -t $$WD/$(word 1, $^) -c $$WD/$(word 2, $^) -f BAM -g hs -n runx1 -q 0.01 --bw 1000 --nomodel --shiftsize=100 --broad 2>&1 | $(LOG)

macs/er_peaks.bed: bwa/32239_TGACCA_C80K5ANXX_6_20150930B_20150930.bwa.sorted.filtered.bam bwa/32238_CGATGT_C80K5ANXX_6_20150930B_20150930.bwa.sorted.filtered.bam
	mkdir -p macs
	WD=$$(pwd) && cd macs && $(MACS2) callpeak -t $$WD/$(word 1, $^) -c $$WD/$(word 2, $^) -f BAM -g hs -n er -q 0.01 --bw 1000 --nomodel --shiftsize=100 --broad 2>&1 | $(LOG)

macs/rhd_peaks.bed: bwa/32241_GCCAAT_C80K5ANXX_6_20150930B_20150930.bwa.sorted.filtered.bam bwa/32240_ACAGTG_C80K5ANXX_6_20150930B_20150930.bwa.sorted.filtered.bam
	mkdir -p macs
	WD=$$(pwd) && cd macs && $(MACS2) callpeak -t $$WD/$(word 1, $^) -c $$WD/$(word 2, $^) -f BAM -g hs -n rhd -q 0.01 --bw 1000 --nomodel --shiftsize=100 --broad 2>&1 | $(LOG)

macs/ChIP24_AT2_ER_peaks.bed: bwa/35116_GCCAAT_C81DHANXX_8_20160104B_20160104.bwa.sorted.filtered.bam bwa/35115_ACAGTG_C81DHANXX_8_20160104B_20160104.bwa.sorted.filtered.bam
	mkdir -p macs
	WD=$$(pwd) && cd macs && $(MACS2) callpeak -t $$WD/$(word 1, $^) -c $$WD/$(word 2, $^) -f BAM -g hs -n ChIP24_AT2_ER -q 0.01 --broad 2>&1 | $(LOG)

macs/ChIP24_REH_ER_peaks.bed: bwa/35118_CTTGTA_C81DHANXX_8_20160104B_20160104.bwa.sorted.filtered.bam bwa/35117_CAGATC_C81DHANXX_8_20160104B_20160104.bwa.sorted.filtered.bam
	mkdir -p macs
	WD=$$(pwd) && cd macs && $(MACS2) callpeak -t $$WD/$(word 1, $^) -c $$WD/$(word 2, $^) -f BAM -g hs -n ChIP24_REH_ER -q 0.01 --broad 2>&1 | $(LOG)

macs/ChIP22_NALM6_RUNX1_peaks.bed: bwa/35120_TGACCA_C8202ANXX_3_20160105B_20160105.bwa.sorted.filtered.bam bwa/35119_CGATGT_C8202ANXX_3_20160105B_20160105.bwa.sorted.filtered.bam
	mkdir -p macs
	WD=$$(pwd) && cd macs && $(MACS2) callpeak -t $$WD/$(word 1, $^) -c $$WD/$(word 2, $^) -f BAM -g hs -n ChIP22_NALM6_RUNX1 -q 0.01 --broad 2>&1 | $(LOG)

macs/ChIP23_NALM6_ER_peaks.bed: bwa/35122_GATCAG_C8202ANXX_3_20160105B_20160105.bwa.sorted.filtered.bam bwa/35121_ACTTGA_C8202ANXX_3_20160105B_20160105.bwa.sorted.filtered.bam
	mkdir -p macs
	WD=$$(pwd) && cd macs && $(MACS2) callpeak -t $$WD/$(word 1, $^) -c $$WD/$(word 2, $^) -f BAM -g hs -n ChIP23_NALM6_ER -q 0.01 --broad 2>&1 | $(LOG)

macs/ChIP23_NALM6_RHD_peaks.bed: bwa/35124_GGCTAC_C8202ANXX_3_20160105B_20160105.bwa.sorted.filtered.bam bwa/35123_TAGCTT_C8202ANXX_3_20160105B_20160105.bwa.sorted.filtered.bam
	mkdir -p macs
	WD=$$(pwd) && cd macs && $(MACS2) callpeak -t $$WD/$(word 1, $^) -c $$WD/$(word 2, $^) -f BAM -g hs -n ChIP23_NALM6_RHD -q 0.01 --bw 1000 --nomodel --shiftsize=100 --broad 2>&1 | $(LOG)

macs/%_top$(TOPN)_peaks.bed: macs/%_peaks.bed
	head -$(TOPN) <(sort -k5,5nr $<) > $@.part
	cut -f 4 $@.part | sed 's/peak/summit/' > $@.summitids
	grep -wf $@.summitids $(word 2, $^) > macs/$*_top$(TOPN)_summits.bed.part
	mv $@.part $@
	mv macs/$*_top$(TOPN)_summits.bed.part macs/$*_top$(TOPN)_summits.bed
	rm $@.summitids

macs/%_runx1Motif_peaks.bed: homer/%_peaks.annotated.with-expr.tsv
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$22 != "" {print $$2,$$3,$$4,$$1,$$5}' $< | sed 's/^chr//' > $@.part
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$22 != "" {print $$2,$$27,$$27+1,$$1,$$5}' $< | sed 's/^chr//' | sed 's/MACS_peak/MACS_summit/' > macs/$*_runx1Motif_summits.bed.part
	mv macs/$*_runx1Motif_summits.bed.part macs/$*_runx1Motif_summits.bed
	mv $@.part $@

macs/%_norunx1Motif_peaks.bed: homer/%_peaks.annotated.with-expr.tsv
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$22 == "" {print $$2,$$3,$$4,$$1,$$5}' $< | sed 's/^chr//' > $@.part
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$22 == "" {print $$2,$$27,$$27+1,$$1,$$5}' $< | sed 's/^chr//' | sed 's/MACS_peak/MACS_summit/' > macs/$*_norunx1Motif_summits.bed.part
	mv macs/$*_norunx1Motif_summits.bed.part macs/$*_norunx1Motif_summits.bed
	mv $@.part $@

macs/%_shared_peaks.bed: homer/%_peaks.annotated.with-expr.tsv
	awk -F "\t" '{OFS="\t"} NR == 1 { next } ($$34 == "shared_better" || $$34 == "shared_worse") {print $$2,$$3,$$4,$$1,$$5}' $< | sed 's/^chr//' > $@.part
	awk -F "\t" '{OFS="\t"} NR == 1 { next } ($$34 == "shared_better" || $$34 == "shared_worse") {print $$2,$$27,$$27+1,$$1,$$5}' $< | sed 's/^chr//' | sed 's/MACS_peak/MACS_summit/' > macs/$*_shared_summits.bed.part
	mv macs/$*_shared_summits.bed.part macs/$*_shared_summits.bed
	mv $@.part $@

macs/%_unique_peaks.bed: homer/%_peaks.annotated.with-expr.tsv
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$34 == "unique" {print $$2,$$3,$$4,$$1,$$5}' $< | sed 's/^chr//' > $@.part
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$34 == "unique" {print $$2,$$27,$$27+1,$$1,$$5}' $< | sed 's/^chr//' | sed 's/MACS_peak/MACS_summit/' > macs/$*_unique_summits.bed.part
	mv macs/$*_unique_summits.bed.part macs/$*_unique_summits.bed
	mv $@.part $@

macs/%_better_peaks.bed: homer/%_peaks.annotated.with-expr.tsv
	awk -F "\t" '{OFS="\t"} NR == 1 { next } ($$34 == "unique" || $$34 == "shared_better") {print $$2,$$3,$$4,$$1,$$5}' $< | sed 's/^chr//' > $@.part
	awk -F "\t" '{OFS="\t"} NR == 1 { next } ($$34 == "unique" || $$34 == "shared_better") {print $$2,$$27,$$27+1,$$1,$$5}' $< | sed 's/^chr//' | sed 's/MACS_peak/MACS_summit/' > macs/$*_better_summits.bed.part
	mv macs/$*_better_summits.bed.part macs/$*_better_summits.bed
	mv $@.part $@

macs/%_worse_peaks.bed: homer/%_peaks.annotated.with-expr.tsv
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$34 == "shared_worse" {print $$2,$$3,$$4,$$1,$$5}' $< | sed 's/^chr//' > $@.part
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$34 == "shared_worse" {print $$2,$$27,$$27+1,$$1,$$5}' $< | sed 's/^chr//' | sed 's/MACS_peak/MACS_summit/' > macs/$*_worse_summits.bed.part
	mv macs/$*_worse_summits.bed.part macs/$*_worse_summits.bed
	mv $@.part $@

macs/%_shuffled_peaks.bed: macs/%_peaks.bed /mnt/projects/generic/data/broad/human_g1k_v37.chromsizes.tsv bedtools/%.zeroInputCoverage.bed
	$(BEDTOOLS) shuffle -chrom -l 4711 -i $< -g $(word 2, $^) -excl $(word 3, $^) > $@.part
	awk -F "\t" '{OFS="\t"} {summit=sprintf("%d", $$2+($$3-$$2)/2) ; print $$1,summit,summit+1,$$4,$$5}' $@.part > macs/$*_shuffled_summits.bed.part
	mv macs/$*_shuffled_summits.bed.part macs/$*_shuffled_summits.bed
	mv $@.part $@

macs/diffbind_ER_vs_RUNX_peaks.bed: /mnt/projects/fiona/scripts/diffBind.R
	Rscript $<
	mv $@.part $@

macs/%_model.pdf: macs/%_model.r
	cd macs && Rscript ../$<

# --------------------------------------------------------------------------------
# HOMER bedgraph files containing normalized read coverage
# --------------------------------------------------------------------------------
.PHONY: tags
tags: homer-tagfiles/ChIP24_AT2_ER_input.ucsc.bedGraph.gz \
	  homer-tagfiles/ChIP24_AT2_ER_chip.ucsc.bedGraph.gz \
	  homer-tagfiles/ChIP24_REH_ER_input.ucsc.bedGraph.gz \
	  homer-tagfiles/ChIP24_REH_ER_chip.ucsc.bedGraph.gz \
	  homer-tagfiles/ChIP22_NALM6_RUNX1_input.ucsc.bedGraph.gz \
	  homer-tagfiles/ChIP22_NALM6_RUNX1_chip.ucsc.bedGraph.gz \
	  homer-tagfiles/ChIP23_NALM6_ER_input.ucsc.bedGraph.gz \
	  homer-tagfiles/ChIP23_NALM6_ER_chip.ucsc.bedGraph.gz \
	  homer-tagfiles/ChIP23_NALM6_RHD_input.ucsc.bedGraph.gz \
	  homer-tagfiles/ChIP23_NALM6_RHD_chip.ucsc.bedGraph.gz

homer-tagfiles/%.ucsc.bedGraph.gz: homer-tagfiles/%/tagInfo.txt
	$(HOMER) makeUCSCfile homer-tagfiles/$* -norm 3e7 -o auto
	zcat homer-tagfiles/$*/$*.ucsc.bedGraph.gz | perl -ne 's/^(\d+|MT|X|Y)/chr$$1/; print $$_;' | gzip -c > $@.part
	rm -f $@ homer-tagfiles/$*/$*.ucsc.bedGraph.gz
	mv $@.part $@

homer-tagfiles/ChIP24_AT2_ER_input/tagInfo.txt: bwa/35115_ACAGTG_C81DHANXX_8_20160104B_20160104.bwa.sorted.filtered.bam
	$(HOMER) makeTagDirectory $(dir $@) $<
homer-tagfiles/ChIP24_AT2_ER_chip/tagInfo.txt: bwa/35116_GCCAAT_C81DHANXX_8_20160104B_20160104.bwa.sorted.filtered.bam
	$(HOMER) makeTagDirectory $(dir $@) $<
homer-tagfiles/ChIP24_REH_ER_input/tagInfo.txt: bwa/35117_CAGATC_C81DHANXX_8_20160104B_20160104.bwa.sorted.filtered.bam
	$(HOMER) makeTagDirectory $(dir $@) $<
homer-tagfiles/ChIP24_REH_ER_chip/tagInfo.txt: bwa/35118_CTTGTA_C81DHANXX_8_20160104B_20160104.bwa.sorted.filtered.bam
	$(HOMER) makeTagDirectory $(dir $@) $<
homer-tagfiles/ChIP22_NALM6_RUNX1_input/tagInfo.txt: bwa/35119_CGATGT_C8202ANXX_3_20160105B_20160105.bwa.sorted.filtered.bam
	$(HOMER) makeTagDirectory $(dir $@) $<
homer-tagfiles/ChIP22_NALM6_RUNX1_chip/tagInfo.txt: bwa/35120_TGACCA_C8202ANXX_3_20160105B_20160105.bwa.sorted.filtered.bam
	$(HOMER) makeTagDirectory $(dir $@) $<
homer-tagfiles/ChIP23_NALM6_ER_input/tagInfo.txt: bwa/35121_ACTTGA_C8202ANXX_3_20160105B_20160105.bwa.sorted.filtered.bam
	$(HOMER) makeTagDirectory $(dir $@) $<
homer-tagfiles/ChIP23_NALM6_ER_chip/tagInfo.txt: bwa/35122_GATCAG_C8202ANXX_3_20160105B_20160105.bwa.sorted.filtered.bam
	$(HOMER) makeTagDirectory $(dir $@) $<
homer-tagfiles/ChIP23_NALM6_RHD_input/tagInfo.txt: bwa/35123_TAGCTT_C8202ANXX_3_20160105B_20160105.bwa.sorted.filtered.bam
	$(HOMER) makeTagDirectory $(dir $@) $<
homer-tagfiles/ChIP23_NALM6_RHD_chip/tagInfo.txt: bwa/35124_GGCTAC_C8202ANXX_3_20160105B_20160105.bwa.sorted.filtered.bam
	$(HOMER) makeTagDirectory $(dir $@) $<

# --------------------------------------------------------------------------------
# peak annotation (Homer)
# --------------------------------------------------------------------------------
.PHONY: homer
homer: homer/runx1_peaks.annotated.with-expr.tsv \
       homer/er_peaks.annotated.with-expr.tsv \
       homer/rhd_peaks.annotated.with-expr.tsv \
       homer/ChIP24_AT2_ER_peaks.annotated.with-expr.tsv \
	   homer/ChIP24_AT2_ER_top$(TOPN)_peaks.annotated.with-expr.tsv \
       homer/ChIP24_AT2_ER_runx1Motif_peaks.annotated.tsv \
       homer/ChIP24_AT2_ER_norunx1Motif_peaks.annotated.tsv \
       homer/ChIP24_AT2_ER_shared_peaks.annotated.tsv \
       homer/ChIP24_AT2_ER_unique_peaks.annotated.tsv \
       homer/ChIP24_AT2_ER_better_peaks.annotated.with-expr.tsv \
       homer/ChIP24_AT2_ER_worse_peaks.annotated.with-expr.tsv \
       homer/ChIP24_AT2_ER_shuffled_peaks.annotated.with-expr.tsv \
       homer/ChIP24_REH_ER_peaks.annotated.with-expr.tsv \
	   homer/ChIP24_REH_ER_top$(TOPN)_peaks.annotated.with-expr.tsv \
       homer/ChIP24_REH_ER_runx1Motif_peaks.annotated.tsv \
       homer/ChIP24_REH_ER_norunx1Motif_peaks.annotated.tsv \
       homer/ChIP24_REH_ER_shared_peaks.annotated.tsv \
       homer/ChIP24_REH_ER_unique_peaks.annotated.tsv \
       homer/ChIP24_REH_ER_better_peaks.annotated.with-expr.tsv \
       homer/ChIP24_REH_ER_worse_peaks.annotated.with-expr.tsv \
       homer/ChIP24_REH_ER_shuffled_peaks.annotated.with-expr.tsv \
       homer/ChIP22_NALM6_RUNX1_peaks.annotated.with-expr.tsv \
	   homer/ChIP22_NALM6_RUNX1_top$(TOPN)_peaks.annotated.with-expr.tsv \
       homer/ChIP22_NALM6_RUNX1_runx1Motif_peaks.annotated.tsv \
       homer/ChIP22_NALM6_RUNX1_norunx1Motif_peaks.annotated.tsv \
       homer/ChIP22_NALM6_RUNX1_shuffled_peaks.annotated.with-expr.tsv \
       homer/ChIP23_NALM6_ER_peaks.annotated.with-expr.tsv \
	   homer/ChIP23_NALM6_ER_top$(TOPN)_peaks.annotated.with-expr.tsv \
       homer/ChIP23_NALM6_ER_runx1Motif_peaks.annotated.tsv \
       homer/ChIP23_NALM6_ER_norunx1Motif_peaks.annotated.tsv \
       homer/ChIP23_NALM6_ER_shared_peaks.annotated.tsv \
       homer/ChIP23_NALM6_ER_unique_peaks.annotated.tsv \
       homer/ChIP23_NALM6_ER_better_peaks.annotated.with-expr.tsv \
       homer/ChIP23_NALM6_ER_worse_peaks.annotated.with-expr.tsv \
       homer/ChIP23_NALM6_ER_shuffled_peaks.annotated.with-expr.tsv \
       homer/ChIP23_NALM6_RHD_peaks.annotated.with-expr.tsv \
       homer/diffbind_ER_vs_RUNX_peaks.annotated.with-expr.tsv

homer/wgEncodeAwgSegmentationCombinedEnhancers.ann.txt: /mnt/projects/fiona/data/enhancers/wgEncodeAwgSegmentationCombinedHelas3.bed \
														/mnt/projects/fiona/data/enhancers/wgEncodeAwgSegmentationCombinedHuvec.bed \
														/mnt/projects/fiona/data/enhancers/wgEncodeAwgSegmentationCombinedGm12878.bed \
														/mnt/projects/fiona/data/enhancers/wgEncodeAwgSegmentationCombinedH1hesc.bed \
														/mnt/projects/fiona/data/enhancers/wgEncodeAwgSegmentationCombinedHepg2.bed \
														/mnt/projects/fiona/data/enhancers/wgEncodeAwgSegmentationCombinedK562.bed
	cat $^ | grep -P "\tE\t" | sort -k 1,1 -k2n,2n | $(BEDTOOLS) merge -nobuf -i stdin \
		| perl -lane '$$id++; print "Enhancer-$$id\t$$F[0]\t$$F[1]\t$$F[2]\t0\tEnhancer"' > $@.part
	mv $@.part $@

motifs/all_motifs.motif: /mnt/projects/fiona/data/runx1.motif \
												 /mnt/projects/fiona/data/ets1.motif \
												 /mnt/projects/fiona/data/ebf.motif \
												 /mnt/projects/fiona/data/gata3.motif \
												 /mnt/projects/fiona/data/ets_runx.motif
	mkdir -p motifs
	cat $^ > $@.part
	mv $@.part $@

homer/%_peaks.ucsc.bed: macs/%_peaks.bed
	mkdir -p homer
	cat $< | perl -ne 'if (/^(\d|X|Y)/) { print "chr$$_" } else { print $$_ }' > $@.part
	mv $@.part $@

Tijssen_2011_TableS1_peak_coordinates.hg19.txt: /mnt/projects/fiona/scripts/liftover_tijssen.R /mnt/projects/fiona/data/Tijssen_2011_TableS1_peak_coordinates.txt
	Rscript /mnt/projects/fiona/scripts/liftover_tijssen.R

homer/%_peaks.annotated.tsv: homer/%_peaks.ucsc.bed motifs/all_motifs.motif homer/wgEncodeAwgSegmentationCombinedEnhancers.ann.txt
	mkdir -p homer/$*
	$(HOMER) annotatePeaks.pl \
		$< \
		hg19 \
		-annStats homer/$*.annStats \
		-gsize 3000000000 \
		-genomeOntology homer/$*/ \
		-go homer/$*/ \
		-cons \
		-CpG \
		-m motifs/all_motifs.motif \
		-mbed homer/$*_peaks.runx1-motif.bed \
		> $@.part
	mv $@.part $@
	cat homer/$*_peaks.runx1-motif.bed | grep -v "^track" | sort -k 1,1 -k2g,2g > homer/$*-peak-motifs.sorted.bed

homer/%_peaks.annotated.with-expr.tsv: homer/%_peaks.annotated.tsv anduril/execute/deseqAnnotated_oeERvsEmptyB1/table.csv anduril/execute/deseqAnnotated_oeRHDvsEmptyB1/table.csv anduril/execute/deseqAnnotated_oeERvsEmptyB2/table.csv anduril/execute/deseqAnnotated_oeRHDvsEmptyB2/table.csv homer/ChIP22_NALM6_RUNX1_peaks.annotated.tsv Tijssen_2011_TableS1_peak_coordinates.hg19.txt /mnt/projects/fiona/scripts/annotate-peaks.R
	Rscript /mnt/projects/fiona/scripts/annotate-peaks.R --homer-peak-file $(word 1, $^) --macs-peak-file macs/$*_peaks.bed --summit-file macs/$*_summits.bed --out-file $@.part
	mv $@.part $@

# --------------------------------------------------------------------------------
# motif discovery (Homer)
# --------------------------------------------------------------------------------
.PHONY: motifs
motifs: motifs/runx1_motifs.homer \
       motifs/er_motifs.homer \
       motifs/rhd_motifs.homer \
       motifs/ChIP24_AT2_ER_motifs.homer \
 	   motifs/ChIP24_AT2_ER_top$(TOPN)_motifs.homer \
       motifs/ChIP24_AT2_ER_runx1Motif_motifs.homer \
       motifs/ChIP24_AT2_ER_norunx1Motif_motifs.homer \
       motifs/ChIP24_AT2_ER_shared_motifs.homer \
       motifs/ChIP24_AT2_ER_unique_motifs.homer \
       motifs/ChIP24_AT2_ER_better_motifs.homer \
       motifs/ChIP24_AT2_ER_worse_motifs.homer \
       motifs/ChIP24_AT2_ER_shuffled_motifs.homer \
       motifs/ChIP24_AT2_ER.motif_hits.tsv \
       motifs/ChIP24_AT2_ER_top$(TOPN).motif_hits.tsv \
       motifs/ChIP24_REH_ER_motifs.homer \
	   motifs/ChIP24_REH_ER_top$(TOPN)_motifs.homer \
       motifs/ChIP24_REH_ER_runx1Motif_motifs.homer \
       motifs/ChIP24_REH_ER_norunx1Motif_motifs.homer \
       motifs/ChIP24_REH_ER_shared_motifs.homer \
       motifs/ChIP24_REH_ER_unique_motifs.homer \
       motifs/ChIP24_REH_ER_better_motifs.homer \
       motifs/ChIP24_REH_ER_worse_motifs.homer \
       motifs/ChIP24_REH_ER_shuffled_motifs.homer \
       motifs/ChIP24_REH_ER.motif_hits.tsv \
       motifs/ChIP24_REH_ER_top$(TOPN).motif_hits.tsv \
       motifs/ChIP23_NALM6_ER_motifs.homer \
	   motifs/ChIP23_NALM6_ER_top$(TOPN)_motifs.homer \
       motifs/ChIP23_NALM6_ER_runx1Motif_motifs.homer \
       motifs/ChIP23_NALM6_ER_norunx1Motif_motifs.homer \
       motifs/ChIP23_NALM6_ER_shared_motifs.homer \
       motifs/ChIP23_NALM6_ER_unique_motifs.homer \
       motifs/ChIP23_NALM6_ER_better_motifs.homer \
       motifs/ChIP23_NALM6_ER_worse_motifs.homer \
       motifs/ChIP23_NALM6_ER_shuffled_motifs.homer \
       motifs/ChIP23_NALM6_ER.motif_hits.tsv \
       motifs/ChIP23_NALM6_ER_top$(TOPN).motif_hits.tsv \
       motifs/ChIP22_NALM6_RUNX1_motifs.homer \
	   motifs/ChIP22_NALM6_RUNX1_top$(TOPN)_motifs.homer \
       motifs/ChIP22_NALM6_RUNX1_runx1Motif_motifs.homer \
       motifs/ChIP22_NALM6_RUNX1_norunx1Motif_motifs.homer \
       motifs/ChIP22_NALM6_RUNX1_shuffled_motifs.homer \
       motifs/ChIP22_NALM6_RUNX1.motif_hits.tsv \
       motifs/ChIP22_NALM6_RUNX1_top$(TOPN).motif_hits.tsv \
       motifs/ChIP23_NALM6_RHD_motifs.homer

motifs/%_summit-regions.bed: macs/%_peaks.bed /mnt/projects/fiona/scripts/summit2region.R
	Rscript /mnt/projects/fiona/scripts/summit2region.R --summit-file macs/%_summits.bed --size 100 --out-file $@.part 2>&1 | $(LOG)
	mv $@.part $@

motifs/%_motifs.homer: motifs/%_summit-regions.bed
	mkdir -p motifs/$*
	$(HOMER) findMotifsGenome.pl $< hg19 motifs/$* -nomotif -size given -mask -p 15 -fdr 100 -dumpFasta > $@.part
	mv $@.part $@

motifs/%.motif_hits.tsv: homer/%_peaks.ucsc.bed motifs/all_motifs.motif
	$(HOMER) findMotifsGenome.pl $< hg19 motifs/$*.tmp -find motifs/all_motifs.motif > $@.part
	rm -rf motifs/$*.tmp
	mv $@.part $@

wgEncodeDacMapabilityConsensusExcludable.bed:
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz
	gunzip $@.gz
	sed 's/^chr//' $@ > wgEncodeDacMapabilityConsensusExcludable.nochr.bed

# --------------------------------------------------------------------------------
# GSEA
# --------------------------------------------------------------------------------

.PHONY: gsea
gsea: gsea/diffbind_ER_vs_RUNX.gsea \
	  gsea/diffbind_ER_vs_RUNX_proteincoding.gsea \
	  gsea/diffbind_ER_vs_RUNX_proteincoding_promoter.gsea \
	  gsea/diffbind_ER_vs_RUNX_proteincoding_runxmotif.gsea

gsea/%.rnk: homer/%_peaks.annotated.tsv
	mkdir -p gsea
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$16 != "" {print $$16,$$6}' $< | sort -u -k1,1 | sort -k2,2nr > $@.part
	mv $@.part $@

gsea/%_proteincoding.rnk: homer/%_peaks.annotated.tsv
	mkdir -p gsea
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$16 != "" && $$19 == "protein-coding" {print $$16,$$6}' $< | sort -u -k1,1 | sort -k2,2nr > $@.part
	mv $@.part $@

gsea/%_proteincoding_promoter.rnk: homer/%_peaks.annotated.tsv
	mkdir -p gsea
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$16 != "" && $$19 == "protein-coding" && $$10 >= -10000 && $$10 <= 5000 {print $$16,$$6}' $< | sort -u -k1,1 | sort -k2,2nr > $@.part
	mv $@.part $@

gsea/%_proteincoding_runxmotif.rnk: homer/%_peaks.annotated.tsv
	mkdir -p gsea
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$16 != "" && $$19 == "protein-coding" && $$22 != "" {print $$16,$$6}' $< | sort -u -k1,1 | sort -k2,2nr > $@.part
	mv $@.part $@

gsea/%.gsea: gsea/%.rnk
	java -cp /data_synology/software/gsea-2.0.13/gsea2-2.0.13.jar -Xmx3048m xtools.gsea.GseaPreranked \
		-rpt_label $* \
		-rnk $< \
		-gmx /mnt/projects/generic/data/msigdb5.0/msigdb.v5.0.symbols.gmt \
		-collapse false -mode Max_probe -norm meandiv -nperm 100 -scoring_scheme weighted -include_only_symbols true -make_sets true \
		-rnd_seed 149 \
		-plot_top_x 300 \
		-set_max 500 \
		-set_min 15 \
		-zip_report false \
		-gui false \
		-out gsea
