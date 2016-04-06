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

all: fastqc bwa flagstat homer

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
	/data_synology/software/FastQC-0.11.2/fastqc -o fastqc/$*.part -f bam $^
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
	/data_synology/software/samtools-0.1.19/samtools index $^ $@.part 2>&1 | $(LOG)
	
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

macs/%_runx1Motif_peaks.bed: homer/%_peaks.annotated.with-expr.tsv
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$22 != "" {print $$3,$$4,$$5,$$2,$$7}' $< | sed 's/^chr//' > $@.part
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$22 != "" {print $$3,$$26,$$26+1,$$2,$$7}' $< | sed 's/^chr//' | sed 's/MACS_peak/MACS_summit/' > macs/$*_runx1Motif_summits.bed.part
	mv macs/$*_runx1Motif_summits.bed.part macs/$*_runx1Motif_summits.bed 
	mv $@.part $@

macs/%_norunx1Motif_peaks.bed: homer/%_peaks.annotated.with-expr.tsv
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$22 == "" {print $$3,$$4,$$5,$$2,$$7}' $< | sed 's/^chr//' > $@.part
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$22 == "" {print $$3,$$26,$$26+1,$$2,$$7}' $< | sed 's/^chr//' | sed 's/MACS_peak/MACS_summit/' > macs/$*_norunx1Motif_summits.bed.part
	mv macs/$*_norunx1Motif_summits.bed.part macs/$*_norunx1Motif_summits.bed 
	mv $@.part $@

macs/%_constitutive_peaks.bed: homer/%_peaks.annotated.with-expr.tsv
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$32 == "constitutive" {print $$3,$$4,$$5,$$2,$$7}' $< | sed 's/^chr//' > $@.part
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$32 == "constitutive" {print $$3,$$26,$$26+1,$$2,$$7}' $< | sed 's/^chr//' | sed 's/MACS_peak/MACS_summit/' > macs/$*_constitutive_summits.bed.part
	mv macs/$*_constitutive_summits.bed.part macs/$*_constitutive_summits.bed 
	mv $@.part $@

macs/%_denovo_peaks.bed: homer/%_peaks.annotated.with-expr.tsv
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$32 == "de novo" {print $$3,$$4,$$5,$$2,$$7}' $< | sed 's/^chr//' > $@.part
	awk -F "\t" '{OFS="\t"} NR == 1 { next } $$32 == "de novo" {print $$3,$$26,$$26+1,$$2,$$7}' $< | sed 's/^chr//' | sed 's/MACS_peak/MACS_summit/' > macs/$*_denovo_summits.bed.part
	mv macs/$*_denovo_summits.bed.part macs/$*_denovo_summits.bed 
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
       homer/ChIP24_REH_ER_peaks.annotated.with-expr.tsv \
       homer/ChIP22_NALM6_RUNX1_peaks.annotated.with-expr.tsv \
       homer/ChIP23_NALM6_ER_peaks.annotated.with-expr.tsv \
       homer/ChIP23_NALM6_RHD_peaks.annotated.with-expr.tsv \
       homer/ChIP24_AT2_ER_runx1Motif_peaks.annotated.tsv \
       homer/ChIP24_REH_ER_runx1Motif_peaks.annotated.tsv \
       homer/ChIP22_NALM6_RUNX1_runx1Motif_peaks.annotated.tsv \
       homer/ChIP23_NALM6_ER_runx1Motif_peaks.annotated.tsv \
       homer/ChIP24_AT2_ER_norunx1Motif_peaks.annotated.tsv \
       homer/ChIP24_REH_ER_norunx1Motif_peaks.annotated.tsv \
       homer/ChIP22_NALM6_RUNX1_norunx1Motif_peaks.annotated.tsv \
       homer/ChIP23_NALM6_ER_norunx1Motif_peaks.annotated.tsv \
       homer/ChIP24_AT2_ER_constitutive_peaks.annotated.tsv \
       homer/ChIP24_REH_ER_constitutive_peaks.annotated.tsv \
       homer/ChIP23_NALM6_ER_constitutive_peaks.annotated.tsv \
       homer/ChIP24_AT2_ER_denovo_peaks.annotated.tsv \
       homer/ChIP24_REH_ER_denovo_peaks.annotated.tsv \
       homer/ChIP23_NALM6_ER_denovo_peaks.annotated.tsv \
       motifs/runx1_motifs.homer \
       motifs/er_motifs.homer \
       motifs/rhd_motifs.homer \
       motifs/ChIP24_AT2_ER_motifs.homer \
       motifs/ChIP24_REH_ER_motifs.homer \
       motifs/ChIP22_NALM6_RUNX1_motifs.homer \
       motifs/ChIP23_NALM6_ER_motifs.homer \
       motifs/ChIP23_NALM6_RHD_motifs.homer \
       motifs/ChIP24_AT2_ER_runx1Motif_motifs.homer \
       motifs/ChIP24_REH_ER_runx1Motif_motifs.homer \
       motifs/ChIP22_NALM6_RUNX1_runx1Motif_motifs.homer \
       motifs/ChIP23_NALM6_ER_runx1Motif_motifs.homer \
       motifs/ChIP24_AT2_ER_norunx1Motif_motifs.homer \
       motifs/ChIP24_REH_ER_norunx1Motif_motifs.homer \
       motifs/ChIP22_NALM6_RUNX1_norunx1Motif_motifs.homer \
       motifs/ChIP23_NALM6_ER_norunx1Motif_motifs.homer \
       motifs/ChIP24_AT2_ER_constitutive_motifs.homer \
       motifs/ChIP24_REH_ER_constitutive_motifs.homer \
       motifs/ChIP23_NALM6_ER_constitutive_motifs.homer \
       motifs/ChIP24_AT2_ER_denovo_motifs.homer \
       motifs/ChIP24_REH_ER_denovo_motifs.homer \
       motifs/ChIP23_NALM6_ER_denovo_motifs.homer \
       motifs/ChIP24_AT2_ER.runx1_hits.tsv      motifs/ChIP24_AT2_ER.ets_hits.tsv      motifs/ChIP24_AT2_ER.ebf_hits.tsv \
       motifs/ChIP24_REH_ER.runx1_hits.tsv      motifs/ChIP24_REH_ER.ets_hits.tsv      motifs/ChIP24_REH_ER.ebf_hits.tsv \
       motifs/ChIP23_NALM6_ER.runx1_hits.tsv    motifs/ChIP23_NALM6_ER.ets_hits.tsv    motifs/ChIP23_NALM6_ER.ebf_hits.tsv \
       motifs/ChIP22_NALM6_RUNX1.runx1_hits.tsv motifs/ChIP22_NALM6_RUNX1.ets_hits.tsv motifs/ChIP22_NALM6_RUNX1.ebf_hits.tsv
       
	
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
		   /mnt/projects/fiona/data/gata3.motif \
		-mbed homer/$*_peaks.runx1-motif.bed \
		> $@.part
	mv $@.part $@
	cat homer/$*_peaks.runx1-motif.bed | grep -v "^track" | sort -k 1,1 -k2g,2g > homer/$*-peak-motifs.sorted.bed
	
motifs/%_motifs.homer: homer/%_peaks.ucsc.bed
	mkdir -p motifs/$*
	$(HOMER) findMotifsGenome.pl $< hg19 motifs/$* -size 200 -mask -p 15 -fdr 100 -dumpFasta > $@.part
	mv $@.part $@

motifs/%.runx1_hits.tsv: homer/%_peaks.ucsc.bed /mnt/projects/fiona/data/runx1.motif
	$(HOMER) findMotifsGenome.pl $(word 1, $^) hg19 motifs/$*.tmp -find $(word 2, $^) > $@.part
	rm -rf motifs/$*.tmp
	mv $@.part $@

motifs/%.ebf_hits.tsv: homer/%_peaks.ucsc.bed /mnt/projects/fiona/data/ebf.motif
	$(HOMER) findMotifsGenome.pl $(word 1, $^) hg19 motifs/$*.tmp -find $(word 2, $^) > $@.part
	rm -rf motifs/$*.tmp
	mv $@.part $@

motifs/%.ets_hits.tsv: homer/%_peaks.ucsc.bed /mnt/projects/fiona/data/ets1.motif
	$(HOMER) findMotifsGenome.pl $(word 1, $^) hg19 motifs/$*.tmp -find $(word 2, $^) > $@.part
	rm -rf motifs/$*.tmp
	mv $@.part $@

homer/%_peaks.annotated.with-expr.tsv: homer/%_peaks.annotated.tsv anduril/execute/deseqAnnotated_oeERvsEmptyB1/table.csv anduril/execute/deseqAnnotated_oeRHDvsEmptyB1/table.csv anduril/execute/deseqAnnotated_oeERvsEmptyB2/table.csv anduril/execute/deseqAnnotated_oeRHDvsEmptyB2/table.csv homer/ChIP22_NALM6_RUNX1_peaks.annotated.tsv /mnt/projects/fiona/scripts/annotate-peaks.R
	Rscript /mnt/projects/fiona/scripts/annotate-peaks.R --peak-file $(word 1, $^) --summit-file macs/$*_summits.bed --out-file $@.part
	mv $@.part $@
	
	 