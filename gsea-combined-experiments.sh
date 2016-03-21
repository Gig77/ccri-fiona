#!/bin/bash

mkdir -p /mnt/projects/fiona/results/gsea_combined_experiments
cd /mnt/projects/fiona/results/gsea_combined_experiments

java -cp /data_synology/software/gsea-2.0.13/gsea2-2.0.13.jar -Xmx20g xtools.gsea.GseaPreranked \
	-rpt_label ETV6_RUNX1 \
	-rnk /mnt/projects/fiona/results/etv6runx1_combined_experiments.rnk \
	-gmx /mnt/projects/generic/data/msigdb5.0/msigdb.v5.0.symbols.gmt,/mnt/projects/generic/data/GeneSigDB/ALL_SIGSv4.nodup.gmt,/mnt/projects/generic/data/DSigDB/DSigDB_v1.0_All.nodup.gmt,/mnt/projects/generic/data/ccri/ccri_custom_gene_sets.gmt,/mnt/projects/generic/data/ccri/ccri_literature_curated_genesets_gsea.gmt,/mnt/projects/iamp/data/anduril/encode_tf_chipseq.ucsc.hg19.gmt,/mnt/projects/generic/data/pazar/pazar.gmt,/mnt/projects/generic/data/opossum3/jaspar_core.gmt \
	-out /mnt/projects/fiona/results/gsea_combined_experiments \
	-plot_top_x 3000 \
	-collapse false \
	-mode Max_probe  \
	-norm meandiv \
	-scoring_scheme weighted \
	-include_only_symbols true \
	-make_sets true \
	-rnd_seed 149 \
	-zip_report false \
	-gui false \
	-nperm 10 \
	-set_max 2000 \
	-set_min 5