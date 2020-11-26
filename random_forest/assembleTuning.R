#!/usr/bin/env Rscript
# Michael VanInsberghe 2020-05-21
# Prepares tuning for ranger model

library(optparse)
library(dplyr)
library(tidyr)
library(data.table)
library(feather)
library(ggplot2)
library(gplots)
library(ggseqlogo)
library(RColorBrewer)
library(cowplot)
library(forcats)
library(ranger)
theme_set(theme_cowplot())

source('/hpc/hub_oudenaarden/mvanins/local/avopipelines/RPF/bin/scRibo/plotFunctions.R')

option_list = list(
	make_option(c("-r", "--reads"), type = "character", default = NULL, help = "feather file containing read information dataframe", metavar="character"),
	make_option(c("-a", "--annotations"), type = "character", default = NULL, help = "feather file containing annotation information", metavar="character"),
	make_option(c("-p", "--prefix"), type = "character", default = NULL, help = "prefix for output filenames", metavar = "character") )

#opt_parser = OptionParser(option_list = option_list)
#args = parse_args(opt_parser)
#
#if(is.null(args$reads) | is.null(args$annotations)){
#	print_help(opt_parser)
#	stop("Please supply reads and annotations", call.= FALSE)
#}
#
#if(is.null(args$prefix)){
#	figurePrefix = strsplit(tools::file_path_sans_ext(basename(args$reads)),"_")[[1]][1]
#}else{
#	figurePrefix = args$prefix
#}
#
#
#print(figurePrefix)
#
#if(grepl("feather",args$annotations)){
#		annot <- read_feather(args$annotations)
#}else{
#		annot <- fread(args$annotations)
#}

set.seed(1234)
figurePrefix <- "FucciModel_training"

annot <- read_feather('gencode.v33.chr_patch_hapl_scaff.annotation.annotations.feather')

canonical_pc <- annot %>%
  filter(set=="canonical", transcript_type == "protein_coding") %>%
  filter(chr != 'chrM') %>%
  filter(transcript_id != "ENST00000437139.7")

# barcodes

barcodes <- read.table('/hpc/hub_oudenaarden/mvanins/local/lib/barcodes/miRv6_barcodemap.tsv',
					   header = FALSE, col.names = c("CB","well_number","well"), sep = "\t", stringsAsFactors = FALSE)

# reads
# HEK
HEKsamples <- data.frame(names = c("HEK293T",
								"HEK293T-Starv1",
								"HEK293T-Starv2"),
					  reads = c("reads/RPFv4-HEK293T_Aligned.toTranscriptome.out.dedup.sorted.csv.gz",
								"reads/RPFv4-HEK293T-Starv1_Aligned.toTranscriptome.out.dedup.sorted.csv.gz",
								"reads/RPFv4-HEK293T-Starv2_Aligned.toTranscriptome.out.dedup.sorted.csv.gz"),
					  stringsAsFactors = FALSE)


HEKReads <- importReads(HEKsamples) %>%
		select(-id) %>%
		filter(transcript_id %in% canonical_pc$transcript_id) %>%
		separate(CB, c(NA,"barcode"), remove = FALSE, sep = "_") %>%
		filter(barcode %in% barcodes[barcodes$well %in% columnFeatures(1:23,LETTERS[1:16]),"CB"]) %>% # remove NTC wells
		select(-barcode)

# Fucci
FucciSamples <- data.frame(names = c("Fucci-Mit01",
								"Fucci-Mit02",
								"Fucci-Mit03",
								"Fucci-Mit04",
								"Fucci-Mit06",
								"Fucci-Mit07",
								"Fucci-Mit08",
								"Fucci-G011",
								"Fucci-G012"),
					  reads = c("reads/RPFv4-Fucci-Mit01_Aligned.toTranscriptome.out.dedup.sorted.csv.gz",
								"reads/RPFv4-Fucci-Mit02_Aligned.toTranscriptome.out.dedup.sorted.csv.gz",
								"reads/RPFv4-Fucci-Mit03_Aligned.toTranscriptome.out.dedup.sorted.csv.gz",
								"reads/RPFv4-Fucci-Mit04_Aligned.toTranscriptome.out.dedup.sorted.csv.gz",
								"reads/RPFv4-Fucci-Mit06_Aligned.toTranscriptome.out.dedup.sorted.csv.gz",
								"reads/RPFv4-Fucci-Mit07_Aligned.toTranscriptome.out.dedup.sorted.csv.gz",
								"reads/RPFv4-Fucci-Mit08_Aligned.toTranscriptome.out.dedup.sorted.csv.gz",
								"reads/RPFv4-Fucci-G011_Aligned.toTranscriptome.out.dedup.sorted.csv.gz",
								"reads/RPFv4-Fucci-G012_Aligned.toTranscriptome.out.dedup.sorted.csv.gz"),
					  stringsAsFactors = FALSE)

FucciReads <- importReads(FucciSamples) %>%
		select(-id) %>%
		filter(transcript_id %in% canonical_pc$transcript_id) %>%
		separate(CB, c(NA,"barcode"), remove = FALSE, sep = "_") %>%
		filter(barcode %in% barcodes[!(barcodes$well %in% columnFeatures(21:24,LETTERS[15:16])),"CB"]) %>% # remove NTC wells
		select(-barcode)

# filter cells 

reads <- bind_rows(HEKReads, FucciReads)

cds_fraction <- canonical_pc %>%
		mutate(l_utr5 = l_utr5 - 25,
			   l_cds = l_cds + 25,
			   l_utr3 = l_utr3 ) %>%
		mutate(l_utr5 = case_when( l_utr5 < 0 ~ 0,
								   l_utr5 >= 0 ~ l_utr5),
			   l_cds = case_when(l_cds < 0 ~ 0,
								 l_cds >= 0 ~ l_cds),
			   l_utr3 = case_when(l_utr3 < 0 ~ 0,
								  l_utr3 >= 0 ~ l_utr3)) %>%
		inner_join(reads, by = "transcript_id") %>%
		filter(length == read_length,
			   length >= 30, 
			   length <= 45,
			   read_strand == "+") %>%
		mutate( region = case_when( cut5 < (l_utr5) ~ "utr",
								((cut5  >= (l_utr5)) & (cut5 < ((l_utr5)+l_cds ) )) ~ "cds", 
								(cut5 >= (l_utr5 + l_cds)) ~ "utr" )) %>%
		group_by(CB,region) %>%
		summarize(n_region = n())%>%
	    ungroup() %>%
	    spread(region, n_region) %>%
		group_by(CB) %>%
		mutate(cds_frac = cds/(sum(cds,utr,na.rm=TRUE)),
			   cell_total = sum(cds,utr, na.rm=TRUE)) %>%
		ungroup() %>%
		filter(cds_frac > 0.75,
			   cell_total > 500)

reads <- reads %>%
		filter(CB %in% cds_fraction$CB) %>%
		filter(read_strand == "+",
			   length == read_length,
			   length >= 30,
			   length <= 45)

# Model Reads

stopReads <- canonical_pc %>%
		inner_join(reads, by="transcript_id") %>%
		mutate(cd5_stop = as.integer(cut5 - (cds_end)),
			   cd3_stop = as.integer(cut3 - (cds_end))) %>%
		filter(cd5_stop >= -18,
			   cd5_stop <= -14,
			   cd3_stop >= 0)

# make metaheatmap plots for each set
pdf(paste0(figurePrefix,"_reads_length_metaheatmap.pdf"), width = 8, height = 8)
for (s in c(FucciSamples$names, HEKsamples$names)){
		sreads <- reads %>%
				filter(sample == s,
					   read_strand == "+")
		plt_allreads_hm <- lengthMetaHeatmap(canonical_pc %>% 
										  inner_join(sreads, by="transcript_id"),
										  sites = c('cut5', 'cut3'), upstream = 60, downstream = 60)
		plt_stopreads_hm <- lengthMetaHeatmap(stopReads %>% filter(sample == s),
											  sites = c('cut5', 'cut3'), upstream = 60, downstream = 60)
		plt_reads_hm <- plot_grid(plt_allreads_hm+ggtitle(s), plt_stopreads_hm, ncol = 1)
		print(plt_reads_hm)
}
dev.off()


modelReads <- stopReads %>%
  mutate(cd5_stop = as.factor(cd5_stop)) %>%
  select(length,
         base5, base5m1, base5m2, base5m3, base5m4, base5m5, base5m6, base5m7, base5m8, base5p1, base5p2, base5p3, base5p4, base5p5, base5p6, base5p7, base5p8,
         base3, base3m1, base3m2, base3m3, base3m4, base3m5, base3m6, base3m7, base3m8, base3p1, base3p2, base3p3, base3p4, base3p5, base3p6, base3p7, base3p8,
         cd5_stop) %>%
  mutate_if(is.character, as.factor)

save(modelReads, file="modelReads.RData")
save.image("modelReads.image.RData")
