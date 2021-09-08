#!/usr/bin/env Rscript
# Michael VanInsberghe 2021-03-15
# Compile, merge, and filter Fucci cells

library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(scales)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(forcats)
library(RColorBrewer)
library(uwot)
library(ggrepel)
theme_set(theme_cowplot())

source('plotFunctions.R')

data_dir <- "../scRiboSeq/data"

# barcodes
barcodes <- read.table('/hpc/hub_oudenaarden/mvanins/local/lib/barcodes/miRv6_barcodemap.tsv',
                       sep="\t",
                       header = FALSE,
                       stringsAsFactors = FALSE,
                       row.names = NULL,
                       col.names = c("CB","number","well"))
rownames(barcodes) <- barcodes$CB

## human annot
annot <- fread(file.path(data_dir, 'hsa_annotations.csv.gz'), stringsAsFactors = FALSE, data.table = FALSE)
canonical_pc <- annot %>%
  filter(set=="canonical", transcript_type == "protein_coding") %>%
  #filter(chr != 'chrM') %>%
  filter(transcript_id != "ENST00000437139.7") %>%
  filter(!grepl("^PCDH", gene_name)) #%>%
  #filter(!grepl("^RPL", gene_name))

codons <- fread(file.path(data_dir, 'hsa_codons.csv.gz'), stringsAsFactors = FALSE, data.table = FALSE)
codons <- codons %>%
  filter( transcript_id %in% canonical_pc$transcript_id)



darnellSamples <- data.frame(names = c("darnell_Arg3h-r1",
                                       "darnell_Arg6h-r1",
                                       "darnell_Leu3h",
                                       "darnell_Rich3h",
                                       "darnell_Arg6h-r2",
                                       "darnell_Leu6h",
                                       "darnell_Rich6h"),
                             counts = c(file.path(data_dir, "bulk/darnell_2018/counts/SRR7073121"), 
                                        file.path(data_dir, "bulk/darnell_2018/counts/SRR7073122"), 
                                        file.path(data_dir, "bulk/darnell_2018/counts/SRR7073123"), 
                                        file.path(data_dir, "bulk/darnell_2018/counts/SRR7073124"), 
                                        file.path(data_dir, "bulk/darnell_2018/counts/SRR7073152"), 
                                        file.path(data_dir, "bulk/darnell_2018/counts/SRR7073153"), 
                                        file.path(data_dir, "bulk/darnell_2018/counts/SRR7073154")),
                             allbiotypes = c(file.path(data_dir, "bulk/darnell_2018/biotypes/SRR7073121"), 
                                             file.path(data_dir, "bulk/darnell_2018/biotypes/SRR7073122"), 
                                             file.path(data_dir, "bulk/darnell_2018/biotypes/SRR7073123"), 
                                             file.path(data_dir, "bulk/darnell_2018/biotypes/SRR7073124"), 
                                             file.path(data_dir, "bulk/darnell_2018/biotypes/SRR7073152"), 
                                             file.path(data_dir, "bulk/darnell_2018/biotypes/SRR7073153"), 
                                             file.path(data_dir, "bulk/darnell_2018/biotypes/SRR7073154")),
                             contamination = c(file.path(data_dir, "bulk/darnell_2018/contamination/SRR7073121_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"), 
                                               file.path(data_dir, "bulk/darnell_2018/contamination/SRR7073122_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"), 
                                               file.path(data_dir, "bulk/darnell_2018/contamination/SRR7073123_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"), 
                                               file.path(data_dir, "bulk/darnell_2018/contamination/SRR7073124_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"), 
                                               file.path(data_dir, "bulk/darnell_2018/contamination/SRR7073152_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"), 
                                               file.path(data_dir, "bulk/darnell_2018/contamination/SRR7073153_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"), 
                                               file.path(data_dir, "bulk/darnell_2018/contamination/SRR7073154_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz")),
                             reads = c(file.path(data_dir, "bulk/darnell_2018/reads/SRR7073121_Aligned.toTranscriptome.out.csv.predicted.csv.gz"), 
                                       file.path(data_dir, "bulk/darnell_2018/reads/SRR7073122_Aligned.toTranscriptome.out.csv.predicted.csv.gz"), 
                                       file.path(data_dir, "bulk/darnell_2018/reads/SRR7073123_Aligned.toTranscriptome.out.csv.predicted.csv.gz"), 
                                       file.path(data_dir, "bulk/darnell_2018/reads/SRR7073124_Aligned.toTranscriptome.out.csv.predicted.csv.gz"), 
                                       file.path(data_dir, "bulk/darnell_2018/reads/SRR7073152_Aligned.toTranscriptome.out.csv.predicted.csv.gz"), 
                                       file.path(data_dir, "bulk/darnell_2018/reads/SRR7073153_Aligned.toTranscriptome.out.csv.predicted.csv.gz"), 
                                       file.path(data_dir, "bulk/darnell_2018/reads/SRR7073154_Aligned.toTranscriptome.out.csv.predicted.csv.gz")),
                             stringsAsFactors = FALSE) 

tanenbaumSamples <- data.frame(names = c("tanenbaum_MG1-r1",
                                         "tanenbaum_MG1-r2",
                                         "tanenbaum_G2-r1",
                                         "tanenbaum_G2-r2",
                                         "tanenbaum_M-r1",
                                         "tanenbaum_M-r2"),
                               counts = c(file.path(data_dir, "bulk/tanenbaum_2015/counts/SRR1976443"), 
                                          file.path(data_dir, "bulk/tanenbaum_2015/counts/SRR1976444"), 
                                          file.path(data_dir, "bulk/tanenbaum_2015/counts/SRR1976445"), 
                                          file.path(data_dir, "bulk/tanenbaum_2015/counts/SRR1976446"), 
                                          file.path(data_dir, "bulk/tanenbaum_2015/counts/SRR1976447"), 
                                          file.path(data_dir, "bulk/tanenbaum_2015/counts/SRR1976448")),
                               allbiotypes = c(file.path(data_dir, "bulk/tanenbaum_2015/biotypes/SRR1976443"), 
                                               file.path(data_dir, "bulk/tanenbaum_2015/biotypes/SRR1976444"), 
                                               file.path(data_dir, "bulk/tanenbaum_2015/biotypes/SRR1976445"), 
                                               file.path(data_dir, "bulk/tanenbaum_2015/biotypes/SRR1976446"), 
                                               file.path(data_dir, "bulk/tanenbaum_2015/biotypes/SRR1976447"), 
                                               file.path(data_dir, "bulk/tanenbaum_2015/biotypes/SRR1976448")),
                               contamination = c(file.path(data_dir, "bulk/tanenbaum_2015/contamination/SRR1976443_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"), 
                                                 file.path(data_dir, "bulk/tanenbaum_2015/contamination/SRR1976444_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"), 
                                                 file.path(data_dir, "bulk/tanenbaum_2015/contamination/SRR1976445_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"), 
                                                 file.path(data_dir, "bulk/tanenbaum_2015/contamination/SRR1976446_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"), 
                                                 file.path(data_dir, "bulk/tanenbaum_2015/contamination/SRR1976447_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"), 
                                                 file.path(data_dir, "bulk/tanenbaum_2015/contamination/SRR1976448_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz")),
                               reads = c(file.path(data_dir, "bulk/tanenbaum_2015/reads/SRR1976443_Aligned.toTranscriptome.out.csv.gz"), 
                                         file.path(data_dir, "bulk/tanenbaum_2015/reads/SRR1976444_Aligned.toTranscriptome.out.csv.gz"), 
                                         file.path(data_dir, "bulk/tanenbaum_2015/reads/SRR1976445_Aligned.toTranscriptome.out.csv.gz"), 
                                         file.path(data_dir, "bulk/tanenbaum_2015/reads/SRR1976446_Aligned.toTranscriptome.out.csv.gz"), 
                                         file.path(data_dir, "bulk/tanenbaum_2015/reads/SRR1976447_Aligned.toTranscriptome.out.csv.gz"), 
                                         file.path(data_dir, "bulk/tanenbaum_2015/reads/SRR1976448_Aligned.toTranscriptome.out.csv.gz")),
                               stringsAsFactors = FALSE) 

ingoliaSamples <- data.frame(names = c("ingolia_HEK293T-low",
                                       "ingolia_HEK293T-med",
                                       "ingolia_HEK293T-high"),
                             counts = c(file.path(data_dir, "bulk/ingolia_2012/counts/SRR493747"),
                                        file.path(data_dir, "bulk/ingolia_2012/counts/SRR493748"),
                                        file.path(data_dir, "bulk/ingolia_2012/counts/SRR493749")),
                             allbiotypes = c(file.path(data_dir, "bulk/ingolia_2012/biotypes/SRR493747"),
                                             file.path(data_dir, "bulk/ingolia_2012/biotypes/SRR493748"),
                                             file.path(data_dir, "bulk/ingolia_2012/biotypes/SRR493749")),
                             contamination = c(file.path(data_dir, "bulk/ingolia_2012/contamination/SRR493747_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                               file.path(data_dir, "bulk/ingolia_2012/contamination/SRR493748_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                               file.path(data_dir, "bulk/ingolia_2012/contamination/SRR493749_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz")),
                             reads = c(file.path(data_dir, "bulk/ingolia_2012/reads/SRR493747_Aligned.toTranscriptome.out.csv.gz"),
                                       file.path(data_dir, "bulk/ingolia_2012/reads/SRR493748_Aligned.toTranscriptome.out.csv.gz"),
                                       file.path(data_dir, "bulk/ingolia_2012/reads/SRR493749_Aligned.toTranscriptome.out.csv.gz")),
                             stringsAsFactors = FALSE)

martinezSamples <- data.frame(names = c("martinez_293T-lores",
                                        "martinez_293T-medres",
                                        "martinez_293T-hires"),
                              counts = c(file.path(data_dir, "bulk/martinez_2020/counts/SRR8449566"),
                                         file.path(data_dir, "bulk/martinez_2020/counts/SRR8449567"),
                                         file.path(data_dir, "bulk/martinez_2020/counts/SRR8449568")),
                              allbiotypes = c(file.path(data_dir, "bulk/martinez_2020/biotypes/SRR8449566"),
                                              file.path(data_dir, "bulk/martinez_2020/biotypes/SRR8449567"),
                                              file.path(data_dir, "bulk/martinez_2020/biotypes/SRR8449568")),
                              contamination = c(file.path(data_dir, "bulk/martinez_2020/contamination/SRR8449566_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                                file.path(data_dir, "bulk/martinez_2020/contamination/SRR8449567_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                                file.path(data_dir, "bulk/martinez_2020/contamination/SRR8449568_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz")),
                              reads = c(file.path(data_dir, "bulk/martinez_2020/reads/SRR8449566_Aligned.toTranscriptome.out.csv.gz"),
                                        file.path(data_dir, "bulk/martinez_2020/reads/SRR8449567_Aligned.toTranscriptome.out.csv.gz"),
                                        file.path(data_dir, "bulk/martinez_2020/reads/SRR8449568_Aligned.toTranscriptome.out.csv.gz")),
                              stringsAsFactors = FALSE)
bulkSamples <- rbind(darnellSamples,
                     tanenbaumSamples,
                     ingoliaSamples, 
                     martinezSamples)

bulkReads <- canonical_pc %>%
  inner_join(fread(file.path(data_dir, "bulk/reads/mergedFilteredReads.csv.gz")), by = "transcript_id") %>%
  separate("id", c("filename", NA), remove = FALSE, sep = '\\.') %>%
  mutate(CB = paste(sample, filename, sep = "_")) %>%
  separate("CB", c("sample"), extra = "drop", remove = FALSE, sep = "_")

bulkCounts <- readMergeCounts(bulkSamples)
bulkBiotypes <- readMergeBiotypes(bulkSamples)

bulkMeta <- data.frame(row.names = colnames(bulkCounts),
                       CB = colnames(bulkCounts),
                       cds = bulkCounts["cds",],
                       utr3 = bulkCounts["utr3",],
                       utr5 = bulkCounts["utr5",],
                       cds_frac = bulkCounts["cds_frac",],
                       cell_total = bulkCounts["cell_total",],
                       stringsAsFactors = FALSE)
bulkCounts <- bulkCounts[!(rownames(bulkCounts) %in% c("cds","cds_frac","cell_total","utr3","utr5")),]


saveRDS(bulkReads, file.path(data_dir, "bulk/compiled/reads.rds"))
saveRDS(bulkCounts, file.path(data_dir, "bulk/compiled/counts.rds"))
saveRDS(bulkMeta, file.path(data_dir, "bulk/compiled/meta.rds"))
saveRDS(bulkBiotypes, file.path(data_dir, "bulk/compiled/biotypes.rds"))
saveRDS(bulkSamples, file.path(data_dir, "bulk/compiled/samples.rds"))

sessionInfo()
