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
hekSamples <- data.frame(names = c("HEK293T",
                                   "HEK293T-Starv1",
                                   "HEK293T-Starv2"),
                         counts = c(file.path(data_dir, "HEK293T/counts/RPFv4-HEK293T"),
                                    file.path(data_dir, "HEK293T/counts/RPFv4-HEK293T-Starv1"),
                                    file.path(data_dir, "HEK293T/counts/RPFv4-HEK293T-Starv2")),
                         allcounts = c(file.path(data_dir, "HEK293T/all/counts/RPFv4-HEK293T"),
                                       file.path(data_dir, "HEK293T/all/counts/RPFv4-HEK293T-Starv1"),
                                       file.path(data_dir, "HEK293T/all/counts/RPFv4-HEK293T-Starv2")),
                         biotypes = c(file.path(data_dir, "HEK293T/biotypes/RPFv4-HEK293T"),
                                      file.path(data_dir, "HEK293T/biotypes/RPFv4-HEK293T-Starv1"),
                                      file.path(data_dir, "HEK293T/biotypes/RPFv4-HEK293T-Starv2")),
                         allbiotypes = c(file.path(data_dir, "HEK293T/all/biotypes/RPFv4-HEK293T"),
                                         file.path(data_dir, "HEK293T/all/biotypes/RPFv4-HEK293T-Starv1"),
                                         file.path(data_dir, "HEK293T/all/biotypes/RPFv4-HEK293T-Starv2")),
                         contamination = c(file.path(data_dir, "HEK293T/contamination/RPFv4-HEK293T_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                           file.path(data_dir, "HEK293T/contamination/RPFv4-HEK293T-Starv1_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                           file.path(data_dir, "HEK293T/contamination/RPFv4-HEK293T-Starv2_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz")),
                         reads = c(file.path(data_dir, "HEK293T/reads/RPFv4-HEK293T_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                   file.path(data_dir, "HEK293T/reads/RPFv4-HEK293T-Starv1_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                   file.path(data_dir, "HEK293T/reads/RPFv4-HEK293T-Starv2_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz")),
                         stringsAsFactors = FALSE)


hekReads <- canonical_pc %>%
  inner_join(fread(file.path(data_dir, "HEK293T/reads/mergedFilteredReads.csv.gz")), by = "transcript_id")
hekCounts <- readMergeCounts(hekSamples)
hekBiotypes <- readMergeBiotypes(hekSamples)
hekMeta <- data.frame(row.names = colnames(hekCounts),
                      CB = colnames(hekCounts),
                      well = barcodes[gsub(".*_(.*?)$","\\1",colnames(hekCounts),perl=TRUE),"well"],
                      plate = gsub("^(.*?)[_].*","\\1",colnames(hekCounts),perl=TRUE),
                      cds = hekCounts["cds",],
                      utr3 = hekCounts["utr3",],
                      utr5 = hekCounts["utr5",],
                      cds_frac = hekCounts["cds_frac",],
                      cell_total = hekCounts["cell_total",],
                      sort_population = rep("Rich",ncol(hekCounts)),
                      stringsAsFactors = FALSE)
hekMeta[ (grepl("Starv", hekMeta$plate) & hekMeta$well %in% columnFeatures(1:5,LETTERS[1:16])), "sort_population" ] = "Arg3h"
hekMeta[ (grepl("Starv", hekMeta$plate) & hekMeta$well %in% columnFeatures(6:10,LETTERS[1:16])), "sort_population" ] = "Arg6h"
hekMeta[ (grepl("Starv", hekMeta$plate) & hekMeta$well %in% columnFeatures(11:15,LETTERS[1:16])), "sort_population" ] = "Leu3h"
hekMeta[ (grepl("Starv", hekMeta$plate) & hekMeta$well %in% columnFeatures(16:20,LETTERS[1:16])), "sort_population" ] = "Leu6h"
hekMeta[ (grepl("Starv", hekMeta$plate) & hekMeta$well %in% columnFeatures(21:23,LETTERS[1:16])), "sort_population" ] = "Rich" 
hekMeta[hekMeta$well %in% columnFeatures(24, LETTERS[1:16]), "sort_population"] = "NTC"

hekCounts <- hekCounts[!(rownames(hekCounts) %in% c("cds","cds_frac","cell_total","utr3","utr5")),]


saveRDS(hekMeta, file.path(data_dir, "HEK293T/compiled/allmeta.rds"))

hekCells <- pull(hekMeta %>% filter(sort_population != "NTC") %>% filter(cds_frac >= 0.75, cell_total > 2300) %>% select(CB))
hekReads <- hekReads %>%
  filter(CB %in% hekCells)
hekCounts <- hekCounts[,hekCells]
hekBiotypes <- hekBiotypes[,hekCells]
hekMeta <- hekMeta[hekCells,]

saveRDS(hekReads, file.path(data_dir, "HEK293T/compiled/reads.rds"))
saveRDS(hekCounts, file.path(data_dir, "HEK293T/compiled/counts.rds"))
saveRDS(hekMeta, file.path(data_dir, "HEK293T/compiled/meta.rds"))
saveRDS(hekBiotypes, file.path(data_dir, "HEK293T/compiled/biotypes.rds"))
saveRDS(hekSamples, file.path(data_dir, "HEK293T/compiled/samples.rds"))
