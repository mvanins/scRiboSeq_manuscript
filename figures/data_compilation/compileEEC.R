#!/usr/bin/env Rscript
# Michael VanInsberghe 2021-03-15
# Compile, merge, and filter EE cells

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
library(mclust)
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

## mouse annot
mm_annot <- fread(file.path(data_dir, 'mmu_annotations.csv.gz'), stringsAsFactors = FALSE, data.table = FALSE)
mm_canonical_pc <- mm_annot %>%
  filter(set=="canonical", transcript_type == "protein_coding") %>%
  filter(chr != 'chrM') %>%
  filter(transcript_id != "ENST00000437139.7") %>%
  filter(!grepl("^PCDH", gene_name)) #%>%
  #filter(!grepl("^RPL", gene_name))

mm_codons <- fread(file.path(data_dir, 'mmu_codons.csv.gz'), stringsAsFactors = FALSE, data.table = FALSE)
mm_codons <- mm_codons %>%
  filter( transcript_id %in% mm_canonical_pc$transcript_id)

eeSamples <- data.frame(names=c("Dis01",
                                "Mid02",
                                "Pro02",
                                "Pro03",
                                "Pro05",
                                "Pro07",
                                "Pro08"),
                        counts = c(file.path(data_dir, "EE/counts/RPFv4-EE-Dist001"),
                                   file.path(data_dir, "EE/counts/RPFv4-EE-Mid002"),
                                   file.path(data_dir, "EE/counts/RPFv4-EE-Prox002"),
                                   file.path(data_dir, "EE/counts/RPFv4-EE-Prox003"),
                                   file.path(data_dir, "EE/counts/RPFv4-EE-Prox005"),
                                   file.path(data_dir, "EE/counts/RPFv4-EE-Prox007"),
                                   file.path(data_dir, "EE/counts/RPFv4-EE-Prox008")),
                        allcounts = c(file.path(data_dir, "EE/all/counts/RPFv4-EE-Dist001"),
                                      file.path(data_dir, "EE/all/counts/RPFv4-EE-Mid002"),
                                      file.path(data_dir, "EE/all/counts/RPFv4-EE-Prox002"),
                                      file.path(data_dir, "EE/all/counts/RPFv4-EE-Prox003"),
                                      file.path(data_dir, "EE/all/counts/RPFv4-EE-Prox005"),
                                      file.path(data_dir, "EE/all/counts/RPFv4-EE-Prox007"),
                                      file.path(data_dir, "EE/all/counts/RPFv4-EE-Prox008")),
                        biotypes = c(file.path(data_dir, "EE/biotypes/RPFv4-EE-Dist001"),
                                     file.path(data_dir, "EE/biotypes/RPFv4-EE-Mid002"),
                                     file.path(data_dir, "EE/biotypes/RPFv4-EE-Prox002"),
                                     file.path(data_dir, "EE/biotypes/RPFv4-EE-Prox003"),
                                     file.path(data_dir, "EE/biotypes/RPFv4-EE-Prox005"),
                                     file.path(data_dir, "EE/biotypes/RPFv4-EE-Prox007"),
                                     file.path(data_dir, "EE/biotypes/RPFv4-EE-Prox008")),
                        allbiotypes = c(file.path(data_dir, "EE/all/biotypes/RPFv4-EE-Dist001"),
                                        file.path(data_dir, "EE/all/biotypes/RPFv4-EE-Mid002"),
                                        file.path(data_dir, "EE/all/biotypes/RPFv4-EE-Prox002"),
                                        file.path(data_dir, "EE/all/biotypes/RPFv4-EE-Prox003"),
                                        file.path(data_dir, "EE/all/biotypes/RPFv4-EE-Prox005"),
                                        file.path(data_dir, "EE/all/biotypes/RPFv4-EE-Prox007"),
                                        file.path(data_dir, "EE/all/biotypes/RPFv4-EE-Prox008")),
                        contamination = c(file.path(data_dir, "EE/contamination/RPFv4-EE-Dist001_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                          file.path(data_dir, "EE/contamination/RPFv4-EE-Mid002_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                          file.path(data_dir, "EE/contamination/RPFv4-EE-Prox002_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                          file.path(data_dir, "EE/contamination/RPFv4-EE-Prox003_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                          file.path(data_dir, "EE/contamination/RPFv4-EE-Prox005_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                          file.path(data_dir, "EE/contamination/RPFv4-EE-Prox007_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                          file.path(data_dir, "EE/contamination/RPFv4-EE-Prox008_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz")),
                        facs = c(file.path(data_dir, "EE/facs/Dist_001_index.csv"),
                                 file.path(data_dir, "EE/facs/Mid_002_index.csv"),
                                 file.path(data_dir, "EE/facs/Prox_002_index.csv"),
                                 file.path(data_dir, "EE/facs/Prox_003_index.csv"),
                                 file.path(data_dir, "EE/facs/Prox_005_index.csv"),
                                 file.path(data_dir, "EE/facs/Prox_007_index.csv"),
                                 file.path(data_dir, "EE/facs/Prox_008_index.csv")),
                        reads = c(file.path(data_dir, "EE/reads/RPFv4-EE-Dist001_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                  file.path(data_dir, "EE/reads/RPFv4-EE-Mid002_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                  file.path(data_dir, "EE/reads/RPFv4-EE-Prox002_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                  file.path(data_dir, "EE/reads/RPFv4-EE-Prox003_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                  file.path(data_dir, "EE/reads/RPFv4-EE-Prox005_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                  file.path(data_dir, "EE/reads/RPFv4-EE-Prox007_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                  file.path(data_dir, "EE/reads/RPFv4-EE-Prox008_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz")),
                        stringsAsFactors = FALSE)

eeFacs <- readFACS(eeSamples, barcodes = barcodes) %>%
  mutate(green = X..488..530.40,
         red = X..561..585.29)
eeCounts <- readMergeCounts(eeSamples)
eeBiotypes <- readMergeBiotypes(eeSamples)

eeReads <- mm_canonical_pc %>%
  inner_join(fread(file.path(data_dir, "EE/reads/mergedFilteredReads.csv.gz")), by = "transcript_id")


eeMeta <- data.frame(row.names = colnames(eeCounts),
                     CB = colnames(eeCounts),
                     well = barcodes[gsub(".*_(.*?)$","\\1",colnames(eeCounts),perl=TRUE),"well"],
                     plate = gsub("^(.*?)[_].*","\\1",colnames(eeCounts),perl=TRUE),
                     cds = eeCounts["cds",],
                     utr3 = eeCounts["utr3",],
                     utr5 = eeCounts["utr5",],
                     cds_frac = eeCounts["cds_frac",],
                     cell_total = eeCounts["cell_total",],
                     sort_population = rep(NA,ncol(eeCounts)),
                     stringsAsFactors = FALSE)

eeMeta[ (grepl("Dis", eeMeta$plate)), "sort_population"] = "Distal"
eeMeta[ (grepl("Mid", eeMeta$plate)), "sort_population"] = "Medial"
eeMeta[ (grepl("Pro", eeMeta$plate)), "sort_population"] = "Proximal"
eeMeta[ (eeMeta$well %in% columnFeatures(21:24, LETTERS[15:16])), "sort_population"] = "NTC"


eeMeta <- eeMeta %>%
  full_join(eeFacs %>% select(-Well, -sample), by = "CB")
rownames(eeMeta) <- eeMeta$CB
eeMeta <- eeMeta[colnames(eeCounts),]

eeCounts <- eeCounts[!(rownames(eeCounts) %in% c("cds","cds_frac","cell_total","utr3","utr5")),]


set.seed(112)
doublet_cluster <- eeMeta %>%
  as_tibble() %>%
  filter(!is.na(FSC)) %>%
  select(CB, Trigger.Pulse.Width, starts_with("FSC", ignore.case = FALSE), starts_with("SSC", ignore.case = FALSE), X..405..460.50) %>%
  column_to_rownames("CB") %>%
  Mclust(., G=1:10)

eeMeta$facs_cluster <- NA
eeMeta[names(doublet_cluster$classification), "facs_cluster"] <- as.factor(doublet_cluster$classification)

plt_ssc_w <- ggplot(eeMeta, aes(x=SSC.Width, y = SSC, colour = as.factor(facs_cluster)))+
  geom_point(size=0.5)
save_plot(file.path(data_dir, "EE/compiled/SSC_threshold.pdf"), plt_ssc_w)

plt_ssc_w <- ggplot(eeMeta, aes(x=SSC.Width, y = SSC, colour = as.factor(facs_cluster)))+
  geom_point(data = eeMeta %>% select(SSC.Width, SSC), colour = "grey60", size = 0.5)+
  geom_point(size=0.5)+
  facet_wrap(~as.factor(facs_cluster))
save_plot(file.path(data_dir, "EE/compiled/SSC_threshold.pdf"), plt_ssc_w, base_width = 10, base_height = 8)

plt_ssc_w <- ggplot(eeMeta, aes(x=FSC.Width, y = FSC, colour = as.factor(facs_cluster)))+
  geom_point(data = eeMeta %>% select(FSC.Width, FSC), colour = "grey60", size = 0.5)+
  geom_point(size=0.5)+
  facet_wrap(~as.factor(facs_cluster))
save_plot(file.path(data_dir, "EE/compiled/FSC_threshold.pdf"), plt_ssc_w, base_width = 10, base_height = 8)



plt_fsc <- ggplot(eeMeta %>% arrange((as.factor(facs_cluster))), aes(x=FSC.Area, y = FSC, colour = as.factor(facs_cluster)))+
  geom_point(data = eeMeta %>% select(-facs_cluster), size = 0.5, colour = "grey60")+
  geom_point(size=0.5)+
  facet_wrap(~as.factor(facs_cluster))+
  coord_equal()
plt_ssc <- ggplot(eeMeta %>% arrange((as.factor(facs_cluster))), aes(x=SSC.Area, y = SSC, colour = as.factor(facs_cluster)))+
  geom_point(data = eeMeta %>% select(-facs_cluster), size = 0.5, colour = "grey60")+
  geom_point(size=0.5)+
  facet_wrap(~as.factor(facs_cluster))+
  coord_equal()
plt_ssc_w <- ggplot(eeMeta %>% arrange((as.factor(facs_cluster))), aes(x=SSC.Width, y = SSC, colour = as.factor(facs_cluster)))+
  geom_point(data = eeMeta %>% select(-facs_cluster), size = 0.5, colour = "grey60")+
  geom_point(size=0.5)+
  facet_wrap(~as.factor(facs_cluster))
plt_fsc_ssc <- ggplot(eeMeta %>% arrange((as.factor(facs_cluster))), aes(x=FSC.Area, y = SSC.Area, colour = as.factor(facs_cluster)))+
  geom_point(data = eeMeta %>% select(-facs_cluster), size = 0.5, colour = "grey60")+
  geom_point(size=0.5)+
  facet_wrap(~as.factor(facs_cluster))+
  coord_equal()
plt_fsc_dapi <- ggplot(eeMeta %>% arrange((as.factor(facs_cluster))), aes(x=log10(X..405..460.50), y = FSC, colour = as.factor(facs_cluster)))+
  geom_point(data = eeMeta %>% select(-facs_cluster), size = 0.5, colour = "grey60")+
  geom_point(size=0.5)+
  facet_wrap(~as.factor(facs_cluster))
plt_fsc_w <- ggplot(eeMeta %>% arrange((as.factor(facs_cluster))), aes(x=FSC.Width, y = FSC, colour = as.factor(facs_cluster)))+
  geom_point(data = eeMeta %>% select(-facs_cluster), size = 0.5, colour = "grey60")+
  geom_point(size=0.5)+
  facet_wrap(~as.factor(facs_cluster))
save_plot(file.path(data_dir, "EE/compiled/FACS.pdf"), plot_grid(plt_fsc, plt_ssc, plt_ssc_w, plt_fsc_ssc, plt_fsc_w, plt_fsc_dapi, ncol = 2 ), base_width = 30, base_height = 30)


doublet_clusters <- c(1, 2, 4, 5, 10)
clusterGated <- eeMeta %>%
  filter(!is.na(FSC)) %>%
  filter(!(facs_cluster %in% doublet_clusters)) %>%
  select(CB) %>%
  pull()

eeMeta$ssc_gated <- "doublet"
eeMeta[clusterGated, "ssc_gated"] <- "single"
saveRDS(eeMeta, file.path(data_dir, "EE/compiled/allmeta.rds"))

plateThresh <- data.frame(plate = c("Dis01", "Mid02", "Pro02", "Pro03", "Pro05", "Pro07", "Pro08"),
                          min_reads = c(2000, 2000,    2000,    1800,    2000,    2000,    2000),
                          min_frac = c(.75, .75, .75, .75, .75, .75, .75),
                          stringsAsFactors = FALSE)
plt_cd_totalQC <- ggplot(eeMeta %>% filter(ssc_gated == "single"), aes(x=log10(cell_total), y=cds_frac, colour = as.factor(facs_cluster)))+
  geom_point(size=0.8)+
  geom_vline(data = plateThresh, aes(xintercept = log10(min_reads)), colour = 'red', linetype = 'dashed')+
  geom_hline(data = plateThresh, aes(yintercept = min_frac), colour = 'red', linetype = 'dashed')+
  facet_wrap(~plate, nrow = 1)
save_plot(file.path(data_dir, "EE/compiled/QC_thresholds.pdf"), plt_cd_totalQC, base_width = 17, base_height = 4.5)

eeCells <- pull(eeMeta %>% 
                 inner_join(plateThresh, by = "plate") %>%
                 filter(!is.na(SSC)) %>%
                 filter(ssc_gated == "single") %>%
                 filter(cds_frac >= min_frac,
                        cell_total > min_reads) %>%
                 filter(!(sort_population %in% c("NTC"))) %>%
                 select(CB))

eeReads <- eeReads %>%
  filter(CB %in% eeCells)
eeCounts <- eeCounts[,eeCells]
eeBiotypes <- eeBiotypes[,eeCells]
eeFacs <- eeFacs %>%
  filter(CB %in% eeCells)
eeMeta <- eeMeta[eeCells,]
eeMeta$facs_cluster <- as.factor(eeMeta$facs_cluster)

saveRDS(eeReads, file.path(data_dir, "EE/compiled/reads.rds"))
saveRDS(eeCounts, file.path(data_dir, "EE/compiled/counts.rds"))
saveRDS(eeMeta, file.path(data_dir, "EE/compiled/meta.rds"))
saveRDS(eeBiotypes, file.path(data_dir, "EE/compiled/biotypes.rds"))
saveRDS(eeSamples, file.path(data_dir, "EE/compiled/samples.rds"))
