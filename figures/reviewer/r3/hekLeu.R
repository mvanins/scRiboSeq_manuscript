#!/usr/bin/env Rscript
# Michael VanInsberghe 2021-04-20
# Figures for Response to R3 leucine pausing

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(data.table)
library(Matrix)
library(Matrix.utils)
library(Seurat)
library(ggplot2)
library(lemon)
library(egg)
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

figurePrefix <- "HEKleu/HEKleu"

data_dir <- "../scRiboSeq/data"

# barcodes
barcodes <- read.table('/hpc/hub_oudenaarden/mvanins/local/lib/barcodes/miRv6_barcodemap.tsv',
                       sep="\t",
                       header = FALSE,
                       stringsAsFactors = FALSE,
                       row.names = NULL,
                       col.names = c("CB","number","well"))
rownames(barcodes) <- barcodes$CB

annot <- fread(file.path(data_dir, 'hsa_annotations.csv.gz'), stringsAsFactors = FALSE, data.table = FALSE)
canonical_pc <- annot %>%
  filter(set=="canonical", transcript_type == "protein_coding") %>%
  filter(transcript_id != "ENST00000437139.7") %>%
  filter(!grepl("^PCDH", gene_name))

codons <- fread(file.path(data_dir, 'hsa_codons.csv.gz'), stringsAsFactors = FALSE, data.table = FALSE)
codons <- codons %>%
  filter( transcript_id %in% canonical_pc$transcript_id)


hekReads <- readRDS(file.path(data_dir, "HEK293T/compiled/reads.rds"))
hekMeta <- readRDS(file.path(data_dir, "HEK293T/compiled/meta.rds"))
hekCounts <- readRDS(file.path(data_dir, "HEK293T/compiled/counts.rds"))
hekBiotypes <- readRDS(file.path(data_dir, "HEK293T/compiled/biotypes.rds"))

# select only starvation experiment
hekMeta <- hekMeta %>%
  filter(plate %in% c("HEK293T-Starv1", "HEK293T-Starv2"))
hekReads <- hekReads %>%
  filter(CB %in% pull(hekMeta %>% select(CB) %>% distinct()))
hekCounts <- hekCounts[, rownames(hekMeta)]
hekBiotypes <- hekBiotypes[, rownames(hekMeta)]

hekLong <- as.data.frame(t(hekCounts))
hekLong$CB <- rownames(hekLong)
hekLong <- hekLong %>%
  pivot_longer(cols = starts_with("ENSG"), names_to = "gene_id", values_to = "counts") 

coltreat <- c("#a6cee3", "#1f78b4",
              "#b2df8a", "#33a02c",
              "grey50")
names(coltreat) <- c("Arg3h", "Arg6h",
                     "Leu3h", "Leu6h",
                     "Rich")



# cluster etc
ENSG_to_name <- data.frame(count_names = rownames(hekCounts), stringsAsFactors = FALSE) %>%
  inner_join(annot, by = c("count_names" = "gene_id")) %>% 
  mutate(gene_out = paste(count_names, gene_name, sep = "-")) %>%
  select(count_names, gene_out) %>%
  distinct() %>%
  column_to_rownames(var = "count_names")
RPFcounts <- hekCounts[!(rownames(hekCounts) %in% c("cds", "cds_frac", "cell_total", "utr3", "utr5")),]
rownames(RPFcounts) <- ENSG_to_name[rownames(RPFcounts),"gene_out"]

P <- CreateSeuratObject(counts = RPFcounts,
                        project = "RPF_starv",
                        assay = "RNA",
                        min.cells = 0,
                        min.features = 0,
                        names.field = 1,
                        names.delim = "-",
                        meta.data = hekMeta)

P <- NormalizeData(P)

P <- FindVariableFeatures(P, selection.method = "vst", nfeatures = 1000)

P <- ScaleData(P, features = rownames(P))
P <- RunPCA(P, features = VariableFeatures(object = P))

P <- JackStraw(P, num.replicate = 100)
P <- ScoreJackStraw(P, dims = 1:20)

plt_jackstraw <- JackStrawPlot(P,dims=1:20)
plt_elbow <- ElbowPlot(P)

save_plot(paste0(figurePrefix,"_qc_plt_dimscores.pdf"), plot_grid(plt_jackstraw, plt_elbow, nrow=1), base_width = 10, base_height = 5)

ndims <- 5 
P <- RunUMAP(P, dims = 1:ndims, verbose = FALSE, seed.use = 121)
P <- FindNeighbors(P, dims = 1:ndims)
P <- FindClusters(P, resolution = 0.4)

plt_umap_plate <- DimPlot(P, reduction = "umap",group.by = "plate",pt.size=1)
plt_umap_cluster <- DimPlot(P, reduction = "umap",pt.size=1)
plt_umap_treatment <- DimPlot(P, reduction = 'umap', group.by = "sort_population", pt.size=1)


save_plot(paste0(figurePrefix,"_plt_umaps_check.pdf"), plot_grid(plt_umap_plate+coord_equal(), plt_umap_treatment+coord_equal() ,plt_umap_cluster+coord_equal(), nrow=2), base_width = 10, base_height = 10)

hekMeta$seurat_clusters <- P$seurat_clusters[rownames(hekMeta)]
hekMeta$UMAP_1 <- as.data.frame(Embeddings(P, reduction="umap"))[rownames(hekMeta),"UMAP_1"]
hekMeta$UMAP_2 <- as.data.frame(Embeddings(P, reduction="umap"))[rownames(hekMeta),"UMAP_2"]



#### is pausing only seen in cell-cycle genes?
## first set up pausing data

codon_reads <- hekReads %>%
  select(-starts_with("base")) %>%
  rowwise() %>% # expand each read to each site (rowwise, mutate, unnest)
    mutate(site = list(c("em2","em1","e","p","a","ap1","ap2")),
           site_position = list(c(psite-9, psite-6, psite-3, psite, psite+3, psite+6, psite+9))) %>%
  unnest(cols = c(site,site_position)) %>%
  select(-psite) %>% # remove p-site column
  inner_join(codons, by = c("transcript_id" = "transcript_id", "site_position" = "codon_position")) # join codon information. Will only return site_position that exactly match codon_position, i.e., in frame

site_order <- data.frame(site = c("em2", "em1", "e", "p", "a", "ap1", "ap2"),
                         site_order = c(1, 2, 3, 4, 5, 6, 7),
                         stringsAsFactors = FALSE)

baseline_counts <- codon_reads %>%
  filter(!(three %in% c("Ter"))) %>%
  mutate(codon_display = paste(gsub("T", "U", codon), three, sep = "_")) %>%
  group_by(CB, codon_display, site) %>%
  summarize(n = n()) %>%
  ungroup()

total_baseline <- baseline_counts %>%
  group_by(CB, codon_display) %>%
  summarize(codon_counts = sum(n)) %>%
  ungroup() %>%
  group_by(CB) %>%
  mutate(cell_total_fraction = codon_counts/sum(codon_counts)) %>%
  ungroup() %>%
  group_by(codon_display) %>%
  summarize(baseline_total_fraction = mean(cell_total_fraction)) %>%
  ungroup()

site_baseline <- baseline_counts %>%
  group_by(CB, site) %>%
  mutate(cell_site_fraction = n/sum(n)) %>%
  ungroup() %>%
  group_by(codon_display, site) %>%
  summarize(baseline_site_fraction = mean(cell_site_fraction)) %>%
  ungroup()

site_change <- codon_reads %>%
  filter(!(three %in% c("Ter"))) %>%
  mutate(codon_display = paste(gsub("T", "U", codon), three, sep = "_")) %>%
  group_by(CB, site, codon_display) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(CB, site) %>% 
  mutate(cell_site_fraction = n/sum(n)) %>%
  ungroup() %>%
  group_by(CB) %>%
  mutate(cell_total_fraction = n/sum(n)) %>%
  ungroup() %>%
  inner_join(hekMeta %>% select(CB, seurat_clusters, sort_population, plate), by = "CB") %>%
  inner_join(total_baseline, by = "codon_display") %>%
  inner_join(site_baseline, by = c("codon_display" = "codon_display", "site" = "site")) %>%
  mutate(total_fc = cell_total_fraction/(baseline_total_fraction/nrow(site_order)),
         site_fc = cell_site_fraction/baseline_site_fraction) %>%
  inner_join(site_order, by = "site") %>%
  mutate(site = fct_reorder(site, site_order))



## define thresholds to find which cells pausing

pausing_signal_Leu <- site_change %>%
  filter(codon_display %in% c("UUA_Leu")) %>%
  filter(site %in% c("p", "a")) %>%
  group_by(CB, codon_display) %>%
  summarize(pausing = mean(site_fc)) %>%
  ungroup() %>%
  inner_join(hekMeta %>% select(CB, sort_population, seurat_clusters), by = "CB")

pausing_signal_Arg <- site_change %>%
  filter(codon_display %in% c("CGC_Arg", "CGU_Arg")) %>%
  filter(site %in% c("e", "p", "a")) %>%
  group_by(CB, codon_display) %>%
  summarize(pausing = mean(site_fc)) %>%
  ungroup() %>%
  inner_join(hekMeta %>% select(CB, sort_population, seurat_clusters), by = "CB")

pausing_signal <- bind_rows(pausing_signal_Leu,
                            pausing_signal_Arg)


pausing_threshold <- pausing_signal %>%
  filter(sort_population == "Rich") %>%
  summarize(mean_pausing = mean(pausing),
            sd_pausing = sd(pausing)) %>%
  mutate(threshold = mean_pausing + 3.2905*sd_pausing) 
  #mutate(threshold = mean_pausing + 4.5*sd_pausing) 

pausing_threshold <- bind_rows(pausing_signal_Leu %>%
                                 filter(sort_population %in% c("Rich", "Arg3h", "Arg6h")),
                               pausing_signal_Arg %>%
                                 filter(sort_population %in% c("Rich", "Leu3h", "Leu6h"))) %>%
  summarize(mean_pausing = mean(pausing),
            sd_pausing = sd(pausing)) %>%
  mutate(threshold = mean_pausing + 4*sd_pausing) 
  

plt_hist_codon <- ggplot(pausing_signal,
                         aes(x=log2(pausing), fill = sort_population))+
  geom_histogram(bins = 50)+
  scale_fill_manual(values = as.character(coltreat))+
  facet_wrap(~codon_display+sort_population, ncol = 5)+
  geom_vline(xintercept = log2(pausing_threshold$threshold), colour = "red", linetype = "dashed")
save_plot(paste0(figurePrefix, "_hist_pausing.pdf"), plt_hist_codon, base_width = 12, base_height = 8)


hekMeta <- hekMeta %>%
  select(-starts_with("ps_")) %>%
  select(-starts_with("pc_")) %>%
  inner_join(pausing_signal %>%
             select(CB, codon_display, pausing) %>%
             pivot_wider(names_from = codon_display, values_from = pausing, names_prefix = "ps_"), by = "CB") %>%
  mutate(pc_UUA_Leu = ps_UUA_Leu >= pausing_threshold$threshold,
         pc_CGC_Arg = ps_CGC_Arg >= pausing_threshold$threshold,
         pc_CGU_Arg = ps_CGU_Arg >= pausing_threshold$threshold)


## set Leu pausing as identity
P <- AddMetaData(P,
                 hekMeta %>% 
                 mutate(leu_pausing = case_when(pc_UUA_Leu == TRUE ~ "leu_pausing",
                                                pc_UUA_Leu == FALSE ~ "rest")) %>%
                 select(CB,leu_pausing) %>%
                 deframe(),
               col.name = "leu_pausing")
Idents(object = P) <- P@meta.data$leu_pausing
leu_markers <- FindAllMarkers(P, min.pct = 0.25, logfc.threshold = 0.25)

leu_markers %>%
  filter(p_val_adj < 0.01) %>%
  count(cluster)

