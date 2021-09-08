#!/usr/bin/env Rscript
# Michael VanInsberghe 2021-04-07
# Figures for Response to R2 #5
# HEK starvation sequence context

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

figurePrefix <- "HEKseq/HEKseq"

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
  #filter(chr != 'chrM') %>%
  filter(transcript_id != "ENST00000437139.7") %>%
  filter(!grepl("^PCDH", gene_name)) #%>%
  #filter(!grepl("^RPL", gene_name))

codons <- fread(file.path(data_dir, 'hsa_codons.csv.gz'), stringsAsFactors = FALSE, data.table = FALSE)
codons <- codons %>%
  filter( transcript_id %in% canonical_pc$transcript_id)


histone_genes <- read.table('histones_group-384.csv', sep = ',', header = TRUE, stringsAsFactors = FALSE)


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


## are histone genes enriched for any codons?
codon_group <- codons %>%
  filter(!(three %in% c("Ter"))) %>%
  mutate(codon_display = paste(three, gsub("T", "U", codon), sep = "_")) %>%
  filter(nchar(codon_display) == 7) %>%
  group_by(transcript_id, codon_display) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  complete(transcript_id, codon_display, fill = list(n=0)) %>%
  group_by(transcript_id) %>%
  mutate(gene_total = sum(n)) %>%
  ungroup() %>%
  inner_join(canonical_pc %>% select(transcript_id, gene_name), by = "transcript_id") %>%
  mutate(codon_group = case_when(gene_name %in% histone_genes$Approved.symbol ~ "histone",
                                 TRUE ~ "rest"))

plt_codon_enrichment <- codon_group %>%
  separate(codon_display, into = c("three", "codon"), remove = FALSE) %>%
  filter(three %in% c("Arg", "Leu")) %>%
  ggplot(aes(x=codon, y=n/gene_total, fill = codon_group))+
  geom_boxplot(outlier.size = 0.8, outlier.stroke = 0)+
  facet_wrap(~three, nrow = 1, scales = "free_x")+
  scale_fill_manual(values = brewer.pal(3, "Greys")[c(2,3)])+
  scale_colour_manual(values = brewer.pal(3, "Greys")[c(2,3)])+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
save_plot(paste0(figurePrefix, "_histone_codon_enrichment.pdf"), plt_codon_enrichment, base_width = 6, base_height = 5)



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

# plot variable features with and without labels
plt_variable <- VariableFeaturePlot(P)
plt_variable_labelled <- LabelPoints(plot = plt_variable, points = head(VariableFeatures(P),10), repel = TRUE)
#save_plot(paste0(figurePrefix, "_variable_features.pdf"), plot_grid(plt_variable, plt_variable_labelled, nrow=1), base_width=16, base_height = 6)

P <- ScaleData(P, features = rownames(P))
P <- RunPCA(P, features = VariableFeatures(object = P))

P <- JackStraw(P, num.replicate = 100)
P <- ScoreJackStraw(P, dims = 1:20)

plt_jackstraw <- JackStrawPlot(P,dims=1:20)
plt_elbow <- ElbowPlot(P)

save_plot(paste0(figurePrefix,"_qc_plt_dimscores.pdf"), plot_grid(plt_jackstraw, plt_elbow, nrow=1), base_width = 10, base_height = 5)

ndims <- 5 
P <- RunUMAP(P, dims = 1:ndims, verbose = FALSE, seed.use = 121)#, n.epochs = 3000)
P <- FindNeighbors(P, dims = 1:ndims)
P <- FindClusters(P, resolution = 0.4)

plt_umap_plate <- DimPlot(P, reduction = "umap",group.by = "plate",pt.size=1)
plt_umap_cluster <- DimPlot(P, reduction = "umap",pt.size=1)
plt_umap_treatment <- DimPlot(P, reduction = 'umap', group.by = "sort_population", pt.size=1)


save_plot(paste0(figurePrefix,"_plt_umaps_check.pdf"), plot_grid(plt_umap_plate+coord_equal(), plt_umap_treatment+coord_equal() ,plt_umap_cluster+coord_equal(), nrow=2), base_width = 10, base_height = 10)

hekMeta$seurat_clusters <- P$seurat_clusters[rownames(hekMeta)]
hekMeta$UMAP_1 <- as.data.frame(Embeddings(P, reduction="umap"))[rownames(hekMeta),"UMAP_1"]
hekMeta$UMAP_2 <- as.data.frame(Embeddings(P, reduction="umap"))[rownames(hekMeta),"UMAP_2"]





### GSEA
library(fgsea)
library(msigdbr)
library(presto)

# find markers


# update ENSG so we can get proper gene symbols again
ENSG_to_name$gene_id <- rownames(ENSG_to_name)
ENSG_to_name <- ENSG_to_name %>%
  inner_join(annot %>% select(gene_id, gene_name) %>% distinct(), by = "gene_id")

presmarkers <- wilcoxauc(P, 'seurat_clusters')


## Reactome
fgsea_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  split(x = .$gene_symbol, f = .$gs_name)

fgseaReactRes <- NULL
set.seed(115)
for(cl in unique(presmarkers$group)){
  cluster_ranks <- presmarkers %>%
    filter(group == cl) %>%
    inner_join(ENSG_to_name %>% select(gene_out, gene_name), by = c("feature" = "gene_out")) %>%
    arrange(desc(auc)) %>%
    select(gene_name, auc) %>%
    deframe()
  fgseaReactRes <- bind_rows(fgseaReactRes, 
                    fgsea(fgsea_reactome, stats = cluster_ranks, nperm = 1e5) %>%
                      as_tibble() %>%
                      mutate(cluster = cl))
}



plt_gsea <- fgseaReactRes %>%
  filter(padj < 0.01) %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = NES) %>%
  ungroup() %>%
  mutate(cluster_pathway = paste(cluster, pathway, sep = "--")) %>%
  mutate(cluster_pathway = fct_reorder(cluster_pathway, NES)) %>%
  mutate(logpadj = log10(padj)) %>%
  select(cluster_pathway, logpadj, NES, cluster) %>%
  pivot_longer(cols = c("logpadj", "NES"), values_to = "value", names_to = "stat") %>%
  ggplot(aes(x = cluster_pathway, y = value, fill = stat))+
  geom_bar(stat = "identity", position = "identity")+
  scale_fill_manual(values = brewer.pal(3, "Greys")[2:3])+
  scale_y_continuous(breaks = scales::pretty_breaks(n=10))+
  facet_wrap(~cluster, scales = "free_y")+
  coord_flip()+
  labs(x = "Reactome term",
       y = "Normalized Enrichment Score")
save_plot(paste0(figurePrefix, "_gsea_cluster_react.pdf"), plt_gsea, base_width = 33, base_height = 10)



## What fraction of RPFs/cell align to histone genes?

histone_fraction <- hekLong %>%
  mutate(gene_group = case_when(gene_id %in% pull(canonical_pc %>% filter(gene_name %in% histone_genes$Approved.symbol) %>% select(gene_id)) ~ "histone",
                                TRUE ~ "rest")) %>%
  group_by(CB, gene_group) %>%
  summarize(cell_counts = sum(counts)) %>%
  ungroup() %>%
  group_by(CB) %>% 
  mutate(cell_frac = cell_counts/sum(cell_counts)) %>%
  ungroup() %>%
  inner_join(hekMeta %>% select(CB, seurat_clusters), by = "CB")

plt_histone_frac <- histone_fraction %>%
  filter(gene_group == "histone") %>%
  ggplot(aes(x=seurat_clusters, y=cell_frac))+
  geom_boxplot()#+
  #stat_summary(fun.data = function(x){ return( c(y=0.1, label = length(x)) ) }, geom = "text", fun = max )
save_plot(paste0(figurePrefix, "_cluster_histone_frac.pdf"), plt_histone_frac, base_width = 3, base_height = 5)

histone_fraction %>%
  filter(gene_group == "histone") %>%
  group_by(seurat_clusters) %>%
  summarize(mean_histone = mean(cell_frac),
            sd_histone = sd(cell_frac),
            n_cells = n(),
            se_histone = sd(cell_frac)/sqrt(n())) %>%
  ungroup()


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
                         #aes(x=log2(pausing), fill = cl))+
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




  
# make figures, highlight cell-cycle genes  
cell_cycle_reactome <- c('REACTOME_CELL_CYCLE',
                         'REACTOME_CELL_CYCLE_CHECKPOINTS',
                         'REACTOME_CELL_CYCLE_MITOTIC',
                         'REACTOME_CYCLIN_A_B1_B2_ASSOCIATED_EVENTS_DURING_G2_M_TRANSITION',
                         'REACTOME_CYCLIN_A_CDK2_ASSOCIATED_EVENTS_AT_S_PHASE_ENTRY',
                         'REACTOME_CYCLIN_D_ASSOCIATED_EVENTS_IN_G1',
                         'REACTOME_DNA_REPAIR',
                         'REACTOME_G0_AND_EARLY_G1',
                         'REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS',
                         'REACTOME_G1_S_SPECIFIC_TRANSCRIPTION',
                         'REACTOME_G2_M_CHECKPOINTS',
                         'REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT',
                         'REACTOME_G2_M_DNA_REPLICATION_CHECKPOINT',
                         'REACTOME_G2_PHASE',
                         'REACTOME_M_PHASE',
                         'REACTOME_MITOTIC_G1_PHASE_AND_G1_S_TRANSITION',
                         'REACTOME_MITOTIC_G2_G2_M_PHASES',
                         'REACTOME_MITOTIC_METAPHASE_AND_ANAPHASE',
                         'REACTOME_MITOTIC_PROMETAPHASE',
                         'REACTOME_MITOTIC_PROPHASE',
                         'REACTOME_MITOTIC_SPINDLE_CHECKPOINT',
                         'REACTOME_MITOTIC_TELOPHASE_CYTOKINESIS',
                         'REACTOME_S_PHASE')

reactome_cc_genes <- pull(msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  filter(gs_name %in% cell_cycle_reactome) %>%
  select(gene_symbol) %>%
  distinct())


gene_codon_background <- codon_reads %>%
  filter(site == "p") %>%
  group_by(gene_id, gene_name, codon) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(gene_name) %>%
  mutate(codon_frac = n/sum(n),
         gene_total = sum(n)) %>%
  ungroup() %>%
  mutate(gene_mean = gene_total / pull(codon_reads %>% select(CB) %>% distinct() %>% summarize(n=n())) ) # divided by number of cells



## this until figure is different for every codon

highlight_codon <- "CGT" # this is used for the plot, groupings are based on pc_CODON_Three groups
cluster_change <- codon_reads %>%
  filter(site == "p") %>%
  filter(gene_id %in% pull(gene_codon_background %>% filter(gene_mean > 1) %>% select(gene_id) %>% distinct())) %>%
  #filter(gene_id %in% pull(gene_codon_background %>% select(gene_id) %>% distinct())) %>%
  inner_join(hekMeta %>% group_by(seurat_clusters, pc_CGU_Arg) %>% mutate(cells_per_cluster = n()) %>% ungroup() %>% select(CB, seurat_clusters, cells_per_cluster, pc_CGU_Arg), by = "CB") %>%
  group_by(gene_id, gene_name, codon, seurat_clusters, cells_per_cluster, pc_CGU_Arg) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(gene_name, seurat_clusters, pc_CGU_Arg) %>%
  mutate(cluster_gene_frac = n/sum(n),
         cluster_gene_total = sum(n)) %>%
  ungroup() %>% 
  mutate(cluster_gene_mean = cluster_gene_total / cells_per_cluster) %>%
  inner_join(gene_codon_background %>% select(gene_id, codon, codon_frac), by = c("gene_id" = "gene_id", "codon" = "codon")) %>%
  mutate(lfc_change = log2(cluster_gene_frac/codon_frac))
shared_genes <- pull(cluster_change %>%
  group_by(gene_id, seurat_clusters, pc_CGU_Arg) %>%
  summarize(asites_per_cluster = sum(n)) %>%
  ungroup() %>%
  complete(gene_id, seurat_clusters, pc_CGU_Arg) %>%
  mutate(asites_per_cluster = replace_na(asites_per_cluster,0)) %>%
  inner_join(hekMeta %>% group_by(seurat_clusters, pc_CGU_Arg) %>% mutate(cells_per_cluster = n()) %>% ungroup() %>% select(seurat_clusters, cells_per_cluster, pc_CGU_Arg) %>% distinct(), by = c("seurat_clusters", "pc_CGU_Arg")) %>%
  mutate(mean_asites = asites_per_cluster / cells_per_cluster) %>%
  filter(mean_asites > 2.5) %>%
  group_by(gene_id) %>%
  summarize(detected_in_clusters = n()) %>%
  ungroup() %>%
  filter(detected_in_clusters > 2) %>%
  select(gene_id) %>%
  distinct())
highlight_genes <- c("GAPDH", "H3C2", "PRKDC")
gene_codon_change_data <- cluster_change %>%
  mutate(gene_class = case_when(gene_name %in% histone_genes$Approved.symbol ~ "histone",
                                gene_name %in% reactome_cc_genes ~ "cell_cycle",
                             TRUE ~ "rest")) %>%
  filter(codon == highlight_codon) %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name %in% highlight_genes ~ gene_name,
                                  TRUE ~ "")) %>%
  arrange(desc(gene_class))
plt_cluster_change_gaa_facet <- ggplot(gene_codon_change_data, aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = gene_class, label = gene_display ))+
  geom_point(size=0.5, alpha = 1)+
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~seurat_clusters+pc_CGU_Arg, nrow=2)+
  geom_text(data = gene_codon_change_data %>% select(seurat_clusters, cells_per_cluster, pc_CGU_Arg) %>% distinct() %>% mutate(gene_class = "rest"), 
            mapping = aes(x = -Inf, y=-Inf, label = cells_per_cluster),
            hjust = -0.1,
            vjust = -1)+
  geom_text_repel(segment.color = "black", box.padding = 0.3, min.segment.length = 0, max.iter = 10000, seed=112, size=3, point.padding = 0.3)+
  xlab("log10 mean cluster expression")+ylab(paste0("log2 f.c. ", highlight_codon, " P-site FOC"))+
  theme(strip.background = element_blank())
save_plot(paste0(figurePrefix, "_plt_cluster_change_", highlight_codon, ".pdf"), plt_cluster_change_gaa_facet, base_width = 8, base_height =5)



highlight_codon <- "CGC" # this is used for the plot, groupings are based on pc_CODON_Three groups
cluster_change <- codon_reads %>%
  filter(site == "p") %>%
  filter(gene_id %in% pull(gene_codon_background %>% filter(gene_mean > 1) %>% select(gene_id) %>% distinct())) %>%
  inner_join(hekMeta %>% group_by(seurat_clusters, pc_CGC_Arg) %>% mutate(cells_per_cluster = n()) %>% ungroup() %>% select(CB, seurat_clusters, cells_per_cluster, pc_CGC_Arg), by = "CB") %>%
  group_by(gene_id, gene_name, codon, seurat_clusters, cells_per_cluster, pc_CGC_Arg) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(gene_name, seurat_clusters, pc_CGC_Arg) %>%
  mutate(cluster_gene_frac = n/sum(n),
         cluster_gene_total = sum(n)) %>%
  ungroup() %>% 
  mutate(cluster_gene_mean = cluster_gene_total / cells_per_cluster) %>%
  inner_join(gene_codon_background %>% select(gene_id, codon, codon_frac), by = c("gene_id" = "gene_id", "codon" = "codon")) %>%
  mutate(lfc_change = log2(cluster_gene_frac/codon_frac))
shared_genes <- pull(cluster_change %>%
  group_by(gene_id, seurat_clusters, pc_CGC_Arg) %>%
  summarize(asites_per_cluster = sum(n)) %>%
  ungroup() %>%
  complete(gene_id, seurat_clusters, pc_CGC_Arg) %>%
  mutate(asites_per_cluster = replace_na(asites_per_cluster,0)) %>%
  inner_join(hekMeta %>% group_by(seurat_clusters, pc_CGC_Arg) %>% mutate(cells_per_cluster = n()) %>% ungroup() %>% select(seurat_clusters, cells_per_cluster, pc_CGC_Arg) %>% distinct(), by = c("seurat_clusters", "pc_CGC_Arg")) %>%
  mutate(mean_asites = asites_per_cluster / cells_per_cluster) %>%
  filter(mean_asites > 2.5) %>%
  group_by(gene_id) %>%
  summarize(detected_in_clusters = n()) %>%
  ungroup() %>%
  filter(detected_in_clusters > 2) %>%
  select(gene_id) %>%
  distinct())
highlight_genes <- c("GAPDH", "H3C2", "PRKDC")
gene_codon_change_data <- cluster_change %>%
  mutate(gene_class = case_when(gene_name %in% histone_genes$Approved.symbol ~ "histone",
                                gene_name %in% reactome_cc_genes ~ "cell_cycle",
                             TRUE ~ "rest")) %>%
  filter(codon == highlight_codon) %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name %in% highlight_genes ~ gene_name,
                                  TRUE ~ "")) %>%
  arrange(desc(gene_class))
plt_cluster_change_gaa_facet <- ggplot(gene_codon_change_data, aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = gene_class, label = gene_display ))+
  geom_point(size=0.5, alpha = 1)+
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~seurat_clusters+pc_CGC_Arg, nrow=2)+
  geom_text(data = gene_codon_change_data %>% select(seurat_clusters, cells_per_cluster, pc_CGC_Arg) %>% distinct() %>% mutate(gene_class = "rest"), 
            mapping = aes(x = -Inf, y=-Inf, label = cells_per_cluster),
            hjust = -0.1,
            vjust = -1)+
  geom_text_repel(segment.color = "black", box.padding = 0.3, min.segment.length = 0, max.iter = 10000, seed=112, size=3, point.padding = 0.3)+
  xlab("log10 mean cluster expression")+ylab(paste0("log2 f.c. ", highlight_codon, " P-site FOC"))+
  theme(strip.background = element_blank())
save_plot(paste0(figurePrefix, "_plt_cluster_change_", highlight_codon, ".pdf"), plt_cluster_change_gaa_facet, base_width = 8, base_height =5)


highlight_codon <- "TTA" # this is used for the plot, groupings are based on pc_CODON_Three groups
cluster_change <- codon_reads %>%
  filter(site == "p") %>%
  filter(gene_id %in% pull(gene_codon_background %>% filter(gene_mean > 1) %>% select(gene_id) %>% distinct())) %>%
  inner_join(hekMeta %>% group_by(seurat_clusters, pc_UUA_Leu) %>% mutate(cells_per_cluster = n()) %>% ungroup() %>% select(CB, seurat_clusters, cells_per_cluster, pc_UUA_Leu), by = "CB") %>%
  group_by(gene_id, gene_name, codon, seurat_clusters, cells_per_cluster, pc_UUA_Leu) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(gene_name, seurat_clusters, pc_UUA_Leu) %>%
  mutate(cluster_gene_frac = n/sum(n),
         cluster_gene_total = sum(n)) %>%
  ungroup() %>% 
  mutate(cluster_gene_mean = cluster_gene_total / cells_per_cluster) %>%
  inner_join(gene_codon_background %>% select(gene_id, codon, codon_frac), by = c("gene_id" = "gene_id", "codon" = "codon")) %>%
  mutate(lfc_change = log2(cluster_gene_frac/codon_frac))
shared_genes <- pull(cluster_change %>%
  group_by(gene_id, seurat_clusters, pc_UUA_Leu) %>%
  summarize(asites_per_cluster = sum(n)) %>%
  ungroup() %>%
  complete(gene_id, seurat_clusters, pc_UUA_Leu) %>%
  mutate(asites_per_cluster = replace_na(asites_per_cluster,0)) %>%
  inner_join(hekMeta %>% group_by(seurat_clusters, pc_UUA_Leu) %>% mutate(cells_per_cluster = n()) %>% ungroup() %>% select(seurat_clusters, cells_per_cluster, pc_UUA_Leu) %>% distinct(), by = c("seurat_clusters", "pc_UUA_Leu")) %>%
  mutate(mean_asites = asites_per_cluster / cells_per_cluster) %>%
  filter(mean_asites > 2.5) %>%
  group_by(gene_id) %>%
  summarize(detected_in_clusters = n()) %>%
  ungroup() %>%
  filter(detected_in_clusters > 2) %>%
  select(gene_id) %>%
  distinct())
highlight_genes <- c("GAPDH", "H3C2", "PRKDC")
gene_codon_change_data <- cluster_change %>%
  mutate(gene_class = case_when(gene_name %in% histone_genes$Approved.symbol ~ "histone",
                                gene_name %in% reactome_cc_genes ~ "cell_cycle",
                             TRUE ~ "rest")) %>%
  filter(codon == highlight_codon) %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name %in% highlight_genes ~ gene_name,
                                  TRUE ~ "")) %>%
  arrange(desc(gene_class))
plt_cluster_change_gaa_facet <- ggplot(gene_codon_change_data, aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = gene_class, label = gene_display ))+
  geom_point(size=0.5, alpha = 1)+
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~seurat_clusters+pc_UUA_Leu, nrow=2)+
  geom_text(data = gene_codon_change_data %>% select(seurat_clusters, cells_per_cluster, pc_UUA_Leu) %>% distinct() %>% mutate(gene_class = "rest"), 
            mapping = aes(x = -Inf, y=-Inf, label = cells_per_cluster),
            hjust = -0.1,
            vjust = -1)+
  geom_text_repel(segment.color = "black", box.padding = 0.3, min.segment.length = 0, max.iter = 10000, seed=112, size=3, point.padding = 0.3)+
  xlab("log10 mean cluster expression")+ylab(paste0("log2 f.c. ", highlight_codon, " P-site FOC"))+
  theme(strip.background = element_blank())
save_plot(paste0(figurePrefix, "_plt_cluster_change_", highlight_codon, ".pdf"), plt_cluster_change_gaa_facet, base_width = 8, base_height =5)



### is there a treatment relationship in the clusters?
library(broom)
treatment_contingency <- hekMeta %>%
  count(sort_population, seurat_clusters) %>%
  spread(seurat_clusters, n) %>%
  column_to_rownames('sort_population') 
fisher.test(treatment_contingency, simulate.p.value=T, B=1e6 ) 


CGC_fisher <- hekMeta %>%
  filter(grepl("^Arg", sort_population)) %>%
  count(seurat_clusters, pc_CGC_Arg) %>%
  spread(seurat_clusters, n) %>%
  column_to_rownames('pc_CGC_Arg') %>%
  fisher.test() %>%
  glance()

CGU_fisher <- hekMeta %>%
  filter(grepl("^Arg", sort_population)) %>%
  count(seurat_clusters, pc_CGU_Arg) %>%
  spread(seurat_clusters, n) %>%
  column_to_rownames('pc_CGU_Arg') %>%
  fisher.test() %>%
  glance()

UUA_fisher <- hekMeta %>%
  filter(grepl("^Leu", sort_population)) %>%
  count(seurat_clusters, pc_UUA_Leu) %>%
  spread(seurat_clusters, n) %>%
  column_to_rownames('pc_UUA_Leu')  %>%
  fisher.test() 
CGC_fisher
CGU_fisher
UUA_fisher

plt_pausing_proportion <- bind_rows(hekMeta %>%
                                      filter(grepl("^Arg", sort_population)) %>%
                                      count(seurat_clusters, pc_CGC_Arg) %>%
                                      rename(pausing = pc_CGC_Arg) %>%
                                      mutate(codon = "CGC"),
                                    hekMeta %>%
                                      filter(grepl("^Arg", sort_population)) %>%
                                      count(seurat_clusters, pc_CGU_Arg) %>%
                                      rename(pausing = pc_CGU_Arg) %>%
                                      mutate(codon = "CGU"),
                                    hekMeta %>%
                                      filter(grepl("^Leu", sort_population)) %>%
                                      count(seurat_clusters, pc_UUA_Leu) %>%
                                      rename(pausing = pc_UUA_Leu) %>%
                                      mutate(codon = "UUA")) %>%
  ggplot(aes(x=seurat_clusters, y=n, fill = pausing))+
  geom_bar(stat = "identity", position = "fill")+
  scale_fill_manual(values = brewer.pal(3, "Greys")[2:3])+
  facet_wrap(~codon, nrow=1)
save_plot(paste0(figurePrefix, "_plt_pausing_cluster_proportions.pdf"), plt_pausing_proportion, base_width = 6, base_height = 4)

pausing_proportion_numbers <- bind_rows(hekMeta %>%
                                      filter(grepl("^Arg", sort_population)) %>%
                                      count(seurat_clusters, pc_CGC_Arg) %>%
                                      rename(pausing = pc_CGC_Arg) %>%
                                      mutate(codon = "CGC"),
                                    hekMeta %>%
                                      filter(grepl("^Arg", sort_population)) %>%
                                      count(seurat_clusters, pc_CGU_Arg) %>%
                                      rename(pausing = pc_CGU_Arg) %>%
                                      mutate(codon = "CGU"),
                                    hekMeta %>%
                                      filter(grepl("^Leu", sort_population)) %>%
                                      count(seurat_clusters, pc_UUA_Leu) %>%
                                      rename(pausing = pc_UUA_Leu) %>%
                                      mutate(codon = "UUA")) 


arg_pausing_cells <- hekMeta %>%
  mutate(pause_arg = (pc_CGC_Arg | pc_CGU_Arg),
         pause_leu = pc_UUA_Leu) %>%
  group_by(seurat_clusters, pause_arg) %>%
  summarize(n = n()) %>%
  ungroup()

featplot <- function(U,aa){
  # draws outline at 98% density around clusters
  plt <- U %>% arrange(.data[[aa]]) %>% ggplot(aes(x=UMAP_1, y=UMAP_2, fill = log2(.data[[aa]]))) +
    geom_point(size=2, shape = 21, colour = "grey60", stroke = 0.1) +
    scale_fill_gradientn(colours=rev(brewer.pal(11,"RdBu")),na.value="grey40", limits = c(-2,2), oob=squish)+
         geom_density_2d(aes(x=UMAP_1, y=UMAP_2, colour = seurat_clusters), data = U, binwidth = 0.005, alpha = 0.8, linetype = "dashed")+
    coord_equal()+
    theme(legend.title = element_blank(),
          axis.line = element_blank(), 
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())+
    ggtitle(aa)+
    expand_limits(x=c(min(U$UMAP_1)-1, max(U$UMAP_1)+1),y=c(min(U$UMAP_2)-1, max(U$UMAP_2)+1))

  gg <- ggplot_build(plt)
  gg$data[[2]] <- subset(gg$data[[2]], level == 0.02)
  p1 <- ggplot_gtable(gg)
  return(p1)
}

plt_CGC_Arg <- featplot(hekMeta, "pc_CGC_Arg")
plt_CGU_Arg <- featplot(hekMeta, "pc_CGU_Arg")
plt_UUA_Leu <- featplot(hekMeta, "pc_UUA_Leu")
save_plot(paste0(figurePrefix, "_plt_aas_umap_t.pdf"),
                 plot_grid(plt_CGC_Arg,
                           plt_CGU_Arg,
                           plt_UUA_Leu, nrow=1),
                 base_width = 16,
                 base_height = 16/3)


display_gene <- "PRKDC"

H3C2codons <- codons %>%
  filter(transcript_id == pull(annot %>% filter(gene_name == display_gene, set == "canonical") %>% select(transcript_id)))

H3C2 <- H3C2codons %>%
  inner_join(codon_reads %>% filter(site %in% c("p")) %>% filter(CB %in% pull(hekMeta %>% filter(seurat_clusters %in% c(0)) %>% select(CB))), by = c("transcript_id" = "transcript_id", "codon_position" = "site_position")) %>%
  group_by(codon_position,CB) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  complete(codon_position = full_seq(c(min(H3C2codons$codon_position), max(H3C2codons$codon_position)),3), CB) %>%
  mutate(n = replace_na(n,0)) %>%
  group_by(CB) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  inner_join(H3C2codons, by="codon_position") %>%
  inner_join(hekMeta, by="CB") %>%
  mutate(cell_frac = n/cds)

H3C2_hm <- H3C2 %>%
  select(CB, codon_position, cell_frac) %>%
  spread(codon_position, cell_frac, fill = 0) %>%
  as.data.frame()
rownames(H3C2_hm) <- H3C2_hm$CB
H3C2_hm <- H3C2_hm[,colnames(H3C2_hm) != "CB"]
H3C2_hm <- H3C2_hm[pull( hekMeta %>% filter(CB %in% rownames(H3C2_hm)) %>% arrange(-ps_CGC_Arg) %>% select(CB)),]
H3C2_hm <- H3C2_hm[!grepl("^NA",rownames(H3C2_hm)),]

argSortSplit <- data.frame(CB = pull(hekMeta %>%
                                     filter(CB %in% rownames(H3C2_hm)) %>%
                                     arrange(-ps_CGC_Arg) %>%
                                     select(CB)),
                           sortSplit = c(rep("0",25),
                                         rep("1", 25),
                                         rep("2", 25),
                                         rep("3", 25),
                                         rep("rest", nrow(H3C2_hm)-100)),
                           stringsAsFactors = FALSE)
rownames(argSortSplit) <- argSortSplit$CB

rownames(hekMeta) <- hekMeta$CB

aa_highlight <- rep("other", nrow(H3C2codons))
aa_highlight[grep("CGC", H3C2codons$codon)] <- "CGC"
aa_highlight[grep("CGT", H3C2codons$codon)] <- "CGT"
aa_highlight[grep("TTA", H3C2codons$codon)] <- "TTA"


pdf(paste0(figurePrefix, "_pausing_hm_", display_gene, ".pdf"), height = 3, width = 8)
aa_show <- HeatmapAnnotation(aa = aa_highlight,
                             col = list(aa = c("other" = "white",
                                               "CGT" = "#377eb8",
                                               "TTA" = "#984ea3",
                                               "CGC" = "#e41a1c" )),
                             annotation_name_side = "left")
aa_annotation <- HeatmapAnnotation(aa = anno_text(H3C2codons$single, gp = gpar(fontsize = 4)))
cell_annotation <- rowAnnotation(treatment = as.character(hekMeta[rownames(H3C2_hm), "sort_population"]),
                                 mean_CGC = log2(hekMeta[rownames(H3C2_hm), "ps_CGC_Arg"]),
                                 col = list(seurat_clusters = c("0" = gg_color_hue(4)[1],
                                                                "1" = gg_color_hue(4)[2],
                                                                "2" = gg_color_hue(4)[3],
                                                                "3" = gg_color_hue(4)[4]),
                                            treatment = c("Arg3h" = as.character(coltreat[1]),
                                                          "Arg6h" = as.character(coltreat[2]),
                                                          "Leu3h" = as.character(coltreat[3]),
                                                          "Leu6h" = as.character(coltreat[4]),
                                                          "Rich" =  as.character(coltreat[5])),
                                            mean_CGC = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))), # log fold change
                                            mean_leu = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))) ) )
length_annotation <- HeatmapAnnotation(pos = anno_mark(at = c(seq(25,ncol(H3C2_hm), by = 25)), labels = as.character(c(seq(25,ncol(H3C2_hm), by = 25))), which = "column", side = "bottom", labels_rot = 0, link_height = unit(3, "mm")))
hm_quantiles <- 0.001
Heatmap(as.matrix(H3C2_hm),
        top_annotation = aa_show,
        right_annotation = cell_annotation,
        bottom_annotation = length_annotation,
        col = colorRamp2(seq(quantile(as.matrix(H3C2_hm), probs = hm_quantiles), quantile(as.matrix(H3C2_hm), probs = 1-hm_quantiles), length = 9), brewer.pal(9, "Greys")),
        na_col = "grey0",
        name = "fraction of reads/cell",
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        clustering_distance_rows = function(x) as.dist(1-cosine_similarity(x)),
        clustering_method_rows = "ward.D2",
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_title = NULL,
        column_title_rot = 90,
        border = TRUE,
        use_raster = FALSE)
dev.off()





### Is the clustering driven by the pausing reads

## cluster based on UTR reads  -- doesn't work
#start_expansion = 25
#stop_expansion = 25
#site = "psite"
#
#read_classification <- hekReads %>%
#  select(-starts_with('base')) %>%
#  mutate(lr_utr5 = l_utr5 - start_expansion,
#         lr_cds = l_cds + start_expansion + stop_expansion,
#         lr_utr3 = l_utr3 - stop_expansion) %>%
#  mutate(lr_utr5 = case_when(lr_utr5 < 0 ~ 0,
#                             lr_utr5 >= 0 ~ as.numeric(lr_utr5)),
#         lr_cds = case_when(lr_cds < 0 ~ 0,
#                            lr_cds >= 0 ~ as.numeric(lr_cds)),
#         lr_utr3 = case_when(lr_utr3 < 0 ~ 0,
#                             lr_utr3 >= 0 ~ as.numeric(lr_utr3))) %>%
#  mutate(read_region = case_when(psite <= lr_utr5 ~ 'utr5',
#                                 ( (psite > (lr_utr5)) & (psite <= (lr_utr5 + lr_cds)) ) ~ 'cds',
#                                 ( (psite > (lr_utr5 + lr_cds)) & (psite <= (lr_utr5 + lr_cds + lr_utr3)) ) ~ 'utr3',
#                                 TRUE ~ 'NA')) %>%
#  mutate(count_id = paste(gene_id, gene_name, read_region, sep = '--')) %>%
#  group_by(CB, count_id) %>%
#  summarize(n=n()) %>%
#  ungroup() 
#
#utrCounts <- read_classification %>%
#  separate(count_id, into = c("gene_id", "gene_name", "read_region"), sep = "--", remove = FALSE) %>%
#  filter(read_region != "cds") %>%
#  mutate(count_id = paste(gene_id, gene_name, sep = "--")) %>%
#  group_by(CB, count_id) %>%
#  summarize(counts = sum(n)) %>%
#  ungroup() %>%
#  rename(n = counts)



## how many reads contain an arginine # expand reads to one line per covered base
expandedReads <- hekReads %>%
  select(-starts_with('base')) %>%
  rowwise() %>%
  mutate(range = list(seq(cut5, cut3))) %>%
  unnest(cols = range) %>%
  inner_join(codons %>% select(-single), by = c("transcript_id" = "transcript_id", "range" = "codon_position")) %>%
  group_by(CB, transcript_id, id) %>%
  #summarize(num_Arg = sum(codon %in% c("CGC", "CGT")),
  #          num_Leu = sum(codon %in% c("TTA"))) %>%
  summarize(num_Arg = sum(three == "Arg"),
            num_Leu = sum(three == "Leu")) %>%
  ungroup() %>%
  mutate(contains_Arg = num_Arg > 0,
         contains_Leu = num_Leu > 0)

plt_arg_fraction <- expandedReads %>%
  group_by(CB) %>%
  summarize(n_reads = n(),
            n_arg = sum(contains_Arg)) %>%
  ungroup() %>%
  inner_join(hekMeta %>% select(CB, sort_population, starts_with('pc_'), seurat_clusters), by = "CB") %>%
  ggplot(aes(x=seurat_clusters, y=n_arg/n_reads, colour = sort_population))+
  geom_point()

save_plot(paste0(figurePrefix, "_plt_arg_fraction.pdf"), plt_arg_fraction, base_width = 5, base_height = 5)


# remove Arg reads and count
noArgCounts <- expandedReads %>%
  filter(!contains_Arg) %>%
  filter(!contains_Leu) %>%
  inner_join(canonical_pc, by = "transcript_id") %>%
  group_by(CB, gene_id) %>%
  summarize(n = n()) %>%
  ungroup()

noArgLeu <- with(noArgCounts %>%
                 mutate_if(is.character, as.factor),
               sparseMatrix(i = as.numeric(gene_id),
                            j = as.numeric(CB),
                            x = n,
                            dimnames = list(levels(gene_id),
                                            levels(CB))))




#geneRegionCounts <- with(read_classification %>% 
#                         mutate_if(is.factor, as.character) %>%
#                         filter(grepl("--utr(5|3)$", count_id)) %>%
#                         mutate_if(is.character, as.factor), 
#                       sparseMatrix(i = as.numeric(count_id),
#                                    j = as.numeric(CB),
#                                    x = n,
#                                    dimnames = list(levels(count_id),
#                                                    levels(CB))))
#
#geneRegionCounts <- with(utrCounts %>% 
#                         mutate_if(is.character, as.factor), 
#                       sparseMatrix(i = as.numeric(count_id),
#                                    j = as.numeric(CB),
#                                    x = n,
#                                    dimnames = list(levels(count_id),
#                                                    levels(CB))))


PR <- CreateSeuratObject(counts = noArgLeu,
                        project = "RPF_starv",
                        assay = "RNA",
                        min.cells = 0,
                        min.features = 0,
                        names.field = 1,
                        names.delim = "-",
                        meta.data = hekMeta)
PR <- NormalizeData(PR)

PR <- FindVariableFeatures(PR, selection.method = "vst", nfeatures = 1000)

# plot variable features with and without labels
plt_variable <- VariableFeaturePlot(PR)
plt_variable_labelled <- LabelPoints(plot = plt_variable, points = head(VariableFeatures(PR),10), repel = TRUE)

PR <- ScaleData(PR, features = rownames(PR))
PR <- RunPCA(PR, features = VariableFeatures(object = PR))

PR <- JackStraw(PR, num.replicate = 100)
PR <- ScoreJackStraw(PR, dims = 1:20)

plt_jackstraw <- JackStrawPlot(PR,dims=1:20)
plt_elbow <- ElbowPlot(PR)

save_plot(paste0(figurePrefix,"_qc_plt_dimscores_region.pdf"), plot_grid(plt_jackstraw, plt_elbow, nrow=1), base_width = 10, base_height = 5)

ndims <- 3
PR <- RunUMAP(PR, dims = 1:ndims, verbose = FALSE, seed.use = 112)#, n.epochs = 3000)
PR <- FindNeighbors(PR, dims = 1:ndims)
PR <- FindClusters(PR, resolution = 0.4) # remove all Arg, Leu codon-spanning reads

plt_reg_umap_plate <- DimPlot(PR, reduction = "umap",group.by = "plate",pt.size=1)
plt_reg_umap_cluster <- DimPlot(PR, reduction = "umap",pt.size=1)
plt_reg_umap_treatment <- DimPlot(PR, reduction = 'umap', group.by = "sort_population", pt.size=1)


save_plot(paste0(figurePrefix,"_plt_reg_umaps_check_region.pdf"), plot_grid(plt_reg_umap_plate+coord_equal(), plt_reg_umap_treatment+coord_equal() ,plt_reg_umap_cluster+coord_equal(), nrow=2), base_width = 10, base_height = 10)

hekMeta$noarg_seurat_clusters <- PR$seurat_clusters[rownames(hekMeta)]


plt_umap_all <- ggplot(hekMeta, aes(x=UMAP_1, y=UMAP_2, colour = seurat_clusters))+
  geom_point(size=0.5)+
  coord_equal()
plt_umap_noarg <- ggplot(hekMeta, aes(x=UMAP_1, y=UMAP_2, colour = noarg_seurat_clusters))+
  geom_point(size=0.5)+
  coord_equal()
save_plot(paste0(figurePrefix, "_plt_noarg_umap_comparison.pdf"), plot_grid(plt_umap_all, plt_umap_noarg, nrow=1, align = 'hv'), base_width = 10, base_height = 4)


getMode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

cluster_translation <- hekMeta %>%
  group_by(noarg_seurat_clusters) %>%
  summarize(translated_cluster = getMode(seurat_clusters)) %>%
  ungroup()

## what is the overlap between clusters, how many kept the same identity?
clusterConcordance <- hekMeta %>%
  inner_join(cluster_translation, by="noarg_seurat_clusters") %>%
  summarize(num_same = sum(seurat_clusters == translated_cluster),
            total = n()) %>%
  mutate(fraction = num_same/total)


readChanges <- expandedReads %>%
  inner_join(canonical_pc) %>%
  summarize(n = n(),
            numArg = sum(contains_Arg),
            numLeu = sum(contains_Leu),
            numBoth = sum(contains_Arg | contains_Leu)) %>%
  mutate(propA = numArg/n,
         propL = numLeu/n,
         propB = numBoth/n)



cluster_proportions <- hekMeta %>%
  group_by(seurat_clusters, sort_population) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters, y=n, fill = sort_population))+
  geom_bar(position = "fill", stat="identity")+
  scale_fill_manual(values = as.character(coltreat))
save_plot(paste0(figurePrefix, "_plt_cluster_proportions.pdf"), cluster_proportions, base_width = 4, base_height = 4)

cluster_proportion_numbers <- hekMeta %>%
  group_by(seurat_clusters, sort_population) %>%
  summarize(n = n()) %>%
  ungroup() 


## Pausing heatmaps separated by treatment

codons_heatmap <- function(tbl, codons, treatments, feature){
  out_colnames <- c(t(outer(codons, levels(tbl$site), FUN = "paste", sep = "_")))
  hm <- tbl %>%
    filter(sort_population %in% treatments,
           codon_display %in% codons) %>%
    mutate(site_out = paste(codon_display, site, sep = "_")) %>%
    select(CB, site_out, {{ feature }}) %>%
    spread(site_out, {{ feature }}, fill = 0) %>%
    as.data.frame()
  rownames(hm) <- hm$CB
  out_colnames <- out_colnames[out_colnames %in% colnames(hm)]
  hm <- hm[,out_colnames]
  d <- as.dist(1-cosine_similarity(as.matrix(hm)))
  hc <- hclust(d, method = "ward.D2")
  hm <- hm[hc$order,]
  return(hm)
}

display_codons <- c("CGC_Arg", "CGU_Arg", "UUA_Leu")

rich_selected <- site_change %>%
  codons_heatmap(codons = display_codons,  
                 treatments = "Rich",
                 feature = site_fc)
leu3h_selected <- site_change %>%
  codons_heatmap(codons = display_codons, 
                 treatments = "Leu3h",
                 feature = site_fc)
leu6h_selected <- site_change %>%
  codons_heatmap(codons = display_codons, 
                 treatments = "Leu6h",
                 feature = site_fc)
arg3h_selected <- site_change %>%
  codons_heatmap(codons = display_codons, 
                 treatments = "Arg3h",
                 feature = site_fc)
arg6h_selected <- site_change %>%
  codons_heatmap(codons = display_codons, 
                 treatments = "Arg6h",
                 feature = site_fc)

selected_hm <- rbind(rich_selected,
                     arg3h_selected,
                     arg6h_selected,
                     leu3h_selected,
                     leu6h_selected)

selected_hm <- log2(as.matrix(selected_hm))
selected_hm[is.infinite(selected_hm)] <- NA

pdf(paste0(figurePrefix, "_selected_codon_stalling_site_l2fc.pdf"), height = 3, width = 8)
site_labels <- HeatmapAnnotation(site = anno_text(rep(c("-2","-1","E","P","A","+1","+2"), times = 3)),
                                 which = "row")
cell_annotation <- HeatmapAnnotation(
                                     clusters = as.character(hekMeta[rownames(selected_hm), "seurat_clusters"]),
                                     col = list(clusters = c("0" = gg_color_hue(4)[1],
                                                             "1" = gg_color_hue(4)[2],
                                                             "2" = gg_color_hue(4)[3],
                                                             "3" = gg_color_hue(4)[4])),
                                     which = "column",
                                     annotation_name_side = "right")
Heatmap(t(selected_hm),
        col = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))), # log fold change
        top_annotation = cell_annotation,
        right_annotation = site_labels,
        na_col = "grey0",
        name = "log2 fold change",
        cluster_columns = TRUE,
        show_column_dend = FALSE,
        clustering_method_columns = "complete",
        cluster_rows = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = factor(hekMeta[rownames(selected_hm), "sort_population"], levels = c("Rich","Arg3h", "Arg6h", "Leu3h", "Leu6h"), ordered = TRUE),
        row_split = factor(rep(display_codons, each = length(levels(site_change$site))), levels = display_codons, ordered = TRUE),
        row_title_rot = 0,
        border = TRUE,
        use_raster = TRUE,
        raster_device = "CairoPNG",
        raster_quality = 15)
dev.off()

