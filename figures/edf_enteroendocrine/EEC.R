#!/usr/bin/env Rscript
# Michael VanInsberghe 2020-05-06
# EE Pausing on EEC data

library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(broom)
library(data.table)
library(feather)
library(Matrix)
library(Matrix.utils)
library(Seurat)
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
library(dendsort)
theme_set(theme_cowplot())

source('plotFunctions.R')

figurePrefix <- "ee_raw/EEpausing"

data_dir <- "../scRiboSeq/data"

# barcodes
barcodes <- read.table('/hpc/hub_oudenaarden/mvanins/local/lib/barcodes/miRv6_barcodemap.tsv',
                       sep="\t",
                       header = FALSE,
                       stringsAsFactors = FALSE,
                       row.names = NULL,
                       col.names = c("CB","number","well"))
rownames(barcodes) <- barcodes$CB

# annot
annot <- fread(file.path(data_dir, 'mmu_annotations.csv.gz'), stringsAsFactors = FALSE, data.table = FALSE)

canonical_pc <- annot %>%
  filter(set=="canonical", transcript_type == "protein_coding") %>%
  filter(chr != 'chrM') %>%
  filter(transcript_id != "ENST00000437139.7") %>%
  filter(!grepl("^Pcdh", gene_name)) #%>%
  #filter(!grepl("^RPL", gene_name))

# codons
codons <- fread(file.path(data_dir, 'mmu_codons.csv.gz'), stringsAsFactors = FALSE, data.table = FALSE) %>%
  filter( transcript_id %in% canonical_pc$transcript_id)

RPFreads <- readRDS(file.path(data_dir, "EE/compiled/reads.rds"))
RPFmeta <- readRDS(file.path(data_dir, "EE/compiled/meta.rds"))
RPFcounts <- readRDS(file.path(data_dir, "EE/compiled/counts.rds"))
RPFbiotypes <- readRDS(file.path(data_dir, "EE/compiled/biotypes.rds"))
RPFallMeta <- readRDS(file.path(data_dir, "EE/compiled/allmeta.rds"))

RPFlong <- as.data.frame(t(RPFcounts))
RPFlong$CB <- rownames(RPFlong)
RPFlong <- RPFlong %>%
  pivot_longer(cols = starts_with("ENSMUSG"), names_to = "gene_id", values_to = "counts") 




ENSG_to_name <- data.frame(count_names = rownames(RPFcounts), stringsAsFactors = FALSE) %>%
  inner_join(annot, by = c("count_names" = "gene_id")) %>% 
  mutate(gene_out = paste(count_names, gene_name, sep = "-")) %>%
  select(count_names, gene_out) %>%
  distinct() %>%
  column_to_rownames(var = "count_names")

rownames(RPFcounts) <- ENSG_to_name[rownames(RPFcounts),"gene_out"]
RPFcounts <- RPFcounts[!grepl("-Pcdh", rownames(RPFcounts)),]

RPFmeta$log_green <- log10(RPFmeta$green)
RPFmeta$log_red <- log10(RPFmeta$red)


P <- CreateSeuratObject(counts = RPFcounts,
                        project = "RPF_EE",
                        assay = "RNA",
                        min.cells = 10,
                        min.features = 10,
                        names.field = 1,
                        names.delim = "-",
                        meta.data = RPFmeta)
P[["percent_mt"]] <- PercentageFeatureSet(P, pattern = "-mt-")

P <- NormalizeData(P, verbose = FALSE)
P <- FindVariableFeatures(P, selection.method = "vst", nfeatures = 2000)

P <- ScaleData(object = P, verbose=FALSE)
P <- RunPCA(object=P, verbose=FALSE)
P <- JackStraw(P, num.replicate = 100, dims = 50)
P <- ScoreJackStraw(P, dims = 1:50)

# plt_jackstraw <- JackStrawPlot(P,dims=1:50)
# plt_elbow <- ElbowPlot(P, ndims=50)
# save_plot(paste0(figurePrefix,"_qc_plt_dimscores.pdf"), plot_grid(plt_elbow, plt_jackstraw, nrow=1), base_width = 15, base_height = 5)

ndims <- 14
P <- RunUMAP(P, dims=1:ndims, verbose=FALSE, n.neighbors = 30L, n.epochs = 500, return.model = TRUE, seed.use=112)
P <- FindNeighbors(P, dims=1:ndims, k.param=25 )
P <- FindClusters(P,  n.start=500, n.iter=500, resolution = 1.5, algorithm = 2)

RPFmeta$seurat_clusters <- P$seurat_clusters[rownames(RPFmeta)]
RPFmeta$percent_mt <- P$percent_mt[rownames(RPFmeta)]
RPFmeta$UMAP_1 <- as.data.frame(Embeddings(P, reduction="umap"))[rownames(RPFmeta),"UMAP_1"]
RPFmeta$UMAP_2 <- as.data.frame(Embeddings(P, reduction="umap"))[rownames(RPFmeta),"UMAP_2"]

redefineClusters <- data.frame(seurat_clusters = levels(RPFmeta$seurat_clusters),
                               output_clusters = c(1, 6, 2, 8, 5, 4, 7, 3))
rownames(redefineClusters) <- redefineClusters$seurat_clusters

RPFmeta$output_clusters <- as.factor(redefineClusters[RPFmeta$seurat_clusters, "output_clusters"])

saveRDS(RPFmeta, file.path(data_dir, "EE/compiled/clustered_meta.rds"))

## dimensionality reduction plots

set.seed(112)
plt_facs_clusters <- ggplot(RPFmeta, aes(x=log_green, y=log_red, colour = output_clusters))+
  geom_point(data = RPFmeta %>% select(log_green, log_red), colour = "grey70", size = 0.8)+
  geom_point(size = 0.8)+
  facet_wrap(~output_clusters, ncol = 2)+
  xlab("log10 mNeonGreen")+
  ylab("log10 dTomato")+
  coord_equal()+
  theme(legend.position = "none",
        strip.background = element_blank())
save_plot(paste0(figurePrefix, "_plt_facs_facet.pdf"), plt_facs_clusters, base_width = 6, base_height = 6.1, limitsize=FALSE)

plt_umap_clusters <- ggplot(RPFmeta[sample(nrow(RPFmeta), replace=FALSE),], aes(x=UMAP_1, y=UMAP_2, colour = output_clusters))+
  geom_point()+
  coord_equal() +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_region<- ggplot(RPFmeta %>% arrange(desc(sort_population)), aes(x=UMAP_1, y=UMAP_2, colour = sort_population))+
  geom_point()+
  scale_colour_manual(values = rev(brewer.pal(3, "Dark2")))+
  coord_equal()+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_red <- ggplot(RPFmeta %>% arrange(log_red), aes(x=UMAP_1, y=UMAP_2, colour = log_red))+
  geom_point()+
  coord_equal()+
  scale_colour_gradientn(colours = brewer.pal(8,"Reds"), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_green <- ggplot(RPFmeta %>% arrange(log_green), aes(x=UMAP_1, y=UMAP_2, colour = log_green))+
  geom_point()+
  coord_equal()+
  scale_colour_gradientn(colours = brewer.pal(8,"Greens"), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
save_plot(paste0(figurePrefix, "_plt_reduction_panels.pdf"), plot_grid(plt_umap_clusters,plt_umap_region,
                                                                       plt_umap_green, plt_umap_red, 
                                                                       nrow = 2, align = "hv"), base_width = 12, base_height = 9)

save_plot(paste0(figurePrefix, "_plt_UMAP_panels.pdf"),
          plot_grid(plt_umap_clusters+theme(legend.position = "bottom"), 
                    plot_grid(plt_umap_region, plt_umap_green, plt_umap_red, ncol = 1),
                    nrow = 1),
          base_width = 6, base_height = 5.15)



## marker genes 
markers <- FindAllMarkers(P, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


plt_first_markers <- FeaturePlot(P, features = c(rownames(RPFcounts)[c( 
                                                     grep("-Dll1$", rownames(RPFcounts)),
                                                     grep("-Neurog3$", rownames(RPFcounts)),
                                                     grep("-Tac1$", rownames(RPFcounts)),
                                                     grep("-Tph1$", rownames(RPFcounts)),
                                                     grep("-Chga$", rownames(RPFcounts)),
                                                     grep("-Chgb$", rownames(RPFcounts)),
                                                     grep("-Reg4$", rownames(RPFcounts)),
                                                     grep("-Isl1", rownames(RPFcounts)),
                                                     grep("-Ghrl$", rownames(RPFcounts)),
                                                     grep("-Cck$", rownames(RPFcounts)),
                                                     grep("-Gcg$", rownames(RPFcounts)),
                                                     grep("-Gip$", rownames(RPFcounts)),
                                                     grep("-Sst$", rownames(RPFcounts)),
                                                     grep("-Agr2$", rownames(RPFcounts)),
                                                     grep("-Spink4$", rownames(RPFcounts))
                                                                       )]), order = TRUE, combine = FALSE)
save_plot(paste0(figurePrefix,"_plt_umaps_first_markers.pdf"), plot_grid(plotlist = reformatFeature(plt_first_markers), ncol=3), base_width = 10, base_height=10)

# marker heatmap
topmarkers <- markers %>% 
  filter(p_val_adj < 0.01) %>% 
  filter(!grepl("-AC[0-9][0-9]+", gene)) %>% 
  filter(!grepl("-AL[0-9][0-9]+", gene))

markerHM <- as.data.frame(GetAssayData(P, slot = "data"))[pull(topmarkers %>% ungroup() %>% select(gene) %>% distinct()), pull(RPFmeta %>% arrange(output_clusters) %>% select(CB))]

markerHM <- sweep(markerHM, 1, apply(markerHM,1,mean), FUN = "-")
markerHM <- sweep(markerHM, 1, apply(markerHM,1,sd), FUN = "/")

log2_fc_RPFpC = log2(RPFmeta[colnames(markerHM), "cds"]/mean(RPFmeta[colnames(markerHM), "cds"]))

pdf(paste0(figurePrefix, "_marker_hm_clusters.pdf"), height = 8, width = 8)
cell_annotation <- HeatmapAnnotation(clusters = as.character(RPFmeta[colnames(markerHM), "output_clusters"]),
                                     region = as.character(RPFmeta[colnames(markerHM), "sort_population"]),
                                     mNeonGreen = RPFmeta[colnames(markerHM), "log_green"],
                                     dTomato = RPFmeta[colnames(markerHM), "log_red"],
                                     # log2_fc_RPFpC = log2_fc_RPFpC, 
                                     col = list(clusters = c("1" = gg_color_hue(8)[1],
                                                             "2" = gg_color_hue(8)[2],
                                                             "3" = gg_color_hue(8)[3],
                                                             "4" = gg_color_hue(8)[4],
                                                             "5" = gg_color_hue(8)[5],
                                                             "6" = gg_color_hue(8)[6],
                                                             "7" = gg_color_hue(8)[7],
                                                             "8" = gg_color_hue(8)[8]),
                                                region = c("Distal" = brewer.pal(3, "Dark2")[3],
                                                           "Medial" = brewer.pal(3, "Dark2")[2],
                                                           "Proximal" = brewer.pal(3, "Dark2")[1]),
                                                dTomato = colorRamp2(seq(min(RPFmeta$log_red, na.rm = TRUE), max(RPFmeta$log_red, na.rm = TRUE), length.out = 8), brewer.pal(8, 'Reds')),
                                                mNeonGreen = colorRamp2(seq(min(RPFmeta$log_green, na.rm = TRUE), max(RPFmeta$log_green, na.rm = TRUE), length.out = 8), brewer.pal(8, 'Greens'))),
                                     annotation_name_side = "right")
cc_annotation <- rowAnnotation(cc = anno_mark(at = c(
                                                     grep("-Dll1$", rownames(markerHM)),
                                                     grep("-Neurog3$", rownames(markerHM)),
                                                     grep("-Pigr$", rownames(markerHM)),
                                                     grep("-Pax6$", rownames(markerHM)),
                                                     grep("-Neurod1$", rownames(markerHM)),
                                                     grep("-Chgb$", rownames(markerHM)),
                                                     grep("-Mcm6$", rownames(markerHM)),
                                                     grep("-Mki67$", rownames(markerHM)),
                                                     grep("-Tac1$", rownames(markerHM)),
                                                     grep("-Tph1$", rownames(markerHM)),
                                                     grep("-Chga$", rownames(markerHM)),
                                                     grep("-Reg4$", rownames(markerHM)),
                                                     grep("-Isl1", rownames(markerHM)),
                                                     grep("-Cck$", rownames(markerHM)),
                                                     grep("-Gcg$", rownames(markerHM)),
                                                     grep("-Gip$", rownames(markerHM)),
                                                     grep("-Sst$", rownames(markerHM)),
                                                     grep("-Ghrl$", rownames(markerHM)),
                                                     grep("-Agr2$", rownames(markerHM)),
                                                     grep("-Spink4$", rownames(markerHM))
                                                     ),
                                              labels = c(
                                                     "Dll1",
                                                     "Neurog3",
                                                     "Pigr",
                                                     "Pax6",
                                                     "Neurod1",
                                                     "Chgb",
                                                     "Mcm6",
                                                     "Mki67",
                                                     "Tac1",
                                                     "Tph1",
                                                     "Chga",
                                                     "Reg4",
                                                     "Isl1",
                                                     "Cck",
                                                     "Gcg",
                                                     "Gip",
                                                     "Sst",
                                                     "Ghrl",
                                                     "Agr2",
                                                     "Spink4"
                                                     ) ) )
row_dend <- dendsort(hclust(as.dist(1-cosine_similarity(as.matrix(markerHM))), method = "ward.D2"), type = "min")
Heatmap(as.matrix(markerHM),
        top_annotation = cell_annotation,
        right_annotation = cc_annotation,
        na_col = "grey0",
        name = "z-score expression",
        cluster_rows = row_dend,
        show_row_dend = FALSE,
        clustering_distance_rows = function(x) as.dist(1-cosine_similarity(x)),
        clustering_method_rows = "ward.D2",
        cluster_columns = TRUE,
        clustering_distance_columns = function(x) as.dist(1-cosine_similarity(x)),
        clustering_method_columns = "ward.D2",
        show_column_dend = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = RPFmeta[colnames(markerHM), "output_clusters"],
        column_title_rot = 90,
        border = FALSE,
        use_raster = TRUE,
        raster_device = "CairoPNG",
        raster_quality = 20)
dev.off()


### pausing

codon_reads <- RPFreads %>%
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
  inner_join(RPFmeta %>% select(CB, red, green, output_clusters, sort_population, plate ), by = "CB") %>%
  inner_join(total_baseline, by = "codon_display") %>%
  inner_join(site_baseline, by = c("codon_display" = "codon_display", "site" = "site")) %>%
  mutate(total_fc = cell_total_fraction/(baseline_total_fraction/nrow(site_order)),
         site_fc = cell_site_fraction/baseline_site_fraction) %>%
  inner_join(site_order, by = "site") %>%
  mutate(site = fct_reorder(site, site_order))



filterGenes <- RPFlong %>%
  filter(counts > 0) %>%
  group_by(CB) %>%
  mutate(cell_fraction = counts/sum(counts)) %>%
  ungroup() %>%
  inner_join(RPFmeta %>% select(CB, output_clusters), by = "CB") %>%
  group_by(gene_id, output_clusters) %>%
  summarize(median_cell_fraction = mean(cell_fraction)) %>%
  ungroup() %>%
  filter(median_cell_fraction > 0.025) %>%
  select(gene_id) %>%
  distinct() %>%
  pull()

depleted_codon_reads <- codon_reads %>%
  filter(!(gene_id %in% filterGenes))

depleted_baseline_counts <- depleted_codon_reads %>%
  filter(!(three %in% c("Ter"))) %>%
  mutate(codon_display = paste(gsub("T", "U", codon), three, sep = "_")) %>%
  group_by(CB, codon_display, site) %>%
  summarize(n = n()) %>%
  ungroup()

depleted_total_baseline <- depleted_baseline_counts %>%
  group_by(CB, codon_display) %>%
  summarize(codon_counts = sum(n)) %>%
  ungroup() %>%
  group_by(CB) %>%
  mutate(cell_total_fraction = codon_counts/sum(codon_counts)) %>%
  ungroup() %>%
  group_by(codon_display) %>%
  summarize(baseline_total_fraction = mean(cell_total_fraction)) %>%
  ungroup()

depleted_site_baseline <- depleted_baseline_counts %>%
  group_by(CB, site) %>%
  mutate(cell_site_fraction = n/sum(n)) %>%
  ungroup() %>%
  group_by(codon_display, site) %>%
  summarize(baseline_site_fraction = mean(cell_site_fraction)) %>%
  ungroup()

depleted_site_change <- depleted_codon_reads %>%
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
  inner_join(RPFmeta %>% select(CB, red, green, output_clusters, sort_population, plate ), by = "CB") %>%
  inner_join(depleted_total_baseline, by = "codon_display") %>%
  inner_join(depleted_site_baseline, by = c("codon_display" = "codon_display", "site" = "site")) %>%
  mutate(total_fc = cell_total_fraction/(baseline_total_fraction/nrow(site_order)),
         site_fc = cell_site_fraction/baseline_site_fraction) %>%
  inner_join(site_order, by = "site") %>%
  mutate(site = fct_reorder(site, site_order))





pausing_signal_CAG <- depleted_site_change %>%
  filter(codon_display %in% c("CAG_Gln")) %>%
  filter(site %in% c("a")) %>%
  group_by(CB, codon_display) %>%
  summarize(pausing = mean(site_fc)) %>%
  ungroup() %>%
  inner_join(RPFmeta %>% select(CB, sort_population, output_clusters), by = "CB")
pausing_signal_GAA <- depleted_site_change %>%
  filter(codon_display %in% c("GAA_Glu")) %>%
  filter(site %in% c("a")) %>%
  group_by(CB, codon_display) %>%
  summarize(pausing = mean(site_fc)) %>%
  ungroup() %>%
  inner_join(RPFmeta %>% select(CB, sort_population, output_clusters), by = "CB")

RPFmeta$cag_clusters <- "background"
RPFmeta$gaa_clusters <- "background"
RPFmeta[pausing_signal_CAG %>% filter(pausing > 2^.8) %>% select(CB) %>% pull(), "cag_clusters"] <- "cag_pausing"
RPFmeta[pausing_signal_GAA %>% filter(pausing > 2^.8) %>% select(CB) %>% pull(), "gaa_clusters"] <- "gaa_pausing"


## Pausing Heatmaps

codons_heatmap <- function(tbl, codons, feature){
  out_colnames <- c(t(outer(codons, levels(tbl$site), FUN = "paste", sep = "_")))
  hm <- tbl %>%
    filter(codon_display %in% codons) %>%
    mutate(site_out = paste(codon_display, site, sep = "_")) %>%
    select(CB, site_out, {{ feature }}) %>%
    spread(site_out, {{ feature }}, fill = 0) %>%
    as.data.frame()
  rownames(hm) <- hm$CB
  out_colnames <- out_colnames[out_colnames %in% colnames(hm)]
  hm <- hm[,out_colnames]
  return(hm)
}


all_codons <- unique(depleted_site_change$codon_display)
all_codons <- all_codons[order(gsub(".*_(.*?)$","\\1",all_codons,perl=TRUE))]
all_codons <- all_codons[nchar(all_codons)==7]


# average fc across all sites to cluster codons
allDeplHMa <- depleted_site_change %>%
  filter(codon_display %in% all_codons) %>%
  filter(site %in% c("em2", "em1", "e", "p", "a", "ap1", "ap2")) %>%
  group_by(CB, codon_display) %>%
  summarize(mean_site_fc = mean(site_fc),
            mean_total_fc = mean(total_fc)) %>%
  ungroup() %>%
  select(CB, codon_display, mean_site_fc) %>%
  spread(codon_display, mean_site_fc, fill = 0) %>%
  as.data.frame()
rownames(allDeplHMa) <- allDeplHMa$CB
allDeplHMa <- t(allDeplHMa[,all_codons])
allDeplHMa_hc <- dendsort(hclust(as.dist(1-cosine_similarity(as.matrix(log2(allDeplHMa)))), method = "ward.D2"), type = "min")

all_codons <- all_codons[allDeplHMa_hc$order]

allDeplHM <- depleted_site_change %>%
  codons_heatmap(codons = all_codons, feature = site_fc)
allDeplHM <- t(allDeplHM)
allDeplHM <- allDeplHM[, pull(RPFmeta %>% filter(!is.na(output_clusters)) %>% arrange(output_clusters) %>% select(CB))]
allDeplHM[allDeplHM ==0] <- NA

pdf(paste0(figurePrefix, "_depleted_all_codon_stalling_sitefc_l2fc.pdf"), height = 8, width = 8)
site_labels <- HeatmapAnnotation(site = anno_mark(at = grep("^CAG_Gln", rownames(allDeplHM)),
                                                  labels = c("-2", "-1", "E", "P", "A", "+1", "+2")), which = "row")
cell_annotation <- HeatmapAnnotation(clusters = as.character(RPFmeta[colnames(allDeplHM), "output_clusters"]),
                                     region = as.character(RPFmeta[colnames(markerHM), "sort_population"]),
                                     col = list(clusters = c("1" = gg_color_hue(8)[1],
                                                             "2" = gg_color_hue(8)[2],
                                                             "3" = gg_color_hue(8)[3],
                                                             "4" = gg_color_hue(8)[4],
                                                             "5" = gg_color_hue(8)[5],
                                                             "6" = gg_color_hue(8)[6],
                                                             "7" = gg_color_hue(8)[7],
                                                             "8" = gg_color_hue(8)[8]),
                                                region = c("Distal" = brewer.pal(3, "Dark2")[3],
                                                           "Medial" = brewer.pal(3, "Dark2")[2],
                                                           "Proximal" = brewer.pal(3, "Dark2")[1])),
                                     annotation_name_side = "left")
Heatmap(log2(as.matrix(allDeplHM)),
        right_annotation = site_labels,
        top_annotation = cell_annotation,
        col = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))), # log fold change
        na_col = "grey0",
        name = "log2 fold change",
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        show_column_dend = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        clustering_distance_rows = function(x) as.dist(1-cosine_similarity(x)),
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_split = factor(rep(all_codons, each = length(levels(depleted_site_change$site))), levels = all_codons, ordered = TRUE),
        row_title_rot = 0,
        border = FALSE,
        use_raster = TRUE,
        raster_device = "CairoPNG",
        raster_quality = 15)
dev.off()



selected_codons <- c("GAG_Glu", "CAG_Gln","GAA_Glu")

selectedDeplHM_site <- depleted_site_change %>%
  codons_heatmap(codons = selected_codons, feature = site_fc)
selectedDeplHM_site <- t(selectedDeplHM_site)
selectedDeplHM_site <- selectedDeplHM_site[, pull(RPFmeta %>% filter(!is.na(output_clusters)) %>% arrange(output_clusters) %>% select(CB))]

selectedDeplHM_total <- depleted_site_change %>%
  codons_heatmap(codons = selected_codons, feature = total_fc)
selectedDeplHM_total <- t(selectedDeplHM_total)
selectedDeplHM_total <- selectedDeplHM_total[, pull(RPFmeta %>% filter(!is.na(output_clusters)) %>% arrange(output_clusters) %>% select(CB))]


pdf(paste0(figurePrefix, "_selected_codon_stalling_sitefc_l2fc.pdf"), height = 3, width = 9)
site_labels <- HeatmapAnnotation(site = anno_mark(at = grep("^CAG_Gln", rownames(selectedDeplHM_site)),
                                                  labels = c("-2", "-1", "E", "P", "A", "+1", "+2")), which = "row")
cell_annotation <- HeatmapAnnotation(clusters = as.character(RPFmeta[colnames(selectedDeplHM_site), "output_clusters"]),
                                     cag_pausing = as.character(RPFmeta[colnames(selectedDeplHM_site), "cag_clusters"]),
                                     gaa_pausing = as.character(RPFmeta[colnames(selectedDeplHM_site), "gaa_clusters"]),
                                     col = list(clusters = c("1" = gg_color_hue(8)[1],
                                                             "2" = gg_color_hue(8)[2],
                                                             "3" = gg_color_hue(8)[3],
                                                             "4" = gg_color_hue(8)[4],
                                                             "5" = gg_color_hue(8)[5],
                                                             "6" = gg_color_hue(8)[6],
                                                             "7" = gg_color_hue(8)[7],
                                                             "8" = gg_color_hue(8)[8]),
                                                cag_pausing = c("background" = "white",
                                                                 "cag_pausing" = "black"),
                                                gaa_pausing = c("background" = "white",
                                                                 "gaa_pausing" = "black")
                                                ),
                                     annotation_name_side = "left")
Heatmap(log2(as.matrix(selectedDeplHM_site)),
        right_annotation = site_labels,
        top_annotation = cell_annotation,
        col = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))), # log fold change
        na_col = "grey0",
        name = "log2 fold change",
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        clustering_distance_columns = function(x) as.dist(1-cosine_similarity(x)),
        clustering_method_columns = "complete",
        show_column_dend = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = factor(RPFmeta[colnames(selectedDeplHM_site), "output_clusters"]),
        row_split = factor(rep(selected_codons, each = length(levels(depleted_site_change$site))), levels = selected_codons, ordered = TRUE),
        row_title_rot = 0,
        border = FALSE,
        use_raster = TRUE,
        raster_device = "CairoPNG",
        raster_quality = 15)
dev.off()




highlight_codons <- c("CAG_Gln", "GAA_Glu")

plt_umap_highlight <- depleted_site_change %>%
  filter(codon_display %in% highlight_codons) %>%
  filter(site %in% c("a")) %>%
  inner_join(RPFmeta %>% select(CB, UMAP_1, UMAP_2), by = "CB") %>%
  mutate(arrange_value = site_fc) %>%
  arrange(arrange_value) %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, colour = log2(site_fc)))+
  geom_point(size=0.5)+
  coord_equal()+
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu")), na.value = "grey40", limits = c(-2,2), oob=squish)+
  facet_wrap(~codon_display)+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank())
save_plot(paste0(figurePrefix, "_highlight_codons.pdf"), plt_umap_highlight, base_width = 6, base_height = 2.7)





## statistics for distribution in clusters
plt_cag_cluster_proportions <- RPFmeta %>%
  count(output_clusters, cag_clusters) %>%
  ggplot(aes(x=output_clusters, y=n, fill = cag_clusters))+
  geom_bar(stat = "identity", position = "fill")
plt_gaa_cluster_proportions <- RPFmeta %>%
  count(output_clusters, gaa_clusters) %>%
  ggplot(aes(x=output_clusters, y=n, fill = gaa_clusters))+
  geom_bar(stat = "identity", position = "fill")
save_plot(paste0(figurePrefix, "_plt_gln_cluster_proportions.pdf"), plot_grid(plt_cag_cluster_proportions, plt_gaa_cluster_proportions, nrow = 1), base_width = 7, base_height = 2)

cag_pausing_contingency <- RPFmeta %>%
  count(output_clusters, cag_clusters) %>%
  pivot_wider(names_from = cag_clusters, values_from = n, values_fill = 0) %>%
  column_to_rownames('output_clusters')
fisher.test(cag_pausing_contingency)
gaa_pausing_contingency <- RPFmeta %>%
  count(output_clusters, gaa_clusters) %>%
  pivot_wider(names_from = gaa_clusters, values_from = n, values_fill = 0) %>%
  column_to_rownames('output_clusters')
fisher.test(gaa_pausing_contingency)






## pausing differential expression, subset by cluster 0 and 1
P <- AddMetaData(P, RPFmeta %>% select(CB, cag_clusters) %>% deframe(), col.name = "cag_clusters")
P <- AddMetaData(P, RPFmeta %>% select(CB, gaa_clusters) %>% deframe(), col.name = "gaa_clusters")

CAGP <- subset(P, subset = seurat_clusters == 0)
Idents(object = CAGP) <- CAGP@meta.data$cag_clusters
cag_pausing_markers <- FindAllMarkers(CAGP, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
cag_pausing_markers %>%
  filter(p_val_adj < 0.01) %>%
  count(cluster)

GAAP <- subset(P, subset = seurat_clusters == 7) ## late EC are cluster 7
Idents(object = GAAP) <- GAAP@meta.data$gaa_clusters
gaa_pausing_markers <- FindAllMarkers(GAAP, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE, min.cells.group = 3, min.cells.feature = 3 )
gaa_pausing_markers %>%
  filter(p_val_adj < 0.01) %>%
  count(cluster)





## scatter gene-cluster global change 

gene_codon_background <- codon_reads %>%
  filter(site == "a") %>%
  group_by(gene_id, gene_name, codon) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(gene_name) %>%
  mutate(codon_frac = n/sum(n),
         gene_total = sum(n)) %>%
  ungroup() %>%
  mutate(gene_mean = gene_total / pull(codon_reads %>% select(CB) %>% distinct() %>% summarize(n=n())) ) # divided by number of cells


# CAG
cluster_change <- codon_reads %>%
  filter(site == "a") %>%
  #filter(gene_id %in% pull(gene_codon_background %>% filter(gene_total > 2500) %>% select(gene_id) %>% distinct())) %>%
  filter(gene_id %in% pull(gene_codon_background %>% filter(gene_mean > 1) %>% select(gene_id) %>% distinct())) %>%
  inner_join(RPFmeta %>% group_by(output_clusters, cag_clusters) %>% mutate(cells_per_cluster = n()) %>% ungroup() %>% select(CB, output_clusters, cag_clusters, cells_per_cluster), by = "CB") %>%
  group_by(gene_id, gene_name, codon, output_clusters, cag_clusters, cells_per_cluster) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(gene_name, output_clusters, cag_clusters) %>%
  mutate(cluster_gene_frac = n/sum(n),
         cluster_gene_total = sum(n)) %>%
  ungroup() %>% 
  mutate(cluster_gene_mean = cluster_gene_total / cells_per_cluster) %>%
  inner_join(gene_codon_background %>% select(gene_id, codon, codon_frac), by = c("gene_id" = "gene_id", "codon" = "codon")) %>%
  mutate(lfc_change = log2(cluster_gene_frac/codon_frac))

shared_genes <- pull(cluster_change %>%
  group_by(gene_id, output_clusters, cag_clusters) %>%
  summarize(asites_per_cluster = sum(n)) %>%
  ungroup() %>%
  complete(gene_id, output_clusters, cag_clusters) %>%
  mutate(asites_per_cluster = replace_na(asites_per_cluster,0)) %>%
  inner_join(RPFmeta %>% group_by(output_clusters, cag_clusters) %>% mutate(cells_per_cluster = n()) %>% ungroup() %>% select(output_clusters, cag_clusters, cells_per_cluster) %>% distinct(), by = c("output_clusters", "cag_clusters")) %>%
  mutate(mean_asites = asites_per_cluster / cells_per_cluster) %>%
  filter(mean_asites > 2.5) %>%
  group_by(gene_id) %>%
  summarize(detected_in_clusters = n()) %>%
  ungroup() %>%
  filter(detected_in_clusters > 6) %>%
  select(gene_id) %>%
  distinct())

plot_codon <- "CAG"
gene_codon_change_data <- cluster_change %>%
  filter(codon == plot_codon) %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name == "Chgb" ~ "Chgb",
                                  TRUE ~ ""))
plt_cluster_change_codon_facet_cag <- ggplot(gene_codon_change_data, aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = output_clusters, label = gene_display ))+
  geom_point(size=1)+
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~output_clusters+cag_clusters, nrow=2)+
  geom_text_repel(segment.color = "black", box.padding = 0.4, min.segment.length = 0)+
  geom_text(data = gene_codon_change_data %>% select(output_clusters, cells_per_cluster, cag_clusters) %>% distinct(), 
            mapping = aes(x = -Inf, y=-Inf, label = cells_per_cluster),
            hjust = -0.1,
            vjust = -1)+
  xlab("log10 mean cluster expression")+ylab(paste0("log2 f.c. ", plot_codon, " A-site FOC"))+
  theme(strip.background = element_blank())
save_plot(paste0(figurePrefix, "_plt_cluster_change_gln_codon", plot_codon, ".pdf"), plt_cluster_change_codon_facet_cag, base_width = 10, base_height =5)

# GAA 
cluster_change <- codon_reads %>%
  filter(site == "a") %>%
  filter(gene_id %in% pull(gene_codon_background %>% filter(gene_mean > 1) %>% select(gene_id) %>% distinct())) %>%
  inner_join(RPFmeta %>% group_by(output_clusters, gaa_clusters) %>% mutate(cells_per_cluster = n()) %>% ungroup() %>% select(CB, output_clusters, gaa_clusters, cells_per_cluster), by = "CB") %>%
  group_by(gene_id, gene_name, codon, output_clusters, gaa_clusters, cells_per_cluster) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(gene_name, output_clusters, gaa_clusters) %>%
  mutate(cluster_gene_frac = n/sum(n),
         cluster_gene_total = sum(n)) %>%
  ungroup() %>% 
  mutate(cluster_gene_mean = cluster_gene_total / cells_per_cluster) %>%
  inner_join(gene_codon_background %>% select(gene_id, codon, codon_frac), by = c("gene_id" = "gene_id", "codon" = "codon")) %>%
  mutate(lfc_change = log2(cluster_gene_frac/codon_frac))

shared_genes <- pull(cluster_change %>%
  group_by(gene_id, output_clusters, gaa_clusters) %>%
  summarize(asites_per_cluster = sum(n)) %>%
  ungroup() %>%
  complete(gene_id, output_clusters, gaa_clusters) %>%
  mutate(asites_per_cluster = replace_na(asites_per_cluster,0)) %>%
  inner_join(RPFmeta %>% group_by(output_clusters, gaa_clusters) %>% mutate(cells_per_cluster = n()) %>% ungroup() %>% select(output_clusters, gaa_clusters, cells_per_cluster) %>% distinct(), by = c("output_clusters", "gaa_clusters")) %>%
  mutate(mean_asites = asites_per_cluster / cells_per_cluster) %>%
  filter(mean_asites > 2.5) %>%
  group_by(gene_id) %>%
  summarize(detected_in_clusters = n()) %>%
  ungroup() %>%
  filter(detected_in_clusters > 6) %>%
  select(gene_id) %>%
  distinct())

plot_codon <- "GAA"
gene_codon_change_data <- cluster_change %>%
  filter(codon == plot_codon) %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name == "Chgb" ~ "Chgb",
                                  TRUE ~ ""))
plt_cluster_change_codon_facet_gaa <- ggplot(gene_codon_change_data, aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = output_clusters, label = gene_display ))+
  geom_point(size=1)+
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~output_clusters+gaa_clusters, nrow=2)+
  geom_text_repel(segment.color = "black", box.padding = 0.4, min.segment.length = 0)+
  geom_text(data = gene_codon_change_data %>% select(output_clusters, cells_per_cluster, gaa_clusters) %>% distinct(), 
            mapping = aes(x = -Inf, y=-Inf, label = cells_per_cluster),
            hjust = -0.1,
            vjust = -1)+
  xlab("log10 mean cluster expression")+ylab(paste0("log2 f.c. ", plot_codon, " A-site FOC"))+
  theme(strip.background = element_blank())
save_plot(paste0(figurePrefix, "_plt_cluster_change_gaa_codon", plot_codon, ".pdf"), plt_cluster_change_codon_facet_gaa, base_width = 10, base_height =5)

save_plot(paste0(figurePrefix, "_plt_cluster_change_codons_cag_gaa.pdf"), plot_grid(plt_cluster_change_codon_facet_cag+theme(legend.position = "none"), plt_cluster_change_codon_facet_gaa+theme(legend.position = "none"), nrow = 1), base_width =11.5, base_height = 4)




## Example gene

asite_pausing <- site_change %>%
  filter(codon_display %in% all_codons) %>%
  filter(site %in% c("a")) %>%
  select(CB, codon_display, site_fc) %>%
  pivot_wider(names_from = codon_display, values_from = site_fc) %>%
  as.data.frame()
rownames(asite_pausing) <- asite_pausing$CB


#for(display_gene in annot %>% filter(gene_id %in% shared_genes) %>% filter(set == "canonical") %>% select(gene_name) %>% distinct() %>% pull()){

display_gene = "Chgb"

EXGcodons <- codons %>%
  filter(transcript_id == pull(annot %>% filter(gene_name == display_gene, set == "canonical") %>% select(transcript_id)))

EXG <- EXGcodons %>%
  inner_join(codon_reads %>% filter(site %in% c("a")), by = c("transcript_id" = "transcript_id", "codon_position" = "site_position")) %>%
  group_by(codon_position,CB) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  complete(codon_position = full_seq(c(min(EXGcodons$codon_position), max(EXGcodons$codon_position)),3), CB = pull(RPFmeta %>% select(CB) %>% distinct())) %>%
  mutate(n = replace_na(n,0)) %>%
  group_by(CB) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  inner_join(EXGcodons, by="codon_position") %>%
  inner_join(RPFmeta, by="CB") %>%
  mutate(cell_frac = n/tot)

EXG_hm <- EXG %>%
  select(CB, codon_position, cell_frac) %>%
  spread(codon_position, cell_frac, fill = 0) %>%
  as.data.frame()
rownames(EXG_hm) <- EXG_hm$CB
EXG_hm <- EXG_hm[,colnames(EXG_hm) != "CB"]
EXG_hm <- EXG_hm[pull(RPFmeta %>% filter(CB %in% rownames(EXG_hm)) %>% arrange(output_clusters) %>% select(CB)), ]

EXG_hm <- EXG_hm[apply(EXG_hm>0,1,sum)>(0.005*ncol(EXG_hm)),]

EXG_hm <- EXG_hm[pull(asite_pausing %>% filter(CB %in% rownames(EXG_hm)) %>% arrange(desc(GAA_Glu)) %>% select(CB)), ]

aa_highlight <- rep("other", nrow(EXGcodons))
aa_highlight[grep("CAG", EXGcodons$codon)] <- "CAG"
aa_highlight[grep("GAA", EXGcodons$codon)] <- "GAA"
pausing_clusters <- RPFmeta %>%
  mutate(pausing_clusters = case_when( (cag_clusters == "cag_pausing" & gaa_clusters == "gaa_pausing") ~ "both",
                                      cag_clusters == "cag_pausing" ~ "cag_pausing",
                                      gaa_clusters == "gaa_pausing" ~ "gaa_pausing",
                                      TRUE ~ "background") ) %>%
  select(CB, pausing_clusters) %>%
  deframe()

pdf(paste0(figurePrefix, "_pausing_hm_", display_gene, ".pdf"), height = 3.5, width = 12)
aa_annotation <- HeatmapAnnotation(aa = anno_text(EXGcodons$single, gp = gpar(fontsize = 2)))
aa_show <- HeatmapAnnotation(aa = aa_highlight,
                             col = list(aa = c("other" = "white",
                                               "CAG" = "#e66101",
                                               "GAA" = "#5e3c99"
                                               )),
                             annotation_name_side = "left")
exprn <- as.numeric(P@assays$RNA[ pull(annot %>% 
                                       filter(gene_name == display_gene) %>% 
                                       mutate(gene_out = paste(gene_id, gene_name, sep = '-')) %>% 
                                       select(gene_out) %>% 
                                       distinct()), rownames(EXG_hm) ] )
cell_annotation <- HeatmapAnnotation(exprn = exprn,
                                     CAG = log2(asite_pausing[rownames(EXG_hm), "CAG_Gln"]),
                                     GAA = log2(asite_pausing[rownames(EXG_hm), "GAA_Glu"]),
                                     clusters = as.character(RPFmeta[rownames(EXG_hm), "output_clusters"]),
                                     region = as.character(RPFmeta[rownames(EXG_hm), "sort_population"]),
                                     col = list(clusters = c("1" = gg_color_hue(8)[1],
                                                             "2" = gg_color_hue(8)[2],
                                                             "3" = gg_color_hue(8)[3],
                                                             "4" = gg_color_hue(8)[4],
                                                             "5" = gg_color_hue(8)[5],
                                                             "6" = gg_color_hue(8)[6],
                                                             "7" = gg_color_hue(8)[7],
                                                             "8" = gg_color_hue(8)[8]),
                                                region = c("Distal" = brewer.pal(3, "Dark2")[3],
                                                           "Medial" = brewer.pal(3, "Dark2")[2],
                                                           "Proximal" = brewer.pal(3, "Dark2")[1]),
                                                red = colorRamp2(seq(0, max(RPFmeta$log_red, na.rm = TRUE), length.out = 8), brewer.pal(8, 'Reds')),
                                                green = colorRamp2(seq(0, max(RPFmeta$log_green, na.rm = TRUE), length.out = 8), brewer.pal(8, 'Greens')),
                                                CAG = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))),
                                                GAA = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))),
                                                exprn = colorRamp2(seq(min(exprn), max(exprn), length.out = 11), rev(brewer.pal(11, "Spectral")))),
                                     which = "row",
                                     annotation_name_side = "top")
length_annotation <- HeatmapAnnotation(pos = anno_mark(at = c(seq(50,ncol(EXG_hm), by = 100)),
                                                       labels = as.character(c(seq(50,ncol(EXG_hm), by = 100))), 
                                                       which = "column", side = "bottom", labels_rot = 0, link_height = unit(3, "mm")))
hm_quantiles <- 0.01
row_dend <- dendsort(hclust(as.dist(1-cosine_similarity(as.matrix(EXG_hm))), method = "ward.D2"), type = "average")
Heatmap(as.matrix(EXG_hm),
        top_annotation = aa_show,
        right_annotation = cell_annotation,
        bottom_annotation = length_annotation,
        col = colorRamp2(seq(quantile(as.matrix(EXG_hm), probs = hm_quantiles), quantile(as.matrix(EXG_hm), probs = 1-hm_quantiles), length = 9), brewer.pal(9, "Greys")),
        na_col = "grey0",
        name = "fraction of readspcell",
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        row_dend_reorder = FALSE,
        clustering_distance_rows = function(x) as.dist(1-cosine_similarity(x)),
        clustering_method_rows = "complete",
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_split = factor(pausing_clusters[rownames(EXG_hm)], levels = c("cag_pausing", "gaa_pausing", "background"), ordered = "TRUE"),
        row_title = NULL,
        column_title_rot = 0,
        border = TRUE,
        use_raster = FALSE,
        raster_device = "CairoPNG",
        raster_quality = 20)
dev.off()


topgenes <- RPFlong %>%
  filter(gene_id %in% shared_genes) %>%
  inner_join(RPFmeta %>% select(CB, output_clusters), by = "CB") %>%
  group_by(gene_id, output_clusters) %>%
  summarize(mean_n = mean(counts)) %>%
  ungroup() %>%
  arrange(desc(mean_n)) %>%
  inner_join(annot %>% filter(set == "canonical") %>% select(gene_id, gene_name) %>% distinct(), by = "gene_id")

topgenes %>% count(gene_id) %>% filter(n>7) %>% inner_join(topgenes %>% select(gene_id, gene_name, mean_n)) %>% arrange(desc(mean_n))




## output for GEO

write.csv(RPFmeta, file = paste0(figurePrefix, "_EE_metadata.csv"))
write.csv(RPFcounts, file = paste0(figurePrefix, "_EE_RPFcounts.csv"))

Distal <- barcodes
Distal$fraction <- "Distal"
Distal[Distal$well %in% columnFeatures(21:24, LETTERS[15:16]), "fraction"] = "NTC"
write.csv(Distal, file = paste0(figurePrefix, "_EE_Distal_layout.csv"), row.names = FALSE, quote = FALSE)

Medial <- barcodes
Medial$fraction <- "Medial"
Medial[Medial$well %in% columnFeatures(21:24, LETTERS[15:16]), "fraction"] = "NTC"
write.csv(Medial, file = paste0(figurePrefix, "_EE_Medial_layout.csv"), row.names = FALSE, quote = FALSE)

Proximal <- barcodes
Proximal$fraction <- "Proximal"
Proximal[Proximal$well %in% columnFeatures(21:24, LETTERS[15:16]), "fraction"] = "NTC"
write.csv(Proximal, file = paste0(figurePrefix, "_EE_Proximal_layout.csv"), row.names = FALSE, quote = FALSE)
