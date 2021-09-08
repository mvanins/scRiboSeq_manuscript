#!/usr/bin/env Rscript
# Michael VanInsberghe 2020-06-03
# Produces raw plots for figure 2 - starvation

library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(data.table)
library(Matrix)
library(Matrix.utils)
library(Seurat)
library(ggplot2)
library(scales)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(forcats)
library(RColorBrewer)
suppressPackageStartupMessages(library(ggseqlogo))
theme_set(theme_cowplot())

source('plotFunctions.R')

data_dir <- "../scRiboSeq/data"

figurePrefix <- "F2_updated/F2_updated"

coltreat <- c("#a6cee3", "#1f78b4",
              "#b2df8a", "#33a02c",
              "grey50")
names(coltreat) <- c("Arg3h", "Arg6h",
                     "Leu3h", "Leu6h",
                     "Rich")

# barcodes
barcodes <- read.table('/hpc/hub_oudenaarden/mvanins/local/lib/barcodes/miRv6_barcodemap.tsv',
                       sep="\t",
                       header = FALSE,
                       stringsAsFactors = FALSE,
                       row.names = 1,
                       col.names = c("","number","well"))
rownames(barcodes) <- barcodes$CB

# annot
annot <- fread(file.path(data_dir, 'hsa_annotations.csv.gz'), stringsAsFactors = FALSE, data.table = FALSE)
canonical_pc <- annot %>%
  filter(set=="canonical", transcript_type == "protein_coding") %>%
  filter(transcript_id != "ENST00000437139.7") %>%
  filter(!grepl("^PCDH", gene_name)) 


# codons
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


#
codon_reads <- hekReads %>%
  select(-starts_with("base")) %>%
  inner_join(hekMeta, by = "CB") %>%
  rowwise() %>% # expand each read to each site (rowwise, mutate, unnest)
    mutate(site = list(c("em2","em1","e","p","a","ap1","ap2")),
           site_position = list(c(psite-9, psite-6, psite-3, psite, psite+3, psite+6, psite+9))) %>%
  unnest(cols = c(site,site_position)) %>%
  select(-psite) %>% # remove p-site column
  inner_join(codons, by = c("transcript_id" = "transcript_id", "site_position" = "codon_position")) # join codon information. Will only return site_position that exactly match codon_position, i.e., in frame

 
### relative change between rich and starvation conditions

baseline <- codon_reads %>%
  mutate(codon_display = paste(gsub("T", "U", codon),three,sep="_")) %>%
  filter(sort_population == "Rich",
         site %in% c("p","a"),
         !(three %in% c("Met","Ter"))) %>%
  group_by(CB,codon_display) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(CB) %>%
  mutate(cell_total = sum(n)) %>%
  ungroup() %>%
  mutate(cell_fraction = n/cell_total) %>%
  group_by(codon_display) %>%
  summarize(baseline_fraction = mean(cell_fraction)) %>%
  ungroup()

rel_change <- codon_reads %>%
  mutate(codon_display = paste(gsub("T", "U", codon),three,sep="_")) %>%
  filter(site %in% c("p","a"),
         !(three %in% c("Met","Ter"))) %>%
  group_by(CB, sort_population, codon_display) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(CB) %>%
  mutate(cell_total = sum(n)) %>%
  ungroup() %>%
  mutate(cell_fraction = n/cell_total) %>%
  group_by(sort_population, codon_display) %>%
  summarize(codon_fraction = mean(cell_fraction)) %>%
  ungroup() %>%
  inner_join(baseline, by = "codon_display") %>%
  mutate(relative_change = 100*(codon_fraction - baseline_fraction)/baseline_fraction,
         fold_change = log2(codon_fraction/baseline_fraction)) %>%
  #   select(hekMeta, codon_display, relative_change) %>%
  select(sort_population, codon_display, fold_change) %>%
  filter(sort_population != "Rich") %>%
  separate(sort_population, into = c("starvation", "time"), sep = 3) %>%
  #   spread(time, relative_change) %>%
  spread(time, fold_change) %>%
  mutate(codon_display = case_when(codon_display %in% c("CGC_Arg", "CGU_Arg", "UUA_Leu") ~ codon_display,
                                   TRUE ~ ""))

plt_change <- ggplot(rel_change, aes(x=`3h`, y=`6h`, colour = starvation, label = codon_display))+
  geom_point()+
  #scale_colour_manual(values = c("#33a02c", "#1f78b4"))+
  scale_colour_manual(values = as.character(coltreat[c(2,4)]))+
  geom_text_repel(segment.color = "black", box.padding = 0.8, seed = 3)+
  coord_equal()+
  scale_x_continuous(breaks=c(-0.25, 0, 0.25, 0.5, 0.75)) +
  scale_y_continuous(breaks=c(-0.25, 0, 0.25, 0.5, 0.75)) +
  xlab("log2 fold change in codon occupancy after 3h starvation [%]")+
  ylab("log2 fold change in codon occupancy after 6h starvation [%]")

save_plot(paste0(figurePrefix, "_fold_change.pdf"), plt_change, base_width = 5, base_height = 5)
# save_plot(paste0(figurePrefix, "_relative_change.pdf"), plt_change, base_width = 6, base_height = 6)
# save_plot(paste0(figurePrefix, "_relative_change_all.pdf"), plt_change, base_width = 12, base_height = 12)

### site-specific heatmaps

site_order <- data.frame(site = c("em2", "em1", "e", "p", "a", "ap1", "ap2"),
                         site_order = c(1, 2, 3, 4, 5, 6, 7),
                         stringsAsFactors = FALSE)

baseline_counts <- codon_reads %>%
  filter(sort_population == "Rich",
         !(three %in% c("Ter"))) %>%
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
  inner_join(hekMeta %>% select(CB, sort_population), by = "CB") %>%
  inner_join(total_baseline, by = "codon_display") %>%
  inner_join(site_baseline, by = c("codon_display" = "codon_display", "site" = "site")) %>%
  mutate(total_fc = cell_total_fraction/(baseline_total_fraction/nrow(site_order)),
         site_fc = cell_site_fraction/baseline_site_fraction) %>%
  inner_join(site_order, by = "site") %>%
  mutate(site = fct_reorder(site, site_order))
  

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

pdf(paste0(figurePrefix, "_selected_codon_stalling_total_l2fc.pdf"), height = 5, width = 6)
site_labels <- HeatmapAnnotation(site = anno_text(rep(c("-2","-1","E","P","A","+1","+2"), times = 3), rot=90))
Heatmap(log2(as.matrix(selected_hm)),
        col = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))), # log fold change
        bottom_annotation = site_labels,
        na_col = "grey0",
        name = "log2 fold change",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_split = factor(hekMeta[rownames(selected_hm), "sort_population"], levels = c("Rich","Arg3h", "Arg6h", "Leu3h", "Leu6h"), ordered = TRUE),
        column_split = factor(rep(display_codons, each = length(levels(site_change$site))), levels = display_codons, ordered = TRUE),
        column_title_rot = 0,
        use_raster = TRUE,
        raster_device = "CairoPNG",
        raster_quality = 15,
        border = TRUE)
dev.off()

all_codons <- unique(site_change$codon_display)
all_codons <- all_codons[order(gsub(".*_(.*?)$","\\1",all_codons,perl=TRUE))]

rich_all <- site_change %>%
  codons_heatmap(codons = all_codons,  
                 treatments = "Rich",
                 feature = total_fc)
leu3h_all <- site_change %>%
  codons_heatmap(codons = all_codons, 
                 treatments = "Leu3h",
                 feature = total_fc)
leu6h_all <- site_change %>%
  codons_heatmap(codons = all_codons, 
                 treatments = "Leu6h",
                 feature = total_fc)
arg3h_all <- site_change %>%
  codons_heatmap(codons = all_codons, 
                 treatments = "Arg3h",
                 feature = total_fc)
arg6h_all <- site_change %>%
  codons_heatmap(codons = all_codons, 
                 treatments = "Arg6h",
                 feature = total_fc)

all_hm <- rbind(rich_all,
                arg3h_all,
                arg6h_all,
                leu3h_all,
                leu6h_all)

pdf(paste0(figurePrefix, "_all_codon_stalling_total_l2fc.pdf"), height = 10, width = 20)
site_labels <- HeatmapAnnotation(site = anno_mark(at = grep("^CGC_Arg", colnames(all_hm)),
                                                  labels = c("-2", "-1", "E", "P", "A", "+1", "+2")), which = "column")
Heatmap(log2(as.matrix(all_hm)),
        bottom_annotation = site_labels,
        col = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))), # log fold change
        #col = colorRamp2(seq(0, 5, length = 11), rev(brewer.pal(11, "RdBu"))),
        na_col = "grey0",
        name = "log2 fold change",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_split = factor(hekMeta[rownames(all_hm), "sort_population"], levels = c("Rich","Arg3h", "Arg6h", "Leu3h", "Leu6h"), ordered = TRUE),
        column_split = factor(rep(all_codons, each = length(levels(site_change$site))), levels = all_codons, ordered = TRUE),
        column_title_rot = 90,
        use_raster = TRUE,
        raster_device = "CairoPNG",
        raster_quality = 15,
        border = TRUE)
dev.off()




### Seurat
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

plt_variable <- VariableFeaturePlot(P)


P <- ScaleData(P, features = rownames(P))
P <- RunPCA(P, features = VariableFeatures(object = P))

P <- JackStraw(P, num.replicate = 100)
P <- ScoreJackStraw(P, dims = 1:20)

plt_jackstraw <- JackStrawPlot(P,dims=1:20)
plt_elbow <- ElbowPlot(P)

#save_plot(paste0(figurePrefix,"_qc_plt_dimscores.pdf"), plot_grid(plt_jackstraw, plt_elbow, nrow=1), base_width = 10, base_height = 5)

ndims <- 5 
P <- RunUMAP(P, dims = 1:ndims, verbose = FALSE, seed.use = 121, n.epochs = 500)#, n.epochs = 3000)
P <- FindNeighbors(P, dims = 1:ndims)
P <- FindClusters(P, resolution = 0.4)

plt_umap_plate <- DimPlot(P, reduction = "umap",group.by = "plate",pt.size=1)
plt_umap_cluster <- DimPlot(P, reduction = "umap",pt.size=1)
plt_umap_treatment <- DimPlot(P, reduction = 'umap', group.by = "sort_population", pt.size=1)


save_plot(paste0(figurePrefix,"_plt_umaps_check.pdf"), plot_grid(plt_umap_plate+coord_equal(), plt_umap_treatment+coord_equal() ,plt_umap_cluster+coord_equal(), nrow=2), base_width = 10, base_height = 10)

hekMeta$seurat_clusters <- P$seurat_clusters[rownames(hekMeta)]
hekMeta$UMAP_1 <- as.data.frame(Embeddings(P, reduction="umap"))[rownames(hekMeta),"UMAP_1"]
hekMeta$UMAP_2 <- as.data.frame(Embeddings(P, reduction="umap"))[rownames(hekMeta),"UMAP_2"]

saveRDS(hekMeta, file.path(data_dir, "HEK293T/compiled/clustered_meta.rds"))


mean_site_fc <- site_change %>% 
  filter((site %in% c("e", "p", "a"))) %>%
  group_by(CB, codon_display) %>%
   summarize(mean_site_fc = mean(site_fc)) %>%
  ungroup() %>%
  spread(codon_display, mean_site_fc, fill = 0) %>%
  as.data.frame()
rownames(mean_site_fc) <- mean_site_fc$CB
mean_site_fc <- mean_site_fc[,!(colnames(mean_site_fc) == "CB")]


arg_pausing <- site_change %>%
  filter(site %in% c("e", "p", "a")) %>%
  filter(codon_display %in% c("CGC_Arg", "CGU_Arg")) %>%
  group_by(CB) %>%
  summarize(mean_arg = mean(site_fc)) %>%
  ungroup()

leu_pausing <- site_change %>%
  filter(site %in% c("p","a")) %>%
  filter(codon_display %in% c("UUA_Leu")) %>%
  group_by(CB) %>%
  summarize(mean_leu = mean(site_fc)) %>%
  ungroup()

pausing <- arg_pausing %>%
  full_join(leu_pausing, by = "CB") %>%
  as.data.frame()
rownames(pausing) <- pausing$CB
pausing <- pausing[, colnames(pausing) != "CB"]

U <- hekMeta
U <- cbind(U, 
           pausing[rownames(U),],
           mean_site_fc[rownames(U),])

plt_utc <- ggplot(U, aes(x=UMAP_1, y=UMAP_2, colour = sort_population))+
  geom_point(size=2)+
  geom_density_2d(aes(x=UMAP_1, y=UMAP_2, colour = seurat_clusters), data = U, binwidth = 0.005, alpha = 0.8, linetype = "dashed")+
  coord_equal()+
  scale_colour_manual(values = c(gg_color_hue(4), as.character(coltreat)))+
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  expand_limits(x=c(min(U$UMAP_1)-1, max(U$UMAP_1)+1),y=c(min(U$UMAP_2)-1, max(U$UMAP_2)+1))
plt_utc_gg <- ggplot_build(plt_utc)
plt_utc_gg$data[[2]] <- subset(plt_utc_gg$data[[2]], level == 0.02)
plt_utc_out <- ggplot_gtable(plt_utc_gg)
save_plot(paste0(figurePrefix, "_plt_hekMeta_cluster_outline.pdf"), plt_utc_out)

  



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

plt_CGC_Arg <- featplot(U, "CGC_Arg")
plt_CGU_Arg <- featplot(U, "CGU_Arg")
plt_UUA_Leu <- featplot(U, "UUA_Leu")

save_plot(paste0(figurePrefix, "_plt_aas_umap.pdf"),
                 plot_grid(plt_CGC_Arg,
                           plt_CGU_Arg,
                           plt_UUA_Leu, nrow=1),
                 base_width = 16,
                 base_height = 16/3)

plt_mean_arg <- featplot(U, "mean_arg")
plt_mean_leu <- featplot(U, "mean_leu")
save_plot(paste0(figurePrefix, "_plt_mean_aas_umap_outline.pdf"),
                 plot_grid(plt_mean_arg,
                           plt_mean_leu, nrow=1),
                 base_width = 10*.7,
                 base_height = 5*.7)

plt_treatment_facet <- ggplot(U, aes(x=UMAP_1, y=UMAP_2, colour = sort_population))+
  geom_point(data = U %>% select(UMAP_1, UMAP_2), colour = "grey80") +
  geom_point()+
  coord_equal()+
  facet_wrap(~sort_population)

#save_plot(paste0(figurePrefix, "_plt_umap_treatment_facet.pdf"), plt_treatment_facet)

# threshold to find which cells are pausing
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
  






plt_hist_codon <- ggplot(pausing_signal, aes(x=log2(pausing), fill = sort_population))+
  geom_histogram(bins = 50)+
  scale_fill_manual(values = as.character(coltreat))+
  facet_wrap(~codon_display+seurat_clusters)+
  geom_vline(xintercept = log2(pausing_threshold$threshold), colour = "red", linetype = "dashed")
save_plot(paste0(figurePrefix, "_hist_pausing.pdf"), plt_hist_codon, base_width = 10, base_height = 5)


## numbers of cells that paused:
total_pausing <- pausing_signal %>%
  separate(codon_display, c("codon", "aa"), remove = FALSE, sep = "_") %>%
  group_by(CB, aa) %>%
  mutate(max_pausing = max(pausing)) %>%
  ungroup()



total_pausing_leu <- pausing_signal %>%
  filter(sort_population %in% c("Leu3h", "Leu6h")) %>%
  filter(codon_display %in% c("UUA_Leu")) %>%
  group_by(CB) %>%
  summarize(paused = max(pausing) > pausing_threshold$threshold) %>%
  ungroup() %>%
  summarize(n_cells = n(),
            paused_cells = sum(paused))

cluster_pausing_leu <- pausing_signal %>%
  filter(sort_population %in% c("Leu3h", "Leu6h")) %>%
  filter(codon_display %in% c("UUA_Leu")) %>%
  group_by(seurat_clusters, CB) %>%
  summarize(paused = max(pausing) > pausing_threshold$threshold) %>%
  ungroup() %>%
  group_by(seurat_clusters) %>%
  summarize(n_cells = n(),
            paused_cells = sum(paused)) %>%
  ungroup()

total_pausing_arg <- pausing_signal %>%
  filter(sort_population %in% c("Arg3h", "Arg6h")) %>%
  filter(codon_display %in% c("CGC_Arg", "CGU_Arg")) %>%
  group_by(CB) %>%
  summarize(paused = max(pausing) > pausing_threshold$threshold) %>%
  ungroup() %>%
  summarize(n_cells = n(),
            paused_cells = sum(paused))

cluster_pausing_arg <- pausing_signal %>%
  filter(sort_population %in% c("Arg3h", "Arg6h")) %>%
  filter(codon_display %in% c("CGC_Arg", "CGU_Arg")) %>%
  group_by(seurat_clusters, CB) %>%
  summarize(paused = max(pausing) > pausing_threshold$threshold) %>%
  ungroup() %>%
  group_by(seurat_clusters) %>%
  summarize(n_cells = n(),
            paused_cells = sum(paused)) %>%
  ungroup()


markers <- FindAllMarkers(P, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% filter(p_val_adj < 0.01) %>% top_n(n = 10, wt = avg_logFC)

## Marker HM figure

HM <- cbind(U, as.data.frame(t(P@assays$RNA[,rownames(U)])))

topMarkers <- markers %>% 
  filter(p_val_adj < 0.01) %>% 
  filter(!grepl("-AC[0-9][0-9]+", gene)) %>% 
  filter(!grepl("-AL[0-9][0-9]+", gene)) %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_logFC) 


U$CB <- rownames(U)
markerHM <- t(  HM[pull(U %>%
                        arrange(seurat_clusters) %>%
                        select(CB)),
              pull(topMarkers %>%
                   ungroup() %>%
                   select(gene) %>%
                   distinct()) ] ) 

markerHM <- sweep(markerHM, 1, apply(markerHM,1,mean), FUN = "-")
markerHM <- sweep(markerHM, 1, apply(markerHM,1,sd), FUN = "/")

# pdf(paste0(figurePrefix, "_cc_marker_hm.pdf"), height = 8, width = 8)

clusterCols <- gg_color_hue(length(unique(U$seurat_clusters)))
names(clusterCols) <- unique(U$seurat_clusters)

pdf(paste0(figurePrefix, "_marker_hm.pdf"), height = 3.3, width = 8)
marker_annotation <- HeatmapAnnotation(seurat_clusters = as.character(U[colnames(markerHM), "seurat_clusters"]),
                                       treatment = as.character(U[colnames(markerHM), "sort_population"]),
                                       mean_arg = log2(U[colnames(markerHM), "mean_arg"]),
                                       mean_leu = log2(U[colnames(markerHM), "mean_leu"]),
                                       col = list(seurat_clusters = c("0" = gg_color_hue(4)[1],
                                                                      "1" = gg_color_hue(4)[2],
                                                                      "2" = gg_color_hue(4)[3],
                                                                      "3" = gg_color_hue(4)[4]),
                                                  treatment = c("Arg3h" = as.character(coltreat[1]),
                                                                "Arg6h" = as.character(coltreat[2]),
                                                                "Leu3h" = as.character(coltreat[3]),
                                                                "Leu6h" = as.character(coltreat[4]),
                                                                "Rich" =  as.character(coltreat[5])),
                                                  mean_arg = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))), # log fold change
                                                  mean_leu = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))) ) )
gene_annotation <- rowAnnotation(foo = anno_mark(at = c(
                                                        grep("-H3C2$", rownames(markerHM)),
                                                        grep("-H2BC5$", rownames(markerHM)),
                                                        grep("-MCM5$", rownames(markerHM)),
                                                        grep("-DTL$", rownames(markerHM)),
                                                        grep("-UNG$", rownames(markerHM)),
                                                        grep("-CCNB1$", rownames(markerHM)),
                                                        grep("-CENPE$", rownames(markerHM)),
                                                        grep("-CDC20$", rownames(markerHM))),
                                                 labels = c(
                                                        "H3C2", 
                                                        "H2BC5", 
                                                        "MCM5", 
                                                        "DTL", 
                                                        "UNG", 
                                                        "CCNB1",
                                                        "CENPE",
                                                        "CDC20")))
Heatmap(as.matrix(markerHM),
        top_annotation = marker_annotation,
        right_annotation = gene_annotation,
        na_col = "grey0",
        name = "z-score expression",
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        clustering_distance_columns = "pearson",
        clustering_method_rows = "ward.D2",
        cluster_columns = TRUE,
        show_column_dend = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_split = as.factor(U[colnames(markerHM), "seurat_clusters"]),
        column_title_rot = 90,
        border = TRUE)
dev.off()



## pausing heatmap

H3C2codons <- codons %>%
  filter(transcript_id == as.character(annot[annot$gene_name == "H3C2", "transcript_id"]))

H3C2 <- H3C2codons %>%
  inner_join(codon_reads %>% filter(site %in% c("p")) %>% filter(CB %in% pull(U %>% filter(seurat_clusters %in% c(0)) %>% select(CB))), by = c("transcript_id" = "transcript_id", "codon_position" = "site_position")) %>%
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
H3C2_hm <- H3C2_hm[pull( U %>% filter(CB %in% rownames(H3C2_hm)) %>% arrange(-mean_arg) %>% select(CB)),]
H3C2_hm <- H3C2_hm[!grepl("^NA",rownames(H3C2_hm)),]

argSortSplit <- data.frame(CB = pull(U %>%
                                     filter(CB %in% rownames(H3C2_hm)) %>%
                                     arrange(-mean_arg) %>%
                                     select(CB)),
                           sortSplit = c(rep("0",25),
                                         rep("1", 25),
                                         rep("2", 25),
                                         rep("3", 25),
                                         rep("rest", nrow(H3C2_hm)-100)),
                           stringsAsFactors = FALSE)
rownames(argSortSplit) <- argSortSplit$CB


pdf(paste0(figurePrefix, "_pausing_hm.pdf"), height = 3, width = 8)
aa_annotation <- HeatmapAnnotation(aa = anno_text(H3C2codons$single, gp = gpar(fontsize = 4)))
cell_annotation <- rowAnnotation(#seurat_clusters = as.character(U[rownames(H3C2_hm), "seurat_clusters"]),
                                 treatment = as.character(U[rownames(H3C2_hm), "sort_population"]),
                                 mean_arg = log2(U[rownames(H3C2_hm), "mean_arg"]),
                                 #mean_leu = log2(U[rownames(H3C2_hm), "mean_leu"]),
                                 col = list(seurat_clusters = c("0" = gg_color_hue(4)[1],
                                                                "1" = gg_color_hue(4)[2],
                                                                "2" = gg_color_hue(4)[3],
                                                                "3" = gg_color_hue(4)[4]),
                                            treatment = c("Arg3h" = as.character(coltreat[1]),
                                                          "Arg6h" = as.character(coltreat[2]),
                                                          "Leu3h" = as.character(coltreat[3]),
                                                          "Leu6h" = as.character(coltreat[4]),
                                                          "Rich" =  as.character(coltreat[5])),
                                            mean_arg = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))), # log fold change
                                            mean_leu = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))) ) )
length_annotation <- HeatmapAnnotation(pos = anno_mark(at = c(seq(25,ncol(H3C2_hm), by = 25)), labels = as.character(c(seq(25,ncol(H3C2_hm), by = 25))), which = "column", side = "bottom", labels_rot = 0, link_height = unit(3, "mm")))
Heatmap(as.matrix(H3C2_hm),
        right_annotation = cell_annotation,
        bottom_annotation = length_annotation,
        col = colorRamp2(seq(0,0.001, length = 9), (brewer.pal(9, "Greys"))),
        na_col = "grey0",
        name = "fraction of reads/cell",
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_title = NULL,
        column_title_rot = 90,
        border = TRUE)
dev.off()


# zoom in barplot
zoom_region <- seq(195-7*3, 198+7*3,by=3)
H3C2z <- H3C2 %>%
  filter(codon_position %in% zoom_region) %>%
  inner_join(argSortSplit, by = "CB") %>%
  mutate(codon_pos_display = paste(codon_position, gsub("T", "U", codon), single, sep="_")) %>%
  group_by(codon_pos_display, sortSplit) %>%
  summarize(mean_cell_frac = mean(cell_frac)) %>%
  ungroup()

plt_H3C2_zoom <- ggplot(H3C2z, aes(x=codon_pos_display, y=mean_cell_frac, fill = sortSplit))+
  geom_bar(stat="identity", position=position_dodge())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_x_discrete(breaks = pull(H3C2z %>% select(codon_pos_display) %>% unique()),
                   labels = pull(H3C2z %>% select(codon_pos_display) %>% unique() %>% separate(codon_pos_display, c("pos", "codon", "single"), sep = "_") %>% mutate(out = paste(codon, single, sep="_")) %>% select(out))) +
  ylab("mean fraction reads/cell")+
  scale_fill_manual(values = c("#d55e00", "#2271B2", "#359B73", "#E69F00", "grey30"))+
  theme(axis.title.x = element_blank())
save_plot(paste0(figurePrefix, "_H3C2_zoom_bar.pdf"), plt_H3C2_zoom, base_width = 6, base_height =3.5)




## old pausing heatmap

H3C2codons <- codons %>%
  filter(transcript_id == as.character(annot[annot$gene_name == "H3C2", "transcript_id"]))

H3C2 <- H3C2codons %>%
  inner_join(codon_reads %>% filter(site %in% c("p")), by = c("transcript_id" = "transcript_id", "codon_position" = "site_position")) %>%
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
  inner_join(HEKStats, by = "CB") %>%
  mutate(cell_frac = n/cell_total)


H3C2_hm <- H3C2 %>%
  select(CB, codon_position, cell_frac) %>%
  spread(codon_position, cell_frac, fill = 0) %>%
  as.data.frame()
rownames(H3C2_hm) <- H3C2_hm$CB
H3C2_hm <- H3C2_hm[,colnames(H3C2_hm) != "CB"]
H3C2_hm <- H3C2_hm[pull(arg_pausing[order(-arg_pausing$mean_arg),],CB),]
H3C2_hm <- H3C2_hm[!grepl("^NA",rownames(H3C2_hm)),]

pdf(paste0(figurePrefix, "_H3C2_heatmap_labelled.pdf"), height = 8, width = 15)
aas <- data.frame(display_codon = paste(gsub("T","U",H3C2codons$codon), H3C2codons$three, sep = "_"),
                  row.names = H3C2codons$codon_position,
                  stringsAsFactors = FALSE)

aas$out_category <- "other"
aas[grep("Arg",aas$display_codon),"out_category"] <- aas[grep("Arg",aas$display_codon),"display_codon"]
aas[grep("Leu",aas$display_codon),"out_category"] <- aas[grep("Leu",aas$display_codon),"display_codon"]

aacols <- rep("white", length(unique(aas$out_category)))
names(aacols) <- unique(aas$out_category)
aacols["CGC_Arg"] = rev(brewer.pal(6, "Reds"))[1]
aacols["CGU_Arg"] = rev(brewer.pal(6, "Reds"))[2] 
aacols["CGG_Arg"] = rev(brewer.pal(6, "Reds"))[3]
aacols["CGA_Arg"] = rev(brewer.pal(6, "Reds"))[4]
aacols["AGA_Arg"] = rev(brewer.pal(6, "Reds"))[5]

aacols["CUG_Leu"] = rev(brewer.pal(4, "Purples"))[1]
aacols["UUG_Leu"] = rev(brewer.pal(4, "Purples"))[2]
aacols["CUU_Leu"] = rev(brewer.pal(4, "Purples"))[3]
aacols["CUC_Leu"] = rev(brewer.pal(4, "Purples"))[4]
#aacols["UUA_Leu"] = "#08519c"
aa_annotation <- HeatmapAnnotation(aa = aas[colnames(H3C2_hm),"out_category"],
                                   col = list(aa = aacols))

Heatmap(as.matrix(H3C2_hm),
        top_annotation = aa_annotation,
        col = colorRamp2(seq(0,0.002, length = 11), rev(brewer.pal(11, "RdYlBu"))),
        na_col = "grey0",
        name = "log2 fold change",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = TRUE,
        show_row_names = FALSE,
        row_split = factor(hekMeta[rownames(H3C2_hm), "hekMeta"], levels = c("Rich","Arg3h", "Arg6h", "Leu3h", "Leu6h"), ordered = TRUE),
        column_title_rot = 90,
        border = TRUE)
dev.off()




plt_H3C2_sc <- ggplot(H3C2, aes(x=codon_position,y=CB,fill=n/cell_total)) +
  geom_tile()+
  scale_fill_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), na.value = "grey40") +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_discrete(expand=c(0,0)) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank()) +
         facet_wrap(~hekMeta, ncol = 1, scales = "free_y")+
         xlab("codon position along ENST00000621411.2-H3C2")

save_plot(paste0(figurePrefix,"_plt_H3C2_psite_codon_sc.pdf"), plt_H3C2_sc, base_width = 6, base_height = 10)

plt_H3C2_pb <- H3C2 %>%
  group_by(codon_position,codon,three,hekMeta) %>%
  summarize(np = sum(n)) %>%
  ungroup() %>%
  group_by(hekMeta) %>%
  mutate(ptot = sum(np)) %>% 
  ungroup() %>%
  #mutate(codon_display = paste(codon,three,sep="_")) %>%
  mutate(codon_display = case_when( three %in% c("Arg","Leu") ~ paste(codon, three, sep="_"),
                                   TRUE ~ "other")) %>%
  mutate(cdl = as.numeric(factor(codon_display, 
                                 levels = c(grep("Arg",unique(paste(H3C2$codon,H3C2$three,sep="_")),value=TRUE),
                                            grep("Leu",unique(paste(H3C2$codon,H3C2$three,sep="_")),value=TRUE),
                                            "other"),ordered=TRUE))) %>% 
  mutate(codon_display = fct_reorder(codon_display, cdl)) %>%
  ggplot(aes(x=codon_position, y = np/ptot, fill=codon_display)) +
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c(brewer.pal(9,"Paired"),"grey50")) +
  facet_wrap(~hekMeta, ncol=1)+
  xlab("codon position along ENST00000621411.2-H3C2")
save_plot(paste0(figurePrefix, "_plt_H3C2_psite_codon_pseudobulk.pdf"), plt_H3C2_pb, base_width = 10, base_height = 10)


metaOut <- cbind(hekMeta,
                 mean_site_fc[rownames(hekMeta),])

write.csv(metaOut, file = paste0(figurePrefix, "_HEK293T_Starv_metadata.csv"))
write.csv(RPFcounts, file = paste0(figurePrefix, "_HEK293T_Starv_RPFcounts.csv"))

## Plate layouts
HEKStarvPL <- barcodes
HEKStarvPL$CB <- rownames(HEKStarvPL)
HEKStarvPL$fraction <- "Rich"
HEKStarvPL[ HEKStarvPL$well %in% columnFeatures(1:5,LETTERS[1:16]) , "fraction" ] = "Arg3h"
HEKStarvPL[ HEKStarvPL$well %in% columnFeatures(6:10,LETTERS[1:16]) , "fraction" ] = "Arg6h"
HEKStarvPL[ HEKStarvPL$well %in% columnFeatures(11:15,LETTERS[1:16]) , "fraction" ] = "Leu3h"
HEKStarvPL[ HEKStarvPL$well %in% columnFeatures(16:20,LETTERS[1:16]) , "fraction" ] = "Leu6h"
HEKStarvPL[ HEKStarvPL$well %in% columnFeatures(21:23,LETTERS[1:16]) , "fraction" ] = "Rich" 
HEKStarvPL[ HEKStarvPL$well %in% columnFeatures(24,LETTERS[1:16]), "fraction" ] = "NTC"
write.csv(HEKStarvPL, file = paste0(figurePrefix, "_HEK293T_Starv_layout.csv"), row.names = FALSE, quote = FALSE)
