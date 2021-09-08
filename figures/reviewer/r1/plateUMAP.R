#!/usr/bin/env Rscript
# Michael VanInsberghe 2021-04-23
# HEK plate layout and UMAP ~ plate figure
# Before running, make sure that the metadata tables have been updated with the UMAP and cluster results

library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(feather)
library(Matrix)
library(Matrix.utils)
library(Seurat)
library(ggplot2)
library(egg)
library(scales)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(forcats)
library(RColorBrewer)
library(uwot)
library(ggrepel)
library(ggseqlogo)
theme_set(theme_cowplot())

source('plotFunctions.R')

figurePrefix <- "plate/plate"

data_dir <- "../scRiboSeq/data"

hekMeta <- readRDS(file.path(data_dir, "HEK293T/compiled/clustered_meta.rds"))
hekAllMeta <- readRDS(file.path(data_dir, "HEK293T/compiled/allmeta.rds")) 
fucciMeta <- readRDS(file.path(data_dir, "Fucci/compiled/clustered_meta.rds"))
eeMeta <- readRDS(file.path(data_dir, "EE/compiled/clustered_meta.rds"))

plt_hek_plate <- hekAllMeta %>%
  filter(plate == "HEK293T-Starv1") %>%
  mutate(column = gsub("[A-Z]*","",well),
         row = gsub("[0-9]*","",well)) %>%
  mutate(column = fct_reorder(as.factor(column), rank(as.numeric(column))),
         row = fct_reorder(as.factor(row), -rank(row))) %>%
  ggplot(aes(x=column, y=row))+
  geom_tile(aes(fill=sort_population))+
  coord_equal()+
  #facet_wrap(~plate)+
  #scale_fill_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), na.value = "grey50")+
  scale_x_discrete(expand = c(0,0), position = "bottom")+
  scale_y_discrete(expand = c(0,0))+
  panel_border(colour="grey20")+
  theme(legend.position = "bottom",
        axis.text = element_text(size=8),
        axis.line = element_blank(),
        axis.title = element_blank())
save_plot(paste0(figurePrefix, "_plt_QC_plate.pdf"), plt_hek_plate, base_width = 4, base_height = 4)


plt_hek_umaps <- hekMeta %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, colour = plate))+
  geom_point(size = 0.25, data = hekMeta %>% select(UMAP_1, UMAP_2), colour = "grey60")+
  geom_point(size = 0.25)+
  coord_equal()+
  facet_wrap(~plate)+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

plt_fucci_umaps <- fucciMeta %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, colour = plate))+
  geom_point(size = 0.25, data = fucciMeta %>% select(UMAP_1, UMAP_2), colour = "grey60")+
  geom_point(size = 0.25)+
  coord_equal()+
  facet_wrap(~plate)+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

plt_ee_umaps <- eeMeta %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, colour = plate))+
  geom_point(size = 0.25, data = eeMeta %>% select(UMAP_1, UMAP_2), colour = "grey60")+
  geom_point(size = 0.25)+
  coord_equal()+
  facet_wrap(~plate, ncol = 2)+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
save_plot(paste0(figurePrefix, "_plate_ee_umaps.pdf"), plt_ee_umaps, base_width = 10, base_height = 10)


plt_plate_umaps <- plot_grid(plt_hek_umaps,
                             plt_fucci_umaps, 
                             plt_ee_umaps,
                             nrow = 3,
                             rel_heights = c(1,5,4),
                             rel_widths = c(2,6,2))

save_plot(paste0(figurePrefix, "_plate_umaps.pdf"), plt_plate_umaps, base_width = 10, base_height = 12)



fucFilMeta <- fucciMeta 
plt_fsc <- ggplot(fucFilMeta, aes(x=FSC.Area, y = FSC, colour = seurat_clusters))+
  geom_point(size=0.8)+
  coord_equal()
plt_ssc <- ggplot(fucFilMeta, aes(x=SSC.Area, y = SSC, colour = seurat_clusters))+
  geom_point(size=0.8)+
  coord_equal()
plt_ssc_w <- ggplot(fucFilMeta, aes(x=SSC.Width, y = SSC, colour = seurat_clusters))+
  geom_point(size=0.8)
plt_fsc_ssc <- ggplot(fucFilMeta, aes(x=FSC.Area, y = SSC.Area, colour = seurat_clusters))+
  geom_point(size=0.8)+
  coord_equal()
plt_fsc_dapi <- ggplot(fucFilMeta, aes(x=log10(X..405..460.50), y = FSC, colour = seurat_clusters))+
  geom_point(size=0.8)
save_plot(paste0(figurePrefix, "_Fucci_facs.pdf"), plot_grid(plt_fsc, plt_ssc, plt_ssc_w, plt_fsc_ssc, plt_fsc_dapi, ncol = 2, align = "hv"), base_width = 10, base_height = 8)



