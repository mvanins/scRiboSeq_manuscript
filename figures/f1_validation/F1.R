#!/usr/bin/env Rscript
# Michael VanInsberghe 2020-05-28
# Figure 1 & supporting scRibo-seq

library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(feather)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(forcats)
library(RColorBrewer)
theme_set(theme_cowplot())

source('plotFunctions.R')
figurePrefix <- "F1_updated/F1_updated"

data_dir <- "../scRiboSeq/data"

annot <- fread('../../data/hsa_annotations.csv.gz', stringsAsFactors = FALSE, data.table = FALSE)

canonical_pc <- annot %>%
  filter(set=="canonical", transcript_type == "protein_coding") %>%
  filter(chr != 'chrM') %>%
  filter(transcript_id != "ENST00000437139.7")


HEKReads <- readRDS(file.path(data_dir, "HEK293T/compiled/reads.rds"))
HEKMeta <- readRDS(file.path(data_dir, "HEK293T/compiled/meta.rds"))
HEKCounts <- readRDS(file.path(data_dir, "HEK293T/compiled/counts.rds"))

FucciReads <- readRDS(file.path(data_dir, "Fucci/compiled/reads.rds"))
FucciMeta <- readRDS(file.path(data_dir, "Fucci/compiled/meta.rds"))
FucciCounts <- readRDS(file.path(data_dir, "Fucci/compiled/counts.rds"))

## Wiggle Heatmaps

HEKwiggle <- bind_rows(wiggleRegion(HEKReads, 
                                    group = 'single',
                                    site = 'cut5',
                                    reference = 'cds_start', 
                                    mindist = -40,
                                    maxdist = 60,
                                    name = "start") %>% mutate(display_rank = dense_rank(csd)),
                       wiggleRegion(HEKReads,
                                    group = 'single',
                                    site = 'cut5',
                                    reference = 'cds_start', 
                                    mindist = 100,
                                    maxdist = 200,
                                    name = "middle") %>% mutate(display_rank = dense_rank(csd)+1E6),
                       wiggleRegion(HEKReads,
                                    group = 'single',
                                    site = 'cut5',
                                    reference = 'cds_end', 
                                    mindist = -60,
                                    maxdist = 40,
                                    name = "stop") %>% mutate(display_rank = dense_rank(csd)+1E7)) %>%
        group_by(CB) %>%
        mutate(cell_total = sum(n),
               cell_mean = mean(n),
               cell_median = median(n)) %>%
        ungroup() %>%
        inner_join(HEKMeta %>% rename(pc_total = cell_total), by = "CB") %>%
        mutate(fc = n/cell_mean,
               pct = n/pc_total) %>%
        group_by(CB) %>%
        select_top_groups(500, order_by = cell_total, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(display_position = dense_rank(display_rank)) %>%
        select(-display_rank)

HEKwiggleHM <- HEKwiggle %>%
  select(display_position, CB, fc) %>%
  spread(display_position, fc, fill = NA) %>%
  as.data.frame() %>%
  column_to_rownames(var = "CB")


Fucciwiggle <- bind_rows(wiggleRegion(FucciReads, 
                                      group = 'single',
                                      site = 'cut5',
                                      reference = 'cds_start', 
                                      mindist = -40,
                                      maxdist = 60,
                                      name = "start") %>% mutate(display_rank = dense_rank(csd)),
                         wiggleRegion(FucciReads,
                                      group = 'single',
                                      site = 'cut5',
                                      reference = 'cds_start', 
                                      mindist = 100,
                                      maxdist = 200,
                                      name = "middle") %>% mutate(display_rank = dense_rank(csd)+1E6),
                         wiggleRegion(FucciReads,
                                      group = 'single',
                                      site = 'cut5',
                                      reference = 'cds_end', 
                                      mindist = -60,
                                      maxdist = 40,
                                      name = "stop") %>% mutate(display_rank = dense_rank(csd)+1E7)) %>%
        group_by(CB) %>%
        mutate(cell_total = sum(n),
               cell_mean = mean(n),
               cell_median = median(n)) %>%
        ungroup() %>%
        inner_join(FucciMeta %>% rename(pc_total = cell_total), by = "CB") %>%
        mutate(fc = n/cell_mean,
               pct = n/pc_total) %>%
        group_by(CB) %>%
        select_top_groups(1000, order_by = cell_total, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(display_position = dense_rank(display_rank)) %>%
        select(-display_rank)

FucciwiggleHM <- Fucciwiggle %>%
  select(display_position, CB, fc) %>%
  spread(display_position, fc, fill = NA) %>%
  as.data.frame() %>%
  column_to_rownames(var = "CB")

wiggleHM <- rbind(HEKwiggleHM,
                  FucciwiggleHM)

pdf(paste0(figurePrefix, "_all_wigglehm.pdf"), height = 2.5, width = 7)
labelby = 20
length_annotation <- HeatmapAnnotation(pos = anno_mark(at = c(seq(1,101,by=labelby),
                                                              seq(102,202, by=labelby),
                                                              seq(203,303, by=labelby)),
                                                       labels = c(as.character(seq(-40,60,labelby)),
                                                                  as.character(seq(100,200,labelby)),
                                                                  as.character(seq(-60,40,labelby))),
                                                       which = "column", side = "bottom", labels_rot = 0, link_height = unit(3, "mm")))
Heatmap(as.matrix(wiggleHM),
        col = colorRamp2(seq(0, 10, length = 11), rev(brewer.pal(11, "RdYlBu"))), # log fold change
        bottom_annotation = length_annotation,
        na_col = "grey0",
        name = "5' cuts",
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_split = factor(c(rep("HEK293T", nrow(HEKwiggleHM)),rep("Fucci",nrow(FucciwiggleHM))), levels = c("HEK293T", "Fucci"), ordered = TRUE),
        column_split = factor(c(rep("start_codon", length(-40:60)), rep("cds", length(100:200)), rep("stop_codon", length(-40:60))), levels = c("start_codon", "cds", "stop_codon"), ordered = TRUE),
        column_title = NULL,
        border = TRUE,
        use_raster = TRUE,
        raster_device = "CairoPNG",
        raster_quality = 15)
dev.off()

### Frame heatmaps

HEKFrame5 <- HEKReads %>%
  mutate( region = case_when(cut5 < (l_utr5) ~ "utr5",
                             ((cut5  >= (l_utr5)) & (cut5 < ((l_utr5)+l_cds ) )) ~ "cds", 
                             (cut5 >= (l_utr5 + l_cds)) ~ "utr3" )) %>%
  filter(region == "cds") %>%
  group_by(CB,frame5) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(CB) %>%
  mutate(cell_tot = sum(n)) %>%
  ungroup() %>%
  mutate(frac = n/cell_tot,
         cell_rank = dense_rank(-cell_tot)) %>%
  mutate(group = "cut5",
         type = "HEK293T",
         plot = "cell") %>%
  rename(frame = frame5) %>%
  filter(cell_rank <= 200)

HEKFrameP <- HEKReads %>%
  mutate( region = case_when(psite < (l_utr5) ~ "utr5",
                             ((psite  >= (l_utr5)) & (cut5 < ((l_utr5)+l_cds ) )) ~ "cds", 
                             (psite >= (l_utr5 + l_cds)) ~ "utr3" )) %>%
  filter(region == "cds") %>%
  mutate(frameP = (psite - cds_start) %% 3) %>% 
  group_by(CB, frameP) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(CB) %>%
  mutate(cell_tot = sum(n)) %>%
  ungroup() %>%
  mutate(frac = n/cell_tot,
         cell_rank = dense_rank(-cell_tot)) %>%
  mutate(group = "psite",
         type= "HEK293T",
         plot = "cell") %>%
  rename(frame = frameP) %>%
  filter(cell_rank <= 200)

FucciFrame5 <- FucciReads %>%
  mutate( region = case_when(cut5 < (l_utr5) ~ "utr5",
                             ((cut5  >= (l_utr5)) & (cut5 < ((l_utr5)+l_cds ) )) ~ "cds", 
                             (cut5 >= (l_utr5 + l_cds)) ~ "utr3" )) %>%
  filter(region == "cds") %>%
  group_by(CB,frame5) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(CB) %>%
  mutate(cell_tot = sum(n)) %>%
  ungroup() %>%
  mutate(frac = n/cell_tot,
         cell_rank = dense_rank(-cell_tot)) %>%
  mutate(group = "cut5",
         type = "Fucci",
         plot = "cell") %>%
  rename(frame = frame5) %>%
  filter(cell_rank <= 200)

FucciFrameP <- FucciReads %>%
  mutate( region = case_when(psite < (l_utr5) ~ "utr5",
                             ((psite  >= (l_utr5)) & (cut5 < ((l_utr5)+l_cds ) )) ~ "cds", 
                             (psite >= (l_utr5 + l_cds)) ~ "utr3" )) %>%
  filter(region == "cds") %>%
  mutate(frameP = (psite - cds_start) %% 3) %>% 
  group_by(CB, frameP) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(CB) %>%
  mutate(cell_tot = sum(n)) %>%
  ungroup() %>%
  mutate(frac = n/cell_tot,
         cell_rank = dense_rank(-cell_tot)) %>%
  mutate(group = "psite",
         type = "Fucci",
         plot = "cell") %>%
  rename(frame = frameP) %>%
  filter(cell_rank <= 200)

FrameCB <- bind_rows(HEKFrame5, HEKFrameP,
                     FucciFrame5, FucciFrameP)


HEKFrame5L <- HEKReads %>%
  group_by(length,frame5) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(length) %>%
  mutate(length_tot = sum(n)) %>%
  ungroup() %>%
  mutate(frac = n/length_tot) %>%
  mutate(group = "cut5",
         type = "HEK293T",
         plot = "length") %>%
  rename(frame = frame5) 

HEKFramePL <- HEKReads %>%
  mutate(frameP = (psite - cds_start) %% 3) %>% 
  group_by(length,frameP) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(length) %>%
  mutate(length_tot = sum(n)) %>%
  ungroup() %>%
  mutate(frac = n/length_tot) %>%
  mutate(group = "psite",
         type = "HEK293T",
         plot = "length") %>%
  rename(frame = frameP) 

FucciFrame5L <- FucciReads %>%
  group_by(length,frame5) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(length) %>%
  mutate(length_tot = sum(n)) %>%
  ungroup() %>%
  mutate(frac = n/length_tot) %>%
  mutate(group = "cut5",
         type = "Fucci",
         plot = "length") %>%
  rename(frame = frame5) 

FucciFramePL <- FucciReads %>%
  mutate(frameP = (psite - cds_start) %% 3) %>% 
  group_by(length,frameP) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(length) %>%
  mutate(length_tot = sum(n)) %>%
  ungroup() %>%
  mutate(frac = n/length_tot) %>%
  mutate(group = "psite",
         type = "Fucci",
         plot = "length") %>%
  rename(frame = frameP) 

FrameL <- bind_rows(HEKFrame5L, HEKFramePL,
                    FucciFrame5L, FucciFramePL)


plt_frame_CB <- ggplot(FrameCB,aes(x=frame,y=CB))+
  geom_tile(aes(fill=frac)) +
  facet_wrap(~type+group, scales="free_y")+
  scale_fill_gradientn(colours=(brewer.pal(9,"Greys")),
                       limits = c( min(c(FrameL$frac,FrameCB$frac)), max(c(FrameL$frac,FrameCB$frac)))) + 
  scale_y_discrete(expand=c(0,0)) +
  scale_x_continuous(breaks=c(0,1,2),expand=c(0,0)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

save_plot(paste0(figurePrefix,"_frame_CB.pdf"),plt_frame_CB, base_width = 4, base_height = 6)


plt_frame_L <- ggplot(FrameL,aes(x=frame,y=length))+
  geom_tile(aes(fill=frac)) +
  facet_wrap(~type+group, scales="free_y")+
  scale_fill_gradientn(colours=(brewer.pal(9,"Greys")),
                       limits = c( min(c(FrameL$frac,FrameCB$frac)), max(c(FrameL$frac,FrameCB$frac)))) + 
  scale_y_continuous(breaks=full_seq(min(FrameL$length):max(FrameL$length),1),expand=c(0,0))+
  scale_x_continuous(breaks=c(0,1,2),expand=c(0,0))+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1))

save_plot(paste0(figurePrefix,"_frame_length.pdf"), plt_frame_L, base_width = 4, base_height = 6)


# frame per length plot with marginal distribution on relative read number on Psite plot

plt_fucci_Lp <- ggplot(FrameL %>% filter(type == "Fucci", group == "psite"), aes(x = frame, y = length)) +
  geom_tile(aes(fill=frac)) +
  scale_fill_gradientn(colours=(brewer.pal(9,"Greys")),
                       limits = c( min(c(FrameL$frac,FrameCB$frac)), max(c(FrameL$frac,FrameCB$frac)))) + 
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks=c(0,1,2),expand=c(0,0))+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1))

plt_fucci_length_hist <- axis_canvas(plt_fucci_Lp, axis = 'y', coord_flip = TRUE)+
  geom_bar(data = FrameL %>% filter(type == "Fucci", group == "psite", frame == 0), aes(x=length, y=length_tot), stat="identity")+
  coord_flip()

plt_fucci_length_combined <- ggdraw(insert_yaxis_grob(plt_fucci_Lp, plt_fucci_length_hist, position = "right", width = grid::unit(.6, "null")))



plt_HEK_Lp <- ggplot(FrameL %>% filter(type == "HEK293T", group == "psite"), aes(x = frame, y = length)) +
  geom_tile(aes(fill=frac)) +
  scale_fill_gradientn(colours=(brewer.pal(9,"Greys")),
                       limits = c( min(c(FrameL$frac,FrameCB$frac)), max(c(FrameL$frac,FrameCB$frac)))) + 
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks=c(0,1,2),expand=c(0,0))+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1))

plt_HEK_length_hist <- axis_canvas(plt_HEK_Lp, axis = 'y', coord_flip = TRUE)+
  geom_bar(data = FrameL %>% filter(type == "HEK293T", group == "psite", frame == 0), aes(x=length, y=length_tot), stat="identity")+
  coord_flip()

plt_HEK_length_combined <- ggdraw(insert_yaxis_grob(plt_HEK_Lp, plt_HEK_length_hist, position = "right", width = grid::unit(.6, "null")))



# add axis to marginal histogram

plt_HEK_length_hist2 <- ggplot(FrameL %>% filter(type == "HEK293T", group == "psite", frame == "0"), 
                               aes(x=length, y = 100*length_tot/sum(length_tot)))+
        geom_bar(stat="identity") + 
        coord_flip() +
        scale_y_continuous(breaks = c(0,5,10), expand=c(0,0))+
        scale_x_continuous(expand=c(0,0))+
        theme(axis.text.y = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank())+
        ylab("% of reads")

plt_HEK_length_combined_label <- plot_grid(plt_HEK_Lp, plt_HEK_length_hist2, align = "h", ncol = 2, rel_widths = c(1,0.4))

plt_fucci_length_hist2 <- ggplot(FrameL %>% filter(type == "Fucci", group == "psite", frame == "0"), 
                               aes(x=length, y = 100*length_tot/sum(length_tot)))+
        geom_bar(stat="identity") + 
        coord_flip() +
        scale_y_continuous(breaks = c(0,5,10), expand=c(0,0))+
        scale_x_continuous(expand=c(0,0))+
        theme(axis.text.y = element_blank(),
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank())+
        ylab("% of reads")

plt_fucci_length_combined_label <- plot_grid(plt_fucci_Lp, plt_fucci_length_hist2, align = "h", ncol = 2, rel_widths = c(1,0.4))
save_plot(paste0(figurePrefix, "_frame_length_psite_lengthdist_label.pdf"), plot_grid(plt_fucci_length_combined_label, plt_HEK_length_combined_label, ncol = 1), base_width = 3, base_height = 6)

## post-correction wiggles


fucciCut5 <- wiggleRegion(FucciReads,
                          group = "single",
                          site = "cut5",
                          reference = "cds_start",
                          mindist = 100,
                          maxdist = 200,
                          name = "cut5") %>%
        group_by(CB) %>%
        mutate(cell_total = sum(n)) %>%
        ungroup() %>%
        group_by(CB) %>% 
        select_top_groups(100, order_by = cell_total, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(CB = fct_reorder(CB, cell_total))
        

fucciPsite <- wiggleRegion(FucciReads,
                           group = "single",
                           site = "psite",
                           reference = "cds_start",
                           mindist = 100,
                           maxdist = 200,
                           name = "psite") %>%
        right_join(fucciCut5 %>% select(CB, cut5_cell_total = cell_total), 
                   by = "CB") %>%
        mutate(CB = fct_reorder(CB, cut5_cell_total))
    
fucciCut5av <- fucciCut5 %>%
  group_by(csd) %>%
  summarize(average = mean(n)) %>%
  ungroup()
      
fucciPsiteav <- fucciPsite %>%
  group_by(csd) %>%
  summarize(average = mean(n)) %>%
  ungroup() 


plt_predict_cut5_hm <- ggplot(fucciCut5, aes(x=csd, y=CB, fill = n))+
  geom_tile() + 
  scale_fill_gradientn(breaks = seq(0, max(c(fucciCut5$n, fucciPsite$n)), 50),
                       limits = c(0, max(c(fucciCut5$n, fucciPsite$n))),
                       colours=rev(brewer.pal(10,"RdYlBu")),na.value=rev(brewer.pal(10,"RdYlBu"))[1])+
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "bottom")+
  xlab("Distance to start codon [nt]")

plt_predict_psite_hm <- ggplot(fucciPsite, aes(x=csd, y=CB, fill = n))+
  geom_tile() + 
  scale_fill_gradientn(breaks = seq(0, max(c(fucciCut5$n, fucciPsite$n)), 50),
                       limits = c(0, max(c(fucciCut5$n, fucciPsite$n))),
                       colours=rev(brewer.pal(10,"RdYlBu")),na.value=rev(brewer.pal(10,"RdYlBu"))[1])+
  scale_x_continuous(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "bottom")+
  xlab("Distance to start codon [nt]")

  
plt_predict_cut5_bar <- ggplot(fucciCut5av, aes(x = csd, y = average)) +
  geom_bar(stat = "identity") + 
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, max(c(fucciCut5av$average, fucciPsiteav$average)), 10),
                     limits = c(0, max(c(fucciCut5av$average, fucciPsiteav$average)))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())+
  ylab("Average cuts per cell")
 
plt_predict_psite_bar <- ggplot(fucciPsiteav, aes(x = csd, y = average)) +
  geom_bar(stat = "identity") + 
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0),
                     breaks = seq(0, max(c(fucciCut5av$average, fucciPsiteav$average)), 10),
                     limits = c(0, max(c(fucciCut5av$average, fucciPsiteav$average)))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank())

save_plot(paste0(figurePrefix, "_corrected_wiggle.pdf"), plot_grid(plot_grid(plt_predict_cut5_bar, plt_predict_cut5_hm, ncol = 1, align = "v", rel_heights = c(0.6, 1)),
                                                                   plot_grid(plt_predict_psite_bar, plt_predict_psite_hm, ncol = 1, align = "v", rel_heights = c(0.6, 1)), ncol = 2),
          base_width = 8, base_height = 4)

#### Region coverage

expanded_pc <- canonical_pc %>%
  mutate(l_utr5 = l_utr5 - 25,
         l_cds = l_cds + 25,
         l_utr3 = as.numeric(l_utr3) ) %>%
  mutate(l_utr5 = case_when( l_utr5 < 0 ~ 0,
                            l_utr5 >= 0 ~ l_utr5),
         l_cds = case_when(l_cds < 0 ~ 0,
                           l_cds >= 0 ~ l_cds),
         l_utr3 = case_when(l_utr3 < 0 ~ 0,
                            l_utr3 >= 0 ~ l_utr3))

pcReads <- bind_rows(HEKReads %>% mutate(CB = paste0("HEK293T-", CB)),
                     FucciReads %>% mutate(CB = paste0("Fucci-", CB))) %>%
  mutate(l_utr5 = l_utr5 - 25,
         l_cds = l_cds + 25,
         l_utr3 = as.numeric(l_utr3) ) %>%
  mutate(l_utr5 = case_when( l_utr5 < 0 ~ 0,
                            l_utr5 >= 0 ~ l_utr5),
         l_cds = case_when(l_cds < 0 ~ 0,
                           l_cds >= 0 ~ l_cds),
         l_utr3 = case_when(l_utr3 < 0 ~ 0,
                            l_utr3 >= 0 ~ l_utr3))
  

all_on_one <- pcReads %>%
  group_by(CB, transcript_id, l_utr5, l_cds, l_utr3) %>%
  summarize(sc_tran_n = n()) %>%
  ungroup() %>% 
  filter(sc_tran_n > 0) %>%
  group_by(CB) %>%
  summarize(bp_utr5 = sum(l_utr5),
            bp_cds = sum(l_cds),
            bp_utr3 = sum(l_utr3)) %>%
  ungroup()

all_separate <- pcReads %>%
  group_by(CB, transcript_id, l_utr5, l_cds, l_utr3) %>%
  mutate(sc_tran_n = n()) %>%
  ungroup() %>% 
  filter(sc_tran_n > 0) %>%
  group_by(CB) %>%
  summarize(bp_utr5_sep = sum(l_utr5),
            bp_cds_sep = sum(l_cds),
            bp_utr3_sep = sum(l_utr3)) %>%
  ungroup()



tpm_pc <- pcReads %>%
  mutate( region = case_when(cut5 < l_utr5 ~"utr5",
                             ( (cut5 >= l_utr5) & (cut5 < (l_utr5+l_cds)) ) ~ "cds",
                             ( cut5 >= (l_utr5 + l_cds) ) ~ "utr3" ) ) %>%
  group_by(CB, region) %>%
  summarize(reads_per_region = n()) %>%
  ungroup() %>%
  group_by(CB) %>%
  mutate(reads_per_cell = sum(reads_per_region)) %>%
  ungroup() %>%
  inner_join(all_on_one, by = "CB") %>%
  inner_join(all_separate, by = "CB") %>%
  mutate(BP_one = case_when(region == 'utr5' ~ bp_utr5,
                            region == 'cds' ~ bp_cds,
                            region == 'utr3' ~ bp_utr3),
         BP_sep = case_when(region == 'utr5' ~ bp_utr5_sep,
                            region == 'cds' ~ bp_cds_sep,
                            region == 'utr3' ~ bp_utr3_sep) ) %>%
  select(-starts_with('bp_', ignore.case = FALSE)) %>%
  mutate(rpk_one = reads_per_region / BP_one,
         rpk_sep = reads_per_region / BP_sep) %>% 
  group_by(CB) %>%
  mutate(norm_one = sum(rpk_one)/100,
         norm_sep = sum(rpk_sep)/100) %>%
  ungroup() %>%
  mutate(tpm_one = rpk_one/norm_one,
         tpm_sep = rpk_sep/norm_sep) %>%
  mutate(regionord = case_when(region == "utr5" ~ 0,
                               region == "cds" ~ 1,
                               region == "utr3" ~ 2)) %>%
  mutate(region = fct_reorder(region,regionord)) %>%
  separate(CB, c("sample",NA), remove = FALSE, sep="_") %>%
  separate(sample, c("group", NA), remove = FALSE, sep = "-") %>%
  mutate(groupord = case_when(group == "Fucci" ~ 1,
                              group == "HEK293T" ~ 0,
                              group == "darnell" ~ 2,
                              group == "martinez" ~ 3,
                              group == "tanenbaum" ~ 4)) %>%
  mutate(group = fct_reorder(group,groupord))

plt_tpm_one <- tpm_pc %>%
  ggplot(aes(x=region, y=tpm_one, fill=group)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#66C2A5", "#FC8D62"))+
  ggtitle("all footprints on one transcript")+
  ylab("5' cuts per region length per hundred")

save_plot(paste0(figurePrefix, "_UTR_tpm_allonone.pdf"), plt_tpm_one, base_width = 4, base_height = 5)

plt_tpm_sep <- tpm_pc %>%
  ggplot( aes(x=region, y=tpm_sep, fill=group)) +
  geom_boxplot()+
  scale_fill_manual(values = c("#1f78b4", "#33a02c"))+
  ggtitle("one transcript per footprint")+
  ylab("5' cuts per region length per hundred")

#save_plot(paste0(figurePrefix,"_UTR_tpm_SClength.pdf"), plot_grid(plt_tpm_one, plt_tpm_sep, nrow=1), base_width = 12, base_height = 8)

plt_tpm_one_separate <- tpm_pc %>%
  ggplot(aes(x=region, y=tpm_one))+
  geom_boxplot()+
  facet_wrap(~group)+
  ylab("5' cuts per region length per hundred")


#save_plot(paste0(figurePrefix, "_UTR_tpm_wrapped.pdf"), plt_tpm_one_separate, base_width = 3, base_height = 6)


### Numbers
# pct of reads in frame
cut5_numbers <- FrameCB %>%
  group_by(group, frame) %>%
  summarize(mean_frame = 100*mean(frac),
            sd_frame = sd(100*frac),
            sem_frame = sd(100*frac)/sqrt(n())) %>%
  ungroup()


### Sequence logo
suppressPackageStartupMessages(library(ggseqlogo))

seqdf <- rbind(HEKReads, FucciReads) %>%
  mutate(motif5 = paste0(base5m8,base5m7,base5m6,base5m5,base5m4,base5m3,base5m2,base5m1,
                       base5,
                       base5p1,base5p2,base5p3,base5p4,base5p5,base5p6,base5p7,base5p8),
         motif3 = paste0(base3m8,base3m7,base3m6,base3m5,base3m4,base3m3,base3m2,base3m1,
                       base3,
                       base3p1,base3p2,base3p3,base3p4,base3p5,base3p6,base3p7,base3p8)) %>%
  select(motif5, motif3) %>%
  mutate(motif5 = gsub("T", "U", motif5),
         motif3 = gsub("T", "U", motif3)) %>%
  filter(nchar(motif5) == 17,
         nchar(motif3) == 17)


plt_icon5_p <- ggplot() + 
  geom_logo(seqdf$motif5,method='prob')+
  scale_x_continuous(breaks=1:17,labels = -8:8)+
  xlab("distance from 5' base")

plt_icon3_p <- ggplot() +
  geom_logo(seqdf$motif3,method='prob')+
  scale_x_continuous(breaks=1:17,labels = -8:8)+
  xlab("distance from 3' base")

plt_icon_figure <- plot_grid(plt_icon5_p, plt_icon3_p, ncol = 2)

save_plot(paste0(figurePrefix, "_sequence_icon.pdf"), plt_icon_figure, base_width = 9, base_height = 2)


### export for GEO
HEK_only <- HEKMeta %>% filter(plate == "HEK293T") %>% select(CB) %>% pull()

# barcodes
barcodes <- read.table('/hpc/hub_oudenaarden/mvanins/local/lib/barcodes/miRv6_barcodemap.tsv',
                       sep="\t",
                       header = FALSE,
                       stringsAsFactors = FALSE,
                       row.names = 1,
                       col.names = c("","number","well"))

treatment <- data.frame(row.names = colnames(HEKCounts),
                        treatment = rep("Rich", ncol(HEKCounts)),
                        well = barcodes[gsub(".*_(.*?)$","\\1",colnames(HEKCounts),perl=TRUE),"well"],
                        plate = gsub("^(.*?)[_].*","\\1",colnames(HEKCounts),perl=TRUE),
                        stringsAsFactors = FALSE)

HEKmeta <- cbind(HEKStats[HEK_only,],
                 treatment[HEK_only,])


# annot
annot <- fread('../../data/hsa_annotations.csv.gz', stringsAsFactors = FALSE, data.table = FALSE)

canonical_pc <- annot %>%
  filter(set=="canonical", transcript_type == "protein_coding") %>%
  filter(chr != 'chrM') %>%
  filter(transcript_id != "ENST00000437139.7")


ENSG_to_name <- data.frame(count_names = rownames(HEKCounts), stringsAsFactors = FALSE) %>%
		inner_join(annot, by = c("count_names" = "gene_id")) %>% 
		mutate(gene_out = paste(count_names, gene_name, sep = "-")) %>%
		select(count_names, gene_out) %>%
		distinct() %>%
		column_to_rownames(var = "count_names")
RPFcounts <- HEKCounts[!(rownames(HEKCounts) %in% c("cds", "cds_frac", "cell_total", "utr3", "utr5")),]
rownames(RPFcounts) <- ENSG_to_name[rownames(RPFcounts),"gene_out"]
RPFcounts <- RPFcounts[,HEK_only]
RPFcounts <- RPFcounts[apply(RPFcounts>0,1,sum)>0,]

write.csv(HEKMeta %>% filter(CB %in% HEK_only), file = paste0(figurePrefix, "_HEK293T_metadata.csv"))
write.csv(as.matrix(RPFcounts[,HEK_only]), file = paste0(figurePrefix, "_HEK293T_RPFcounts.csv"))

## Plate layouts
HEKonlyPL <- barcodes
HEKonlyPL$CB <- rownames(barcodes)
HEKonlyPL$fraction <- "HEK293T"
HEKonlyPL[ HEKonlyPL$well %in% columnFeatures(24,LETTERS[1:16]), "fraction" ] = "NTC"
write.csv(HEKonlyPL, file = paste0(figurePrefix, "_HEK293T_layout.csv"), row.names = FALSE, quote = FALSE)
