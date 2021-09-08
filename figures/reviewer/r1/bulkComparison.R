#!/usr/bin/env Rscript
# Michael VanInsberghe 2021-03-13
# Compare to standard bulk RPF

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

figurePrefix <- "bulk/bulk_comparison"

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


## bulk data
bulkReads <- readRDS(file.path(data_dir, "bulk/compiled/reads.rds")) %>% select(-filename)
bulkMeta <- readRDS(file.path(data_dir, "bulk/compiled/meta.rds"))
bulkCounts <- readRDS(file.path(data_dir, "bulk/compiled/counts.rds"))
bulkBiotypes <- readRDS(file.path(data_dir, "bulk/compiled/biotypes.rds"))

## scRPF data
hekReads <- readRDS(file.path(data_dir, "HEK293T/compiled/reads.rds"))
hekMeta <- readRDS(file.path(data_dir, "HEK293T/compiled/meta.rds"))
hekCounts <- readRDS(file.path(data_dir, "HEK293T/compiled/counts.rds"))
hekBiotypes <- readRDS(file.path(data_dir, "HEK293T/compiled/biotypes.rds"))

fucciReads <- readRDS(file.path(data_dir, "Fucci/compiled/reads.rds"))
fucciMeta <- readRDS(file.path(data_dir, "Fucci/compiled/meta.rds"))
fucciCounts <- readRDS(file.path(data_dir, "Fucci/compiled/counts.rds"))
fucciBiotypes <- readRDS(file.path(data_dir, "Fucci/compiled/biotypes.rds"))

eeReads <- readRDS(file.path(data_dir, "EE/compiled/reads.rds"))
eeMeta <- readRDS(file.path(data_dir, "EE/compiled/meta.rds"))
eeCounts <- readRDS(file.path(data_dir, "EE/compiled/counts.rds"))
eeBiotypes <- readRDS(file.path(data_dir, "EE/compiled/biotypes.rds"))



## Wiggle plots

pseudobulkWiggle <- function(reads, site = "cut5", sample_total, sample_name){
  wiggle <- bind_rows(wiggleRegion(reads,
                                   group = 'bulk',
                                   site = site,
                                   reference = 'cds_start',
                                   mindist = -40,
                                   maxdist = 60,
                                   name = "start") %>% mutate(display_rank = dense_rank(csd)),
                      wiggleRegion(reads,
                                   group = 'bulk',
                                   site = site,
                                   reference = 'cds_start', 
                                   mindist = 100,
                                   maxdist = 200,
                                   name = "middle") %>% mutate(display_rank = dense_rank(csd)+1E6),
                      wiggleRegion(reads,
                                   group = 'bulk',
                                   site = site,
                                   reference = 'cds_end', 
                                   mindist = -60,
                                   maxdist = 40,
                                   name = "stop") %>% mutate(display_rank = dense_rank(csd)+2E6)) %>%
    mutate(library_mean = mean(n),
           pct = n/sample_total,
           fc = n/library_mean) %>%
    mutate(display_position = dense_rank(display_rank)) %>%
    select(-display_rank) %>%
    mutate(CB = sample_name)
  hm <- wiggle %>%
    select(display_position, CB, pct) %>%
    spread(display_position, pct, fill = NA) %>%
    as.data.frame() %>%
    column_to_rownames(var = "CB")
  return(hm)
}


scRiboWiggle5 <- rbind(pseudobulkWiggle(hekReads %>% filter(CB %in% pull(hekMeta %>% filter(sort_population == "Rich") %>% select(CB))),
                                        site = "cut5",
                                        sample_total = sum(hekMeta %>% filter(sort_population == "Rich") %>% select(cell_total)),
                                        sample_name = "HEK293T_Rich"),
                       pseudobulkWiggle(hekReads %>% filter(CB %in% pull(hekMeta %>% filter(grepl("^Arg", sort_population)) %>% select(CB))),
                                        site = "cut5",
                                        sample_total = sum(hekMeta %>% filter(grepl("^Arg", sort_population)) %>% select(cell_total)),
                                        sample_name = "HEK293T_Arg"),
                       pseudobulkWiggle(hekReads %>% filter(CB %in% pull(hekMeta %>% filter(grepl("^Leu", sort_population)) %>% select(CB))),
                                        site = "cut5",
                                        sample_total = sum(hekMeta %>% filter(grepl("^Leu", sort_population)) %>% select(cell_total)),
                                        sample_name = "HEK293T_Leu"),
                       pseudobulkWiggle(fucciReads %>% filter(CB %in% pull(fucciMeta %>% filter(sort_population %in% c("Interphase", "cycling_pre", "cycling_post")) %>% select(CB))),
                                        site = "cut5", 
                                        sample_total = sum(fucciMeta %>% filter(sort_population %in% c("Interphase", "cycling_pre", "cycling_post")) %>% select(cell_total)),
                                        sample_name = "Fucci_Interphase"),
                       pseudobulkWiggle(fucciReads %>% filter(CB %in% pull(fucciMeta %>% filter(sort_population %in% c("Mit")) %>% select(CB))),
                                        site = "cut5",
                                        sample_total = sum(fucciMeta %>% filter(sort_population %in% c("Mit")) %>% select(cell_total)),
                                        sample_name = "Fucci_Mit"),
                       pseudobulkWiggle(fucciReads %>% filter(CB %in% pull(fucciMeta %>% filter(sort_population %in% c("G0")) %>% select(CB))),
                                        site = "cut5",
                                        sample_total = sum(fucciMeta %>% filter(sort_population %in% c("G0")) %>% select(cell_total)),
                                        sample_name = "Fucci_G0"),
                       pseudobulkWiggle(eeReads,
                                        site = "cut5",
                                        sample_total = sum(eeMeta %>% select(cell_total)),
                                        sample_name = "EE"))


scRiboWiggleP <- rbind(pseudobulkWiggle(hekReads %>% filter(CB %in% pull(hekMeta %>% filter(sort_population == "Rich") %>% select(CB))),
                                        site = "psite",
                                        sample_total = sum(hekMeta %>% filter(sort_population == "Rich") %>% select(cell_total)),
                                        sample_name = "HEK293T_Rich"),
                       pseudobulkWiggle(hekReads %>% filter(CB %in% pull(hekMeta %>% filter(grepl("^Arg", sort_population)) %>% select(CB))),
                                        site = "psite",
                                        sample_total = sum(hekMeta %>% filter(grepl("^Arg", sort_population)) %>% select(cell_total)),
                                        sample_name = "HEK293T_Arg"),
                       pseudobulkWiggle(hekReads %>% filter(CB %in% pull(hekMeta %>% filter(grepl("^Leu", sort_population)) %>% select(CB))),
                                        site = "psite",
                                        sample_total = sum(hekMeta %>% filter(grepl("^Leu", sort_population)) %>% select(cell_total)),
                                        sample_name = "HEK293T_Leu"),
                       pseudobulkWiggle(fucciReads %>% filter(CB %in% pull(fucciMeta %>% filter(sort_population %in% c("Interphase", "cycling_pre", "cycling_post")) %>% select(CB))),
                                        site = "psite", 
                                        sample_total = sum(fucciMeta %>% filter(sort_population %in% c("Interphase", "cycling_pre", "cycling_post")) %>% select(cell_total)),
                                        sample_name = "Fucci_Interphase"),
                       pseudobulkWiggle(fucciReads %>% filter(CB %in% pull(fucciMeta %>% filter(sort_population %in% c("Mit")) %>% select(CB))),
                                        site = "psite",
                                        sample_total = sum(fucciMeta %>% filter(sort_population %in% c("Mit")) %>% select(cell_total)),
                                        sample_name = "Fucci_Mit"),
                       pseudobulkWiggle(fucciReads %>% filter(CB %in% pull(fucciMeta %>% filter(sort_population %in% c("G0")) %>% select(CB))),
                                        site = "psite",
                                        sample_total = sum(fucciMeta %>% filter(sort_population %in% c("G0")) %>% select(cell_total)),
                                        sample_name = "Fucci_G0"),
                       pseudobulkWiggle(eeReads,
                                        site = "psite",
                                        sample_total = sum(eeMeta %>% select(cell_total)),
                                        sample_name = "EE"))


bulkWiggle <- bind_rows(wiggleRegion(bulkReads,
                                     group = "single",
                                     site = "cut5",
                                     reference = 'cds_start',
                                     mindist = -40,
                                     maxdist = 60,
                                     name = "start") %>% mutate(display_rank = dense_rank(csd)),
                        wiggleRegion(bulkReads,
                                     group = 'single',
                                     site = "cut5",
                                     reference = 'cds_start', 
                                     mindist = 100,
                                     maxdist = 200,
                                     name = "middle") %>% mutate(display_rank = dense_rank(csd)+1E6),
                        wiggleRegion(bulkReads,
                                     group = 'single',
                                     site = "cut5",
                                     reference = 'cds_end', 
                                     mindist = -60,
                                     maxdist = 40,
                                     name = "stop") %>% mutate(display_rank = dense_rank(csd)+2E6)) %>%
  group_by(CB) %>%
  mutate(library_mean = mean(n)) %>%
  ungroup() %>%
  inner_join(bulkMeta %>% select(CB, cell_total), by = "CB") %>%
  mutate(fc = n/library_mean,
         pct = n/cell_total,
         display_position = dense_rank(display_rank)) %>%
  select(-display_rank) %>%
  select(display_position, CB, pct) %>%
  spread(display_position, pct, fill = NA) %>%
  as.data.frame() %>%
  column_to_rownames(var = "CB")

darnellWiggleP <- bind_rows(wiggleRegion(bulkReads %>% filter(sample == "darnell"),
                                     group = "single",
                                     site = "psite",
                                     reference = 'cds_start',
                                     mindist = -40,
                                     maxdist = 60,
                                     name = "start") %>% mutate(display_rank = dense_rank(csd)),
                        wiggleRegion(bulkReads %>% filter(sample == "darnell"),
                                     group = 'single',
                                     site = "psite",
                                     reference = 'cds_start', 
                                     mindist = 100,
                                     maxdist = 200,
                                     name = "middle") %>% mutate(display_rank = dense_rank(csd)+1E6),
                        wiggleRegion(bulkReads %>% filter(sample == "darnell"),
                                     group = 'single',
                                     site = "psite",
                                     reference = 'cds_end', 
                                     mindist = -60,
                                     maxdist = 40,
                                     name = "stop") %>% mutate(display_rank = dense_rank(csd)+2E6)) %>%
  group_by(CB) %>%
  mutate(library_mean = mean(n)) %>%
  ungroup() %>%
  inner_join(bulkMeta %>% select(CB, cell_total), by = "CB") %>%
  mutate(fc = n/library_mean,
         pct = n/cell_total,
         display_position = dense_rank(display_rank)) %>%
  select(-display_rank) %>%
  select(display_position, CB, pct) %>%
  spread(display_position, pct, fill = NA) %>%
  as.data.frame() %>%
  column_to_rownames(var = "CB")


rownames(bulkWiggle) <- gsub("_SRR[0-9]+$","", rownames(bulkWiggle))
rownames(darnellWiggleP) <- gsub("_SRR[0-9]+$","", rownames(darnellWiggleP))

bulkWiggle <- bulkWiggle[c("darnell_Rich3h", "darnell_Rich6h", "darnell_Arg3h-r1", "darnell_Arg6h-r1", "darnell_Arg6h-r2", "darnell_Leu3h", "darnell_Leu6h",
                           "ingolia_HEK293T-low", "ingolia_HEK293T-med", "ingolia_HEK293T-high", 
                           "martinez_293T-lores", "martinez_293T-medres", "martinez_293T-hires", 
                           "tanenbaum_M-r1", "tanenbaum_M-r2", "tanenbaum_MG1-r1", "tanenbaum_MG1-r2", "tanenbaum_G2-r1", "tanenbaum_G2-r2"),]
darnellWiggleP <- darnellWiggleP[c("darnell_Rich3h", "darnell_Rich6h", "darnell_Arg3h-r1", "darnell_Arg6h-r1", "darnell_Arg6h-r2", "darnell_Leu3h", "darnell_Leu6h"),]

rownames(scRiboWiggleP) <- paste(rownames(scRiboWiggleP), "psite", sep = "_")
rownames(darnellWiggleP) <- paste(rownames(darnellWiggleP), "psite", sep = "_")


wiggleHM <- 100*rbind(scRiboWiggle5,
                      bulkWiggle,
                      scRiboWiggleP,
                      darnellWiggleP)

hm_quantiles <- 0.01
pdf(paste0(figurePrefix, "_pseudobulk_wigglehm.pdf"), height = 6, width = 10)
labelby = 20
length_annotation <- HeatmapAnnotation(pos = anno_mark(at = c(seq(1,101,by=labelby),
                                                              seq(102,202, by=labelby),
                                                              seq(203,303, by=labelby)),
                                                       labels = c(as.character(seq(-40,60,labelby)),
                                                                  as.character(seq(100,200,labelby)),
                                                                  as.character(seq(-60,40,labelby))),
                                                       which = "column", side = "bottom", labels_rot = 0, link_height = unit(3, "mm")))
Heatmap(as.matrix(wiggleHM),
        col = colorRamp2(
                         seq(quantile(as.matrix(wiggleHM), probs = hm_quantiles), quantile(as.matrix(wiggleHM), probs = 1-hm_quantiles), length = 11),
                         rev(brewer.pal(11, "RdYlBu"))),
        bottom_annotation = length_annotation,
        na_col = "grey0",
        name = 'pct of reads',
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        row_split = factor(c(rep("scRibo-cut5", nrow(scRiboWiggle5)),
                             rep("darnell", sum(grepl("^darnell", rownames(bulkWiggle)))),
                             rep("ingolia", sum(grepl("^ingolia", rownames(bulkWiggle)))),
                             rep("martinez", sum(grepl("^martin", rownames(bulkWiggle)))),
                             rep("tanenbaum", sum(grepl("^tanenbaum", rownames(bulkWiggle)))),
                             rep("scRibo-psite", nrow(scRiboWiggleP)),
                             rep("darnell-psite", nrow(darnellWiggleP))),
                           levels = c("scRibo-cut5", "darnell", "ingolia", "martinez", "tanenbaum", "scRibo-psite", "darnell-psite"),
                           ordered = TRUE),
        row_title_rot = 0,
        column_split = factor(c(rep("start_codon", length(-40:60)), rep("cds", length(100:200)), rep("stop_codon", length(-40:60))), levels = c("start_codon", "cds", "stop_codon"), ordered = TRUE),
        column_title = NULL,
        border = TRUE,
        use_raster = TRUE,
        raster_device = "CairoPNG",
        raster_quality = 15)
dev.off()



## Region mapping

expanded_canonical <- rbind(canonical_pc,
                            mm_canonical_pc) %>%
                              mutate(lr_utr5 = l_utr5 - 25,
                                     lr_cds = l_cds + 25,
                                     lr_utr3 = as.numeric(l_utr3)) %>%
                              mutate(lr_utr5 = case_when(lr_utr5 < 0 ~ 0,
                                                         lr_utr5 >= 0 ~ lr_utr5),
                                     lr_cds = case_when(lr_cds < 0 ~ 0,
                                                        lr_cds >= 0 ~ lr_cds),
                                     lr_utr3 = case_when(lr_utr3 < 0 ~ 0,
                                                         lr_utr3 >= 0 ~ lr_utr3)) %>%
                              select(transcript_id, starts_with("lr_"))

## count bp per region for each library
# either by assuming all RPFs originate from the same transcript,
# or that they originate from separate transcripts

regionBases <- rbind(hekReads,
                     fucciReads,
                     eeReads,
                     bulkReads) %>%
  select(CB, transcript_id) %>%
  inner_join(expanded_canonical, by = "transcript_id") %>%
  group_by(CB) %>%
  summarize(utr5 = sum(lr_utr5),
            cds = sum(lr_cds),
            utr3 = sum(lr_utr3)) %>%
  ungroup() %>%
  pivot_longer(cols = c("utr5", "cds", "utr3"), names_to = "region", values_to = "bp_sep") %>%
  inner_join(rbind(hekReads,
                   fucciReads,
                   eeReads,
                   bulkReads) %>%
               select(CB, transcript_id) %>%
               distinct() %>%
               inner_join(expanded_canonical, by = "transcript_id") %>%
               group_by(CB) %>%
               summarize(utr5 = sum(lr_utr5),
                         cds = sum(lr_cds),
                         utr3 = sum(lr_utr3)) %>%
               ungroup() %>%
               pivot_longer(cols = c("utr5", "cds", "utr3"), names_to = "region", values_to = "bp_one"),
             by = c("CB" = "CB", "region" = "region"))



regionTPM <- rbind(hekReads,
                   fucciReads,
                   eeReads,
                   bulkReads) %>%
  select(CB, transcript_id, cut5) %>%
  inner_join(expanded_canonical, by = "transcript_id") %>%
  mutate(read_region = case_when(cut5 < lr_utr5 ~ "utr5",
                                 ( (cut5 >= lr_utr5) & (cut5 < (lr_utr5 + lr_cds)) ) ~ "cds",
                                 ( cut5 >= (lr_utr5 + lr_cds)) ~ "utr3")) %>%
  group_by(CB, read_region) %>%
  summarize(reads_per_region = n()) %>%
  ungroup() %>%
  inner_join(regionBases, by = c("CB" = "CB", "read_region" = "region")) %>%
  mutate(rph_one = reads_per_region / bp_one,
         rph_sep = reads_per_region / bp_sep) %>%
  group_by(CB) %>%
  mutate(norm_one = sum(rph_one)/100,
         norm_sep = sum(rph_sep)/100) %>%
  ungroup() %>%
  mutate(tph_one = rph_one/norm_one,
         tph_sep = rph_sep/norm_sep) %>%
  mutate(regionorder = case_when(read_region == "utr5" ~ 0,
                                 read_region == "cds"  ~ 1,
                                 read_region == "utr3" ~ 2)) %>%
  mutate(read_region = fct_reorder(read_region, regionorder))


regionDisplay <- rbind(hekMeta %>% mutate(display = "scRibo-HEK293T") %>% select(CB, display),
                       fucciMeta %>% mutate(display = "scRibo-RPE1") %>% select(CB, display),
                       eeMeta %>% mutate(display = "scRibo-mEE") %>% select(CB, display),
                       bulkMeta %>% separate(CB, into = c("display"), sep = "_", extra = "drop", remove = FALSE) %>% select(CB, display)) %>%
  mutate(displayorder = case_when(display == "scRibo-HEK293T" ~ 0,
                                  display == "scRibo-RPE1"    ~ 1,
                                  display == "scRibo-mEE"     ~ 2,
                                  display == "darnell"        ~ 3,
                                  display == "ingolia"        ~ 4,
                                  display == "martinez"       ~ 5,
                                  display == "tanenbaum"      ~ 6)) %>% 
  mutate(display = fct_reorder(display, displayorder))

region_colours <- brewer.pal(7, "Set2")


# https://stackoverflow.com/questions/15660829/how-to-add-a-number-of-observations-per-group-and-use-group-mean-in-ggplot2-boxp
plt_region <- regionTPM %>%
  separate(CB, into = c("prefix"), extra = "drop", remove = FALSE) %>%
  inner_join(regionDisplay, by = "CB") %>%
  ggplot(aes(x=read_region, y=tph_one, fill = display))+
  geom_boxplot(aes(fill=display))+
  geom_point(position = position_jitterdodge(jitter.width = 0), size = 0.5, aes(colour = display)) +
  scale_colour_manual(values = region_colours)+
  scale_fill_manual(values = region_colours)+
  stat_summary(fun.data = function(x){ return( c(y=1, label = length(x)) ) }, geom = "text", fun = max, position = position_jitterdodge(jitter.width = 0, jitter.height = 0))
save_plot(paste0(figurePrefix, "_plt_region_tpm.pdf"), plt_region, base_width = 6, base_height = 5)




## length-frame heatmaps

frame_length_plot <- function(reads, frame, fill_range = NULL){
  frameL <- reads %>%
    group_by(length, .data[[frame]]) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    group_by(length) %>%
    mutate(length_total = sum(n)) %>%
    ungroup() 

  if(missing(fill_range)){
    fill_range = c(min(frameL$n/frameL$length_total), max(frameL$n/frameL$length_total))
  }
  plt_frame_hm <- ggplot(frameL, aes(x=.data[[frame]], y=length))+
    geom_tile(aes(fill = n/length_total))+
    scale_fill_gradientn(colours=(brewer.pal(9,"Greys")), limits = fill_range)+
    scale_y_continuous(expand=c(0,0))+
    scale_x_continuous(breaks=c(0,1,2),expand=c(0,0))+
    labs(fill = "fraction of reads/length")+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = "bottom")
  plt_len_hist <- frameL %>%
    select(length, length_total) %>%
    distinct() %>%
    mutate(total = sum(length_total)) %>%
    ggplot(aes(x=length, y=100*length_total/total))+
    geom_bar(stat = "identity", position = "dodge")+
    scale_x_continuous(expand=c(0,0))+
    coord_flip()+
    ylab("% of reads")+
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  plt_out <- ggarrange(plt_frame_hm, plt_len_hist,
                       nrow = 1, ncol = 2, widths = c(1,1), heights = c(1), draw = FALSE )

  return(plt_out)
}

multi_frame_length_plot <- function(reads, first_frame, second_frame, fill_range = NULL){
  frame1 <- reads %>%
    group_by(length, .data[[first_frame]]) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    group_by(length) %>%
    mutate(length_total = sum(n)) %>%
    ungroup() 
  frame2 <- reads %>%
    group_by(length, .data[[second_frame]]) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    group_by(length) %>%
    mutate(length_total = sum(n)) %>%
    ungroup() 

  if(missing(fill_range)){
    fill_range = c(min(c(frame1$n/frame1$length_total, frame2$n/frame2$length_total)),
                   max(c(frame1$n/frame1$length_total, frame2$n/frame2$length_total)))
  }

  plt_frame1_hm <- ggplot(frame1, aes(x=.data[[first_frame]], y=length))+
    geom_tile(aes(fill = n/length_total))+
    scale_fill_gradientn(colours=(brewer.pal(9,"Greys")), limits = fill_range)+
    scale_y_continuous(expand=c(0,0))+
    scale_x_continuous(breaks=c(0,1,2),expand=c(0,0))+
    labs(fill = "fraction of reads/length")+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = "bottom")
  plt_frame2_hm <- ggplot(frame2, aes(x=.data[[second_frame]], y=length))+
    geom_tile(aes(fill = n/length_total))+
    scale_fill_gradientn(colours=(brewer.pal(9,"Greys")), limits = fill_range)+
    scale_y_continuous(expand=c(0,0))+
    scale_x_continuous(breaks=c(0,1,2),expand=c(0,0))+
    labs(fill = "fraction of reads/length")+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position = "none",
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  plt_len_hist <- frame1 %>%
    select(length, length_total) %>%
    distinct() %>%
    mutate(total = sum(length_total)) %>%
    ggplot(aes(x=length, y=100*length_total/total))+
    geom_bar(stat = "identity", position = "dodge")+
    scale_x_continuous(expand=c(0,0))+
    coord_flip()+
    ylab("% of reads")+
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  plt_out <- ggarrange(plt_frame1_hm, plt_frame2_hm, plt_len_hist,
                       nrow = 1, ncol = 3, widths = c(1,1,1), heights = c(1), draw = FALSE )

  return(plt_out)
}

# find ranges
frame5 <- rbind(bulkReads, eeReads, fucciReads, hekReads) %>%
  inner_join(regionDisplay, by = "CB") %>%
  group_by(display, length, frame5) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(display, length) %>%
  mutate(length_total = sum(n)) %>%
  ungroup()
frameP <- rbind(bulkReads %>% filter(sample == "darnell"), eeReads, fucciReads, hekReads) %>%
  mutate(frameP = (psite-cds_start) %% 3) %>%
  inner_join(regionDisplay, by = "CB") %>%
  group_by(display, length, frameP) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(display, length) %>%
  mutate(length_total = sum(n)) %>%
  ungroup()

fill_range <- c(min(c(frame5$n/frame5$length_total, frameP$n/frameP$length_total)),
                max(c(frame5$n/frame5$length_total, frameP$n/frameP$length_total)))


plt_lhm_sc_HEK <- multi_frame_length_plot(reads = hekReads %>% mutate(frameP = (psite-cds_start) %% 3), first_frame = "frame5", second_frame = "frameP", fill_range = fill_range)
plt_lhm_sc_Fucci <- multi_frame_length_plot(reads = fucciReads %>% mutate(frameP = (psite-cds_start) %% 3), first_frame = "frame5", second_frame = "frameP", fill_range = fill_range)
plt_lhm_sc_ee <- multi_frame_length_plot(reads = eeReads %>% mutate(frameP = (psite-cds_start) %% 3), first_frame = "frame5", second_frame = "frameP", fill_range = fill_range)
plt_lhm_darnell <- multi_frame_length_plot(reads = bulkReads %>% filter(sample == "darnell") %>% mutate(frameP = (psite-cds_start) %% 3), first_frame = "frame5", second_frame = "frameP", fill_range = fill_range)
save_plot(paste0(figurePrefix, "_lengh_hm_scRibo_HEK.pdf"), plt_lhm_sc_HEK, base_width = 3, base_height = 3)
save_plot(paste0(figurePrefix, "_lengh_hm_scRibo_Fucci.pdf"), plt_lhm_sc_Fucci, base_width = 3, base_height = 3)
save_plot(paste0(figurePrefix, "_lengh_hm_scRibo_ee.pdf"), plt_lhm_sc_ee, base_width = 3, base_height = 3)
save_plot(paste0(figurePrefix, "_lengh_hm_darnell.pdf"), plt_lhm_darnell, base_width = 3, base_height = 3)

plt_lhm_ingolia <- frame_length_plot(reads = bulkReads %>% filter(sample == "ingolia"), frame = "frame5", fill_range = fill_range)
plt_lhm_martinez <- frame_length_plot(reads = bulkReads %>% filter(sample == "martinez"), frame = "frame5", fill_range = fill_range)
plt_lhm_tanenbaum <- frame_length_plot(reads = bulkReads %>% filter(sample == "tanenbaum"), frame = "frame5", fill_range = fill_range)
save_plot(paste0(figurePrefix, "_lengh_hm_ingolia.pdf"), plt_lhm_ingolia, base_width = 2.3, base_height = 3)
save_plot(paste0(figurePrefix, "_lengh_hm_martinez.pdf"), plt_lhm_martinez, base_width = 2.3, base_height = 3)
save_plot(paste0(figurePrefix, "_lengh_hm_tanenbaum.pdf"), plt_lhm_tanenbaum, base_width = 2.3, base_height = 3)


## contamination

biotypes <- as.data.frame(t(Reduce(reduceCountMatrix, list(hekBiotypes, fucciBiotypes, eeBiotypes, bulkBiotypes))))
biotypes$CB <- rownames(biotypes)
biotypes <- biotypes %>%
  pivot_longer(cols = starts_with(c("allbiot_", "biot_", "contam_")), names_to = "biotype", values_to = "counts") %>%
  inner_join(regionDisplay, by = "CB") 


plt_contam <- biotypes %>%
  filter(grepl("^contam_", biotype)) %>%
  group_by(CB) %>%
  mutate(cell_total = sum(counts)) %>%
  ungroup() %>%
  filter(biotype != "contam_Unaligned") %>%
  ggplot(aes(x=biotype, y=100*counts/cell_total, fill = display))+
  geom_boxplot()+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = 0), size = 0.5, aes(colour = display)) +
  scale_colour_manual(values = region_colours)+
  scale_fill_manual(values = region_colours)
save_plot(paste0(figurePrefix, "_plt_contam_fraction.pdf"), plt_contam, base_width = 6, base_height = 5)

allbiotypes_to_merge <- biotypes %>%
  filter(!grepl("rRNA", biotype)) %>%
  filter(grepl("^allbiot_", biotype)) %>%
  group_by(biotype) %>%
  summarize(biotype_total = sum(counts)) %>%
  ungroup() %>%
  arrange(-biotype_total) %>%
  mutate(biotype_fraction = biotype_total/sum(biotype_total)) %>%
  mutate(biotype_group = case_when(biotype_fraction >= 0.02 ~ biotype,
                                 biotype_fraction < 0.02 ~ "allbiot_other") ) %>%
  select(biotype, biotype_group)

## fraction of all reads that align to top biotypes
## total reads come from contam_ 
## merge infrequent biotypes

plt_biot <- biotypes %>%
  filter(grepl("^allbiot_", biotype)) %>%
  inner_join(allbiotypes_to_merge, by = "biotype") %>%
  group_by(CB, biotype_group) %>%
  summarize(group_counts = sum(counts)) %>%
  ungroup() %>% 
  inner_join(biotypes %>% 
             filter(grepl("contam_", biotype)) %>%
             group_by(CB) %>%
             summarize(total_reads = sum(counts)) %>%
             ungroup() %>%
             select(CB, total_reads), by = "CB") %>%
  inner_join(regionDisplay, by = "CB") %>%
  group_by(biotype_group) %>%
  mutate(biotype_average = mean(group_counts/total_reads)) %>%
  ungroup() %>%
  mutate(biotype_group = fct_reorder(biotype_group, -biotype_average)) %>%
  ggplot(aes(x=biotype_group, y=100*group_counts/total_reads, fill = display))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = 0), size = 0.5, aes(colour = display)) +
  scale_colour_manual(values = region_colours)+
  scale_fill_manual(values = region_colours)
save_plot(paste0(figurePrefix, "_plt_biotype_fraction.pdf"), plt_biot, base_width = 6, base_height = 5)


plt_box_mapping <- plot_grid(plt_contam+theme(legend.position = "bottom"),
                             plt_region+theme(legend.position = "bottom"),
                             plt_biot+theme(legend.position = "bottom"),
                             nrow = 1, align = "h", rel_widths = c(2,3,4))
save_plot(paste0(figurePrefix, "_plt_box_mapping.pdf"), plt_box_mapping, base_width = 10, base_height = 4)

## duplicate rate for scRibo-seq libraries
## calculate for protein-coding genes

plt_duplicate <- biotypes %>%
  filter(CB %in% c(hekMeta$CB, fucciMeta$CB, eeMeta$CB)) %>%
  filter(!grepl("^contam_", biotype)) %>%
  filter(grepl("protein_coding", biotype)) %>%
  pivot_wider(names_from = biotype, values_from = counts) %>%
  mutate(duplicate_fraction = (allbiot_protein_coding - biot_protein_coding)/allbiot_protein_coding) %>%
  ggplot(aes(x=1, y=100*duplicate_fraction, fill = display))+
  geom_boxplot()+
  scale_colour_manual(values = region_colours)+
  scale_fill_manual(values = region_colours)+
  ylab("% duplicate protein-coding")
  #stat_summary(fun.data = function(x){ return( c(y=1, label = length(x)) ) }, geom = "text", fun = max, position = position_jitterdodge(jitter.width = 0, jitter.height = 0))
save_plot(paste0(figurePrefix, "_plt_duplicates.pdf"), plt_duplicate, base_width = 4, base_height = 5)



### count
counts <- as.data.frame((Reduce(reduceCountMatrix, list(as(rowSums(hekCounts[, pull(hekMeta %>% filter(sort_population == "Rich") %>% select(CB))]), "sparseMatrix"),
                                                        as(rowSums(bulkCounts[,grep('^ingolia', colnames(bulkCounts))]), "sparseMatrix"),
                                                        as(rowSums(bulkCounts[,grep('^martinez', colnames(bulkCounts))]), "sparseMatrix"),
                                                        as(rowSums(bulkCounts[,grep('^darnell_Rich', colnames(bulkCounts))]), "sparseMatrix")
                                                         ))))
colnames(counts) <- c("scRibo", "ingolia", "martinez", "darnell")


normCounts <- 1E6*sweep(counts, 2, apply(counts,2,sum), FUN = "/")

#https://gastonsanchez.wordpress.com/2012/08/27/scatterplot-matrices-with-ggplot/
makePairs <- function(data){
  grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
  grid <- subset(grid, x != y)
  all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
                             xcol <- grid[i, "x"]
                             ycol <- grid[i, "y"]
                             data.frame(xvar = names(data)[ycol], yvar = names(data)[xcol], 
                                        x = data[, xcol], y = data[, ycol], data)
                                 }))
  all$xvar <- factor(all$xvar, levels = names(data))
  all$yvar <- factor(all$yvar, levels = names(data))
  densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
                                   data.frame(xvar = names(data)[i], yvar = names(data)[i], x = data[, i])
                                 }))
  list(all=all, densities=densities)
}

# expand iris data frame for pairs plot
gg1 = makePairs(normCounts)

megaCounts <- data.frame(gg1$all)

plt_HEK_counts_comparison <- ggplot(megaCounts, aes(x = log10(x+1), y =log10(y+1))) + 
  facet_grid(xvar ~ yvar) + 
  coord_equal()+
  geom_point(size=0.3, na.rm = TRUE, alpha=0.05)


save_plot(paste0(figurePrefix, "_plt_genecounts_HEK.pdf"), plt_HEK_counts_comparison, base_width = 6, base_height = 6)


cor(normCounts, method = "spearman")

## numbers for text 
# contamination

contam <- biotypes %>%
  filter(grepl("^contam_", biotype)) %>%
  group_by(CB) %>%
  mutate(cell_total = sum(counts)) %>%
  mutate(cell_pct = 100*counts/cell_total) %>%
  ungroup() %>%
  filter(biotype != "contam_Unaligned") %>%
  filter(grepl("^scRibo", display))
contam %>%
  group_by(biotype) %>%
  summarize(mn = mean(cell_pct),
            sem = sd(cell_pct)/sqrt(n())) %>%
  ungroup()

contam %>%
  group_by(biotype, display) %>%
  summarize(mn = mean(cell_pct),
            sem = sd(cell_pct)/sqrt(n())) %>%
  ungroup()
