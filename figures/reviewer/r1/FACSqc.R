#!/usr/bin/env Rscript
# Michael VanInsberghe 2021-03-30
# QC plots for FACS comparisons

library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
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
theme_set(theme_cowplot())

source('plotFunctions.R')

figurePrefix <- "facs/facsqc"

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

fucciReads <- readRDS(file.path(data_dir, "Fucci/compiled/reads.rds"))
fucciMeta <- readRDS(file.path(data_dir, "Fucci/compiled/meta.rds"))
fucciAllMeta <- readRDS(file.path(data_dir, "Fucci/compiled/allmeta.rds"))
fucciCounts <- readRDS(file.path(data_dir, "Fucci/compiled/counts.rds"))
fucciBiotypes <- readRDS(file.path(data_dir, "Fucci/compiled/biotypes.rds"))



plt_QC <- fucciAllMeta %>%
  filter(grepl("^SR2", plate)) %>%
  filter(grepl("^cycling", sort_population)) %>%
  mutate(plate = gsub("SR20", "", plate)) %>%
  ggplot(aes(x=log10(cell_total), y=cds_frac, colour = sort_population))+
  geom_point(size=0.5)+
  geom_hline(yintercept = 0.87, colour = 'red', linetype = 'dashed')+
  geom_vline(xintercept = log10(4500), colour = 'red', linetype = 'dashed')+
  facet_wrap(~plate+sort_population)
save_plot(paste0(figurePrefix, "_plt_QC.pdf"), plt_QC, base_width = 8, base_height = 6)


# fraction of cells passing QC
kept_cells <- pull(fucciMeta %>% select(CB))

plate_pass <- fucciAllMeta %>%
  filter(grepl("^SR2", plate)) %>%
  select(CB, plate, sort_population) %>%
  mutate(pass_qc = case_when(CB %in% kept_cells ~ "pass",
                             TRUE ~ "fail")) %>%
  filter(sort_population != "NTC") %>%
  filter(!(sort_population %in% c("KO14", "KO13", "KO15", "KO17",
                                  "9h_post", "5h_post", "2h_post", "20h_post", "10h_post", "6h_post"))) %>%
  group_by(plate, sort_population) %>%
  summarize(pass = sum(pass_qc == "pass"),
            fail = sum(pass_qc == "fail")) %>%
  ungroup() %>%
  mutate(pass_frac = pass/(pass+fail)) %>%
  mutate(plate = gsub("SR20", "", plate))

plt_plate_pass <- ggplot(plate_pass, aes(x=plate, y=pass_frac, colour = sort_population, group = sort_population))+
  geom_point()+
  geom_line()
  #theme(axis.text.x=element_text(angle = 90, hjust = 0))
save_plot(paste0(figurePrefix, "_plt_plate_QC_pass.pdf"), plt_plate_pass, base_width = 5, base_height = 4)



# length-frame

plt_length <- fucciReads %>%
  select(-starts_with('base')) %>%
  filter(CB %in% pull(fucciMeta %>% filter(grepl("^SR2", plate)) %>% filter(grepl("^cycling_", sort_population)) %>% select(CB))) %>%
  inner_join(fucciMeta %>% select(CB, plate), by = "CB") %>%
  mutate(plate = gsub("SR20", "", plate)) %>%
  ggplot(aes(x=length))+
  geom_histogram(binwidth = 1)+
  facet_wrap(~plate, ncol=1, scales = "free_y")
plt_frame <- fucciReads %>%
  select(-starts_with('base')) %>%
  filter(CB %in% pull(fucciMeta %>% filter(grepl("^SR2", plate)) %>% filter(grepl("^cycling_", sort_population)) %>% select(CB))) %>%
  inner_join(fucciMeta %>% select(CB, plate), by = "CB") %>%
  mutate(frameP = (psite-cds_start) %% 3) %>%
  mutate(plate = gsub("SR20", "", plate)) %>%
  group_by(CB, plate) %>%
  summarize(frameP0 = sum(frameP ==0)/n()) %>%
  ungroup() %>%
  ggplot(aes(x=plate, y=frameP0))+
  geom_boxplot()
save_plot(paste0(figurePrefix, "_plt_length_frame.pdf"), plot_grid(plt_length, plt_frame, nrow=1), base_width = 7, base_height = 6)


# contamination and biotypes
region_colours <- brewer.pal(7, "Set2")

biotypes <- as.data.frame(t(fucciBiotypes))
biotypes$CB <- rownames(biotypes)
biotypes <- biotypes %>%
  pivot_longer(cols = starts_with(c("allbiot_", "biot_", "contam_")), names_to = "biotype", values_to = "counts")

plt_contam <- biotypes %>%
  filter(grepl("^contam_", biotype)) %>%
  filter(CB %in% pull(fucciMeta %>% filter(grepl("^SR2", plate)) %>% filter(grepl("^cycling_", sort_population)) %>% select(CB))) %>%
  group_by(CB) %>%
  mutate(cell_total = sum(counts)) %>%
  ungroup() %>%
  filter(biotype != "contam_Unaligned") %>%
  inner_join(fucciMeta %>% select(CB, plate, sort_population), by = "CB") %>%
  mutate(plate = gsub("SR20", "", plate)) %>%
  ggplot(aes(x=plate, y=100*counts/cell_total, fill = sort_population))+
  geom_boxplot()+
  facet_wrap(~biotype, nrow = 1)+
  stat_summary(fun.data = function(x){ return( c(y=1, label = length(x)) ) }, geom = "text", fun = max, position = position_jitterdodge(jitter.width = 0, jitter.height = 0))
save_plot(paste0(figurePrefix, "_plt_contam_fraction.pdf"), plt_contam, base_width = 6, base_height = 5)

allbiotypes_to_merge <- biotypes %>%
  filter(CB %in% pull(fucciMeta %>% filter(grepl("^SR2", plate)) %>% filter(grepl("^cycling_", sort_population)) %>% select(CB))) %>%
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
  filter(CB %in% pull(fucciMeta %>% filter(grepl("^SR2", plate)) %>% filter(grepl("^cycling_", sort_population)) %>% select(CB))) %>%
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
  inner_join(fucciMeta %>% select(CB, plate, sort_population), by = "CB") %>%
  group_by(biotype_group) %>%
  mutate(biotype_average = mean(group_counts/total_reads)) %>%
  ungroup() %>%
  mutate(biotype_group = fct_reorder(biotype_group, -biotype_average)) %>%
  mutate(plate = gsub("SR20", "", plate)) %>%
  ggplot(aes(x=plate, y=100*group_counts/total_reads, fill = sort_population))+
  geom_boxplot()+
  facet_wrap(~biotype_group, nrow=1)+
  stat_summary(fun.data = function(x){ return( c(y=1, label = length(x)) ) }, geom = "text", fun = max, position = position_jitterdodge(jitter.width = 0, jitter.height = 0))
save_plot(paste0(figurePrefix, "_plt_biotype_fraction.pdf"), plt_biot, base_width = 7, base_height = 5)

## region

expanded_canonical <- canonical_pc %>% 
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

regionBases <- fucciReads %>%
  filter(CB %in% pull(fucciMeta %>% filter(grepl("^SR2", plate)) %>% filter(grepl("^cycling_", sort_population)) %>% select(CB))) %>%
  select(CB, transcript_id) %>%
  inner_join(expanded_canonical, by = "transcript_id") %>%
  group_by(CB) %>%
  summarize(utr5 = sum(lr_utr5),
            cds = sum(lr_cds),
            utr3 = sum(lr_utr3)) %>%
  ungroup() %>%
  pivot_longer(cols = c("utr5", "cds", "utr3"), names_to = "region", values_to = "bp_sep") %>%
  inner_join(fucciReads %>% 
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



regionTPM <- fucciReads %>% 
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


plt_region <- regionTPM %>%
  separate(CB, into = c("prefix"), extra = "drop", remove = FALSE) %>%
  inner_join(fucciMeta %>% select(CB, plate, sort_population), by = "CB") %>%
  mutate(plate = gsub("SR20", "", plate)) %>%
  ggplot(aes(x=plate, y=tph_one, fill = sort_population))+
  geom_boxplot()+
  facet_wrap(~read_region, nrow=1)+
  stat_summary(fun.data = function(x){ return( c(y=1, label = length(x)) ) }, geom = "text", fun = max, position = position_jitterdodge(jitter.width = 0, jitter.height = 0))
save_plot(paste0(figurePrefix, "_plt_region_tpm.pdf"), plt_region, base_width = 6, base_height = 5)




test <- grid.arrange(grobs = lapply( list(plt_contam + theme(legend.position = "none"),
                                          plt_biot + theme(legend.position = "none"),
                                          plt_region + theme(legend.position = "right")),
                                    set_panel_size,
                                    width = unit(1.1, "in"),
                                    height = unit(3.0, "in")),
                     nrow = 1)

save_plot(paste0(figurePrefix, "_plt_FACS_qc_all_box.pdf"), test, base_width = 20, base_height = 5)


## duplicates

plt_duplicate <- biotypes %>%
  filter(CB %in% pull(fucciMeta %>% filter(grepl("^SR2", plate)) %>% filter(grepl("^cycling_", sort_population)) %>% select(CB))) %>%
  filter(!grepl("^contam_", biotype)) %>%
  filter(grepl("protein_coding", biotype)) %>%
  pivot_wider(names_from = biotype, values_from = counts) %>%
  mutate(duplicate_fraction = (allbiot_protein_coding - biot_protein_coding)/allbiot_protein_coding) %>%
  inner_join(fucciMeta %>% select(CB, plate, sort_population), by = "CB") %>%
  mutate(plate = gsub("SR20", "", plate)) %>%
  ggplot(aes(x=plate, y=100*duplicate_fraction, fill = sort_population))+
  geom_boxplot()+
  scale_colour_manual(values = region_colours)+
  ylab("% duplicate protein-coding")
save_plot(paste0(figurePrefix, "_plt_duplicates.pdf"), plt_duplicate, base_width = 4, base_height = 5)


## count profiles

fucciLong <- as.data.frame(t(fucciCounts))
fucciLong$CB <- rownames(fucciLong)
fucciLong <- fucciLong %>%
  pivot_longer(cols = starts_with("ENSG"), names_to = "gene_id", values_to = "counts") 

group_counts <- fucciLong %>%
  filter(CB %in% pull(fucciMeta %>% filter(grepl("^SR2", plate)) %>% filter(grepl("^cycling_", sort_population)) %>% select(CB))) %>%
  filter(counts > 0) %>%
  inner_join(fucciMeta %>% select(CB, plate, sort_population), by = "CB") %>%
  group_by(plate, sort_population, gene_id) %>%
  summarize(group_counts = sum(counts)) %>%
  ungroup() %>%
  group_by(plate, sort_population) %>%
  mutate(total_counts = sum(group_counts)) %>%
  ungroup() %>%
  mutate(id = paste(plate, sort_population, sep = "__"),
         rpm = 1E6*group_counts/total_counts) %>%
  select(id,gene_id, rpm) %>%
  pivot_wider(names_from = "id", values_from = "rpm", values_fill = 0) %>%
  as.data.frame()
rownames(group_counts) <- group_counts$gene_id
group_counts <- group_counts[, colnames(group_counts) != "gene_id"]

CM <- cor(group_counts, method = "spearman")
diag(CM) <- NA

cors <- CM[upper.tri(CM)]
mean(cors)
sd(cors)/sqrt(length(cors))

pdf(paste0(figurePrefix, "_gene_count_cor.pdf"), height = 4, width = 7)
top_ann <- HeatmapAnnotation(plate = gsub("^(.*?)__.*", "\\1", colnames(CM), perl = TRUE),
                             sort_population = gsub("^.*__(.*?)$", "\\1", colnames(CM), perl = TRUE),
                             col = list(plate = c("SR201" = brewer.pal(6, "Set2")[1],
                                                  "SR202" = brewer.pal(6, "Set2")[2],
                                                  "SR203" = brewer.pal(6, "Set2")[3],
                                                  "SR204" = brewer.pal(6, "Set2")[4],
                                                  "SR205" = brewer.pal(6, "Set2")[5],
                                                  "SR206" = brewer.pal(6, "Set2")[6]),
                                        sort_population = c("cycling_post" = gg_color_hue(2)[1],
                                                            "cycling_pre" = gg_color_hue(2)[2])))
row_ann <- HeatmapAnnotation(plate = gsub("^(.*?)__.*", "\\1", colnames(CM), perl = TRUE),
                             sort_population = gsub("^.*__(.*?)$", "\\1", colnames(CM), perl = TRUE),
                             col = list(plate = c("SR201" = brewer.pal(6, "Set2")[1],
                                                  "SR202" = brewer.pal(6, "Set2")[2],
                                                  "SR203" = brewer.pal(6, "Set2")[3],
                                                  "SR204" = brewer.pal(6, "Set2")[4],
                                                  "SR205" = brewer.pal(6, "Set2")[5],
                                                  "SR206" = brewer.pal(6, "Set2")[6]),
                                        sort_population = c("cycling_post" = gg_color_hue(2)[1],
                                                            "cycling_pre" = gg_color_hue(2)[2])),
                             which = "row")
Heatmap(as.matrix(CM),
        right_annotation = row_ann,
        top_annotation = top_ann,
        #col = colorRamp2(seq(-1, 1, length = 11), rev(brewer.pal(11, "RdBu"))),
        na_col = "grey10",
        name = "spearman",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        show_column_names = FALSE,
        border = TRUE,
        use_raster = FALSE)
dev.off()

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

gg1 = makePairs(group_counts)

megaCounts <- data.frame(gg1$all)

plt_counts_comparison <- ggplot(megaCounts, aes(x = log10(x+1), y =log10(y+1))) + 
  facet_grid(xvar ~ yvar) + 
  coord_equal()+
  geom_point(size=1, na.rm = TRUE, alpha=0.05)+
  geom_abline(slope=1, intercept = 0, colour = 'black', linetype = 'dashed')
save_plot(paste0(figurePrefix, "_plt_genecounts.png"), plt_counts_comparison, base_width = 20, base_height = 20)

