#!/usr/bin/env Rscript
# Michael VanInsberghe 2021-03-22
# single-cell QC metrics

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
theme_set(theme_cowplot())

source('plotFunctions.R')

figurePrefix <- "scmetrics/scmetrics"

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


## scRPF data
#hekReads <- readRDS(file.path(data_dir, "HEK293T/compiled/reads.rds"))
hekMeta <- readRDS(file.path(data_dir, "HEK293T/compiled/meta.rds"))
hekAllMeta <- readRDS(file.path(data_dir, "HEK293T/compiled/allmeta.rds"))
hekCounts <- readRDS(file.path(data_dir, "HEK293T/compiled/counts.rds"))
hekBiotypes <- readRDS(file.path(data_dir, "HEK293T/compiled/biotypes.rds"))

#fucciReads <- readRDS(file.path(data_dir, "Fucci/compiled/reads.rds"))
fucciMeta <- readRDS(file.path(data_dir, "Fucci/compiled/meta.rds"))
fucciAllMeta <- readRDS(file.path(data_dir, "Fucci/compiled/allmeta.rds"))
fucciCounts <- readRDS(file.path(data_dir, "Fucci/compiled/counts.rds"))
fucciBiotypes <- readRDS(file.path(data_dir, "Fucci/compiled/biotypes.rds"))

#eeReads <- readRDS(file.path(data_dir, "EE/compiled/reads.rds"))
eeMeta <- readRDS(file.path(data_dir, "EE/compiled/meta.rds"))
eeAllMeta <- readRDS(file.path(data_dir, "EE/compiled/allmeta.rds"))
eeCounts <- readRDS(file.path(data_dir, "EE/compiled/counts.rds"))
eeBiotypes <- readRDS(file.path(data_dir, "EE/compiled/biotypes.rds"))


regionDisplay <- rbind(hekMeta %>% mutate(display = "scRibo-HEK293T") %>% select(CB, display),
                       fucciMeta %>% mutate(display = "scRibo-RPE1") %>% select(CB, display),
                       eeMeta %>% mutate(display = "scRibo-mEE") %>% select(CB, display)) %>%
  mutate(displayorder = case_when(display == "scRibo-HEK293T" ~ 0,
                                  display == "scRibo-RPE1"    ~ 1,
                                  display == "scRibo-mEE"     ~ 2,
                                  display == "darnell"        ~ 3,
                                  display == "ingolia"        ~ 4,
                                  display == "martinez"       ~ 5,
                                  display == "tanenbaum"      ~ 6)) %>% 
  mutate(display = fct_reorder(display, displayorder))

region_colours <- brewer.pal(7, "Set2")


## metrics per cell
# reads per cell
reads_per_cell <- bind_rows(hekMeta  %>% select(CB, cds),
                     fucciMeta %>% select(CB, cds),
                     eeMeta %>% select(CB, cds)) %>%
  inner_join(regionDisplay, by = "CB") 
plt_cds <- reads_per_cell %>%
  ggplot(aes(x=log10(cds), fill = display))+
  geom_histogram(bins = 50) +
  scale_fill_manual(values = region_colours)+
  facet_wrap(~display, scales = "free_y", ncol = 3)+
  geom_text(data = reads_per_cell %>% group_by(display) %>% summarize(mn = mean(cds), sem = sd(cds)/sqrt(n())) %>% ungroup() %>% mutate(txt = sprintf("%.0f ± %.0f", mn, sem)),
            mapping = aes(x= -Inf, y=-Inf, label = txt),
            hjust = -0.1, vjust = -1)+
  scale_x_continuous(breaks = scales::pretty_breaks(n=4))+
  theme(legend.position = "none")
#save_plot(paste0(figurePrefix, "_cds_reads_per_cell.pdf"), plt_cds, base_width = 6, base_height = 2) 

reads_per_cell %>%
  summarize(mn = mean(cds),
            sem = sd(cds)/sqrt(n()))

reads_per_cell %>%
  group_by(display) %>%
  summarize(mn = mean(cds),
            sem = sd(cds)/sqrt(n()))

# genes per cell
hekMeta$genes <- apply(hekCounts[, hekMeta$CB]>0, 2, sum)
fucciMeta$genes <- apply(fucciCounts[, fucciMeta$CB]>0, 2, sum)
eeMeta$genes <- apply(eeCounts[, eeMeta$CB]>0, 2, sum)

genes_per_cell<- bind_rows(hekMeta  %>% select(CB, genes),
                     fucciMeta %>% select(CB, genes),
                     eeMeta %>% select(CB, genes)) %>%
  inner_join(regionDisplay, by = "CB") %>%
  mutate(display = fct_reorder(display, displayorder))
plt_genes <- genes_per_cell %>%
  ggplot(aes(x=log10(genes), fill = display))+
  geom_histogram(bins = 40) +
  scale_fill_manual(values = region_colours)+
  facet_wrap(~display, scales = "free_y", ncol = 3)+
  geom_text(data = genes_per_cell %>% group_by(display) %>% summarize(mn = mean(genes), sem = sd(genes)/sqrt(n())) %>% ungroup() %>% mutate(txt = sprintf("%.0f ± %.0f", mn, sem)),
            mapping = aes(x= -Inf, y=-Inf, label = txt),
            hjust = -0.1, vjust = -1)+
  scale_x_continuous(breaks = scales::pretty_breaks(n=3))+
  theme(legend.position = "none")

#save_plot(paste0(figurePrefix, "_genes_per_cell.pdf"), plt_genes, base_width = 6, base_height = 2) 
save_plot(paste0(figurePrefix, "_QC_histograms.pdf"), plot_grid(plt_cds, plt_genes, nrow=2, align = "v"), base_width = 6.5, base_height = 4)

genes_per_cell %>%
  summarize(mn = mean(genes),
            sdn = sd(genes),
            sem = sd(genes)/sqrt(n()))

genes_per_cell %>%
  group_by(display) %>%
  summarize(mn = mean(genes),
            sdn = sd(genes),
            sem = sd(genes)/sqrt(n()))




## duplicate rate for scRibo-seq libraries
## calculate for protein-coding genes

biotypes <- as.data.frame(t(Reduce(reduceCountMatrix, list(hekBiotypes, fucciBiotypes, eeBiotypes))))
biotypes$CB <- rownames(biotypes)
biotypes <- biotypes %>%
  pivot_longer(cols = starts_with(c("allbiot_", "biot_", "contam_")), names_to = "biotype", values_to = "counts") %>%
  inner_join(regionDisplay, by = "CB") 

duplicate_per_cell <- biotypes %>%
  filter(CB %in% c(hekMeta$CB, fucciMeta$CB, eeMeta$CB)) %>%
  filter(!grepl("^contam_", biotype)) %>%
  filter(grepl("protein_coding", biotype)) %>%
  pivot_wider(names_from = biotype, values_from = counts) %>%
  mutate(duplicate_fraction = (allbiot_protein_coding - biot_protein_coding)/allbiot_protein_coding) 

plt_duplicate <- duplicate_per_cell %>%
  ggplot(aes(x=1, y=100*duplicate_fraction, fill = display))+
  geom_boxplot()+
  scale_colour_manual(values = region_colours)+
  scale_fill_manual(values = region_colours)+
  ylab("% duplicate protein-coding")
  #stat_summary(fun.data = function(x){ return( c(y=1, label = length(x)) ) }, geom = "text", fun = max, position = position_jitterdodge(jitter.width = 0, jitter.height = 0))
save_plot(paste0(figurePrefix, "_plt_duplicates.pdf"), plt_duplicate, base_width = 4, base_height = 4)

plt_duphist <- duplicate_per_cell %>%
  ggplot(aes(x=100*duplicate_fraction, fill = display))+
  scale_fill_manual(values = region_colours)+
  geom_histogram(bins = 40) +
  facet_wrap(~display, scales = "free_y", ncol = 3)+
  geom_text(data = duplicate_per_cell %>% group_by(display) %>% summarize(mn = mean(100*duplicate_fraction), sem = sd(100*duplicate_fraction)/sqrt(n())) %>% ungroup() %>% mutate(txt = sprintf("%.2f ± %.2f", mn, sem)),
            mapping = aes(x= -Inf, y=-Inf, label = txt),
            hjust = -0.1, vjust = -1)+
  theme(legend.position = "none")

save_plot(paste0(figurePrefix, "_QC_histograms.pdf"), plot_grid(plt_cds, plt_genes, plt_duphist, ncol=1, align = "v"), base_width = 6.5, base_height = 6)


duplicate_per_cell %>%
  summarize(mn = mean(100*duplicate_fraction),
            sem = sd(100*duplicate_fraction)/sqrt(n()))
## QC plots
selectedPlates <- data.frame(plate = c("HEK293T", "HEK293T-Starv1", "HEK293T-Starv2",
                                        "Mit01", "SR201","MGI02",
                                        "Pro03", "Mid02", "Dis01"),
                             plateorder = seq(1,9),
                             cds_cutoff = c(0.75, 0.75, 0.75,
                                            0.87, 0.87, 0.87,
                                            0.75, 0.75, 0.75),
                             reads_cutoff = c(1000, 1000, 1000,
                                              4500, 4500, 4500,
                                              1800, 2000, 2000),
                             stringsAsFactors = FALSE) %>%
  mutate(plate = fct_reorder(plate, plateorder))

plt_qc_example <- bind_rows(hekAllMeta %>% select(CB, cds_frac, cell_total, plate, sort_population),
                            fucciAllMeta %>% select(CB, cds_frac, cell_total, plate, sort_population),
                            eeAllMeta %>% 
                              mutate(ssc_gated = case_when(sort_population == "NTC" ~ "NTC",
                                                                       sort_population != "NTC" ~ ssc_gated)) %>%
                              filter(ssc_gated != "doublet") %>%
                              select(CB, cds_frac, cell_total, plate, sort_population)) %>%
  mutate(type = case_when(sort_population == "NTC" ~ "empty",
                          TRUE ~ "cell")) %>%
  inner_join(selectedPlates, by = "plate") %>%
  mutate(plate = fct_reorder(plate, plateorder)) %>%
  arrange(type) %>%
  ggplot(aes(x=log10(cell_total), y=cds_frac, colour = type))+
  geom_point(size = 0.5)+
  scale_colour_manual(values = c("grey30", "#fb9a99"))+
  geom_vline(data = selectedPlates, aes(xintercept = log10(reads_cutoff)), colour = "grey60", linetype = "dashed")+
  geom_hline(data = selectedPlates, aes(yintercept = cds_cutoff), colour = "grey60", linetype = "dashed")+
  facet_wrap(~plate, ncol=3)
save_plot(paste0(figurePrefix, "_qc_example_selected.pdf"), plt_qc_example, base_width = 6.25, base_height = 5.5)



## fraction of cells passing QC

kept_cells <- pull(bind_rows(hekMeta %>% select(CB),
                             fucciMeta %>% select(CB),
                             eeMeta %>% select(CB)) %>% select(CB))

plate_pass <- bind_rows(hekAllMeta %>% select(CB, plate, sort_population),
                            fucciAllMeta %>% select(CB, plate, sort_population),
                            eeAllMeta %>% filter(ssc_gated == "single") %>% select(CB, plate, sort_population)) %>%
  mutate(pass_qc = case_when(CB %in% kept_cells ~ "pass",
                             TRUE ~ "fail")) %>%
  filter(sort_population != "NTC") %>%
  filter(!(sort_population %in% c("KO14", "KO13", "KO15", "KO17",
                                  "9h_post", "5h_post", "2h_post", "20h_post", "10h_post", "6h_post"))) %>%
  group_by(plate) %>%
  summarize(pass = sum(pass_qc == "pass"),
            fail = sum(pass_qc == "fail")) %>%
  ungroup() %>%
  mutate(pass_frac = pass/(pass+fail)) %>%
  inner_join(bind_rows(hekMeta %>% select(CB, plate),
                       fucciMeta %>% select(CB, plate), 
                       eeMeta %>% select(CB, plate)) %>%
             inner_join(regionDisplay, by = "CB") %>%
             select(plate, display) %>%
             distinct(), by = "plate") 

library(ggbeeswarm)
plt_plate_pass <-  ggplot(plate_pass, aes(x=display, y=pass_frac, colour = display))+
  geom_beeswarm()+
  scale_colour_manual(values = region_colours)+
  scale_fill_manual(values = region_colours)+
  stat_summary(fun.data = function(x){ return( c(y=0.75, label = length(x)) ) }, geom = "text", fun = max, position = position_jitterdodge(jitter.width = 0, jitter.height = 0))
save_plot(paste0(figurePrefix, "_plt_QC_pass.pdf"), plt_plate_pass, base_width = 6, base_height = 4)

## plate pass numbers for text
# On average, 
plate_pass %>%
  summarize(mn = mean(pass_frac))

plate_pass %>%
  group_by(display) %>%
  summarize(mn = mean(100*pass_frac),
            sd = sd(100*pass_frac),
            sem = sd(100*pass_frac)/sqrt(n())) %>%
  ungroup()
