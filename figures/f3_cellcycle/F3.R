#!/usr/bin/env Rscript
# Michael VanInsberghe 2020-10-27
# Fucci RPF cell cycle GAA pausing, et al. figures

library(dplyr)
library(tidyr)
library(tibble)
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
theme_set(theme_cowplot())

source('/hpc/hub_oudenaarden/mvanins/local/avopipelines/RPF/bin/scRibo/plotFunctions.R')

figurePrefix <- "raw/F3"

# barcodes
barcodes <- read.table('/hpc/hub_oudenaarden/mvanins/local/lib/barcodes/miRv6_barcodemap.tsv',
                       sep="\t",
                       header = FALSE,
                       stringsAsFactors = FALSE,
                       row.names = NULL,
                       col.names = c("CB","number","well"))
rownames(barcodes) <- barcodes$CB

# annot
annot <- read_feather('../../data/hsa_annotations.feather')

canonical_pc <- annot %>%
  filter(set=="canonical", transcript_type == "protein_coding") %>%
  #filter(chr != 'chrM') %>%
  filter(transcript_id != "ENST00000437139.7") %>%
  filter(!grepl("^PCDH", gene_name)) #%>%
  #filter(!grepl("^RPL", gene_name))

codons <- read_feather('../../data/hsa_codons.feather')
codons <- codons %>%
  filter( transcript_id %in% canonical_pc$transcript_id)

## RPF sample definition
RPFsamples <- data.frame(names=c("G011",
                                 "G012",
                                 "Mit01",
                                 "Mit02",
                                 "Mit03",
                                 "Mit04",
                                 "Mit06",
                                 "Mit07",
                                 "Mit08",
                                 "KO01",
                                 "KO02",
                                 "KO03",
                                 "KO04",
                                 "SR201",
                                 "SR202",
                                 "SR203",
                                 "SR204",
                                 "SR205",
                                 "SR206"),
                         counts = c("../../data/Fucci/counts/RPFv4-Fucci-G011",
                                    "../../data/Fucci/counts/RPFv4-Fucci-G012",
                                    "../../data/Fucci/counts/RPFv4-Fucci-Mit01",
                                    "../../data/Fucci/counts/RPFv4-Fucci-Mit02",
                                    "../../data/Fucci/counts/RPFv4-Fucci-Mit03",
                                    "../../data/Fucci/counts/RPFv4-Fucci-Mit04",
                                    "../../data/Fucci/counts/RPFv4-Fucci-Mit06",
                                    "../../data/Fucci/counts/RPFv4-Fucci-Mit07",
                                    "../../data/Fucci/counts/RPFv4-Fucci-Mit08",
                                    "../../data/Fucci/counts/RPFv4-CDKN1AKO01",
                                    "../../data/Fucci/counts/RPFv4-CDKN1AKO02",
                                    "../../data/Fucci/counts/RPFv4-CDKN1AKO03",
                                    "../../data/Fucci/counts/RPFv4-CDKN1AKO04",
                                    "../../data/Fucci/counts/RPFv4-SR201",
                                    "../../data/Fucci/counts/RPFv4-SR202",
                                    "../../data/Fucci/counts/RPFv4-SR203",
                                    "../../data/Fucci/counts/RPFv4-SR204",
                                    "../../data/Fucci/counts/RPFv4-SR205",
                                    "../../data/Fucci/counts/RPFv4-SR206"),
                         facs = c("../../data/Fucci/facs/G011.csv",
                                  "../../data/Fucci/facs/G012.csv",
                                  "../../data/Fucci/facs/Mit01.csv",
                                  "../../data/Fucci/facs/Mit02.csv",
                                  "../../data/Fucci/facs/Mit03.csv",
                                  "../../data/Fucci/facs/Mit04.csv",
                                  "../../data/Fucci/facs/Mit06.csv",
                                  "../../data/Fucci/facs/Mit07.csv",
                                  "../../data/Fucci/facs/Mit08.csv",
                                  "../../data/Fucci/facs/KO01.csv",
                                  "../../data/Fucci/facs/KO02.csv",
                                  "../../data/Fucci/facs/KO03.csv",
                                  "../../data/Fucci/facs/KO04.csv",
                                  "../../data/Fucci/facs/SR201.csv",
                                  "../../data/Fucci/facs/SR202.csv",
                                  "../../data/Fucci/facs/SR203.csv",
                                  "../../data/Fucci/facs/SR204.csv",
                                  "../../data/Fucci/facs/SR205.csv",
                                  "../../data/Fucci/facs/SR206.csv"),
                         reads = c("../../data/Fucci/reads/RPFv4-Fucci-G011_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-Fucci-G012_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-Fucci-Mit01_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-Fucci-Mit02_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-Fucci-Mit03_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-Fucci-Mit04_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-Fucci-Mit06_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-Fucci-Mit07_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-Fucci-Mit08_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-CDKN1AKO01_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-CDKN1AKO02_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-CDKN1AKO03_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-CDKN1AKO04_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-SR201_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-SR202_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-SR203_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-SR204_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-SR205_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz",
                                   "../../data/Fucci/reads/RPFv4-SR206_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                         stringsAsFactors = FALSE)

## cells meeting cutoff filter, but are much lower than cells in similar phase of cell cycle 
s_low <- c("G011_GTACTCCGTC", "Mit01_ATCCTAACAG", "Mit01_TTGTCTCTAC", "Mit02_GACCTACGAG", "Mit02_TAAGCCTCTA", "Mit02_TGGTGTTGAC", "Mit03_ACAGAGACCG",
           "Mit03_AGTTCTTCCG", "Mit03_CCGACAAGGA", "Mit03_TCAGAGGCAC", "Mit04_AACCGTAACA", "Mit04_CTAATCGCGG", "Mit04_GTCTTGGAAT", "Mit04_TCTCACTCTA",
           "Mit06_CCGTGAGAAC", "Mit06_CGATACTGAA", "Mit06_GAACCTCCTA", "Mit06_GCAGCGTTAG", "Mit06_TAAGCCTCTA", "Mit07_CAACTGTCAT", "Mit07_TTGCTCACGG",
           "Mit08_AATTCGGTGG", "Mit08_CAACTGTCAT", "KO01_CCGACAAGGA", "SR201_CGAGACTTCA", "SR201_GAAGTCGTGG", "SR201_TAGGTGGTGG", "SR201_TCAAGGTCGC",
           "SR202_ACGCGCAACA", "SR202_AGTCACAACA", "SR202_CGTCGGTAAG", "SR202_TTGGTGTGTC", "SR204_TATCCGATCT")

RPFcounts <- readMergeCounts(RPFsamples)

RPFfacs <- readFACS(RPFsamples, barcodes = barcodes) %>%
  mutate(green = X..488..530.40,
         red = X..561..585.29)

RPFreads <- canonical_pc %>%
  inner_join(fread("../../data/Fucci/reads/mergedFilteredReads.csv.gz"), by = "transcript_id") %>%
  select(-id)

## Assemble metadata information
RPFmeta <- data.frame(row.names = colnames(RPFcounts),
                      CB = colnames(RPFcounts),
                      well = barcodes[gsub(".*_(.*?)$","\\1",colnames(RPFcounts),perl=TRUE),"well"],
                      plate = gsub("^(.*?)[_].*","\\1",colnames(RPFcounts),perl=TRUE),
                      cds = RPFcounts["cds",],
                      utr3 = RPFcounts["utr3",],
                      utr5 = RPFcounts["utr5",],
                      cds_frac = RPFcounts["cds_frac",],
                      cell_total = RPFcounts["cell_total",],
                      sort_population = rep("Interphase",ncol(RPFcounts)),
                      type = rep("RPF", ncol(RPFcounts)),
                      stringsAsFactors = FALSE)

RPFmeta[ (grepl("G0", RPFmeta$plate) & RPFmeta$well %in% columnFeatures(1:8, LETTERS[1:16])), "sort_population"] = "G0"
RPFmeta[ (grepl("Mit", RPFmeta$plate) & RPFmeta$well %in% columnFeatures(1:8, LETTERS[1:16])), "sort_population"] = "Mit"
RPFmeta[ (grepl("Mit|G0", RPFmeta$plate) & RPFmeta$well %in% columnFeatures(21:24, LETTERS[15:16])), "sort_population"] = "NTC"

RPFmeta[ (grepl("^KO", RPFmeta$plate) & RPFmeta$well %in% columnFeatures(1:4, LETTERS[1:16]) ), "sort_population"] = "KO17"
RPFmeta[ (grepl("^KO", RPFmeta$plate) & RPFmeta$well %in% columnFeatures(5:8, LETTERS[1:16]) ), "sort_population"] = "KO15"
RPFmeta[ (grepl("^KO", RPFmeta$plate) & RPFmeta$well %in% columnFeatures(9:12, LETTERS[1:16]) ), "sort_population"] = "KO14"
RPFmeta[ (grepl("^KO", RPFmeta$plate) & RPFmeta$well %in% columnFeatures(13:16, LETTERS[1:16]) ), "sort_population"] = "KO13"
RPFmeta[ (grepl("^KO", RPFmeta$plate) & RPFmeta$well %in% columnFeatures(17:24, LETTERS[1:16]) ), "sort_population"] = "WT"
RPFmeta[ (grepl("^KO", RPFmeta$plate) & RPFmeta$well %in% columnFeatures(22:24, LETTERS[15:16]) ), "sort_population"] = "NTC"

RPFmeta[ (grepl("^SR2", RPFmeta$plate) & RPFmeta$well %in% c(columnFeatures(1:4, LETTERS[1:16]), columnFeatures(5, LETTERS[1:8])) ), "sort_population"] = "2h_post"
RPFmeta[ (grepl("^SR2", RPFmeta$plate) & RPFmeta$well %in% c(columnFeatures(6:9, LETTERS[1:16]), columnFeatures(5, LETTERS[9:16])) ), "sort_population"] = "6h_post"
RPFmeta[ (grepl("^SR2", RPFmeta$plate) & RPFmeta$well %in% c(columnFeatures(10:13, LETTERS[1:16]), columnFeatures(14, LETTERS[1:8])) ), "sort_population"] = "10h_post"
RPFmeta[ (grepl("^SR2", RPFmeta$plate) & RPFmeta$well %in% c(columnFeatures(15:18, LETTERS[1:16]), columnFeatures(14, LETTERS[9:16])) ), "sort_population"] = "20h_post"
RPFmeta[ (grepl("^SR2", RPFmeta$plate) & RPFmeta$well %in% columnFeatures(19:21, LETTERS[1:16]) ), "sort_population"] = "cycling_post"
RPFmeta[ (grepl("^SR2", RPFmeta$plate) & RPFmeta$well %in% columnFeatures(22:24, LETTERS[1:14]) ), "sort_population"] = "cycling_pre"
RPFmeta[ (grepl("^SR2", RPFmeta$plate) & RPFmeta$well %in% columnFeatures(22:24, LETTERS[15:16]) ), "sort_population"] = "NTC"

RPFmeta[ (grepl("^KO", RPFmeta$plate)), "batch"] = "RPF_KO"
RPFmeta[ (grepl("Mit|G0", RPFmeta$plate)), "batch"] = "RPF_First"
RPFmeta[ (grepl("^SR2", RPFmeta$plate)), "batch"] = "RPF_SR2"


RPFmeta <- RPFmeta %>%
  full_join(RPFfacs %>% select(-Well, -sample), by = "CB")
rownames(RPFmeta) <- RPFmeta$CB
RPFmeta <- RPFmeta[colnames(RPFcounts),]

RPFcounts <- RPFcounts[!(rownames(RPFcounts) %in% c("cds","cds_frac","cell_total","utr3","utr5")),]


ENSG_to_name <- data.frame(count_names = rownames(RPFcounts), stringsAsFactors = FALSE) %>%
  inner_join(annot, by = c("count_names" = "gene_id")) %>% 
  mutate(gene_out = paste(count_names, gene_name, sep = "-")) %>%
  select(count_names, gene_out) %>%
  distinct() %>%
  column_to_rownames(var = "count_names")

rownames(RPFcounts) <- ENSG_to_name[rownames(RPFcounts),"gene_out"]
RPFcounts <- RPFcounts[!grepl("-PCDH", rownames(RPFcounts)),]


## QC check
# plt_cd_totalQC <- ggplot(RPFmeta %>% filter(cell_total > 100), aes(x=log10(cell_total), y=cds_frac))+
#   geom_point()+
#   geom_vline(xintercept = log10(4500), colour = 'red', linetype = 'dashed')+
#   geom_hline(yintercept = .85, colour = 'red', linetype = 'dashed')+
#   facet_wrap(~plate, ncol = 7)
# save_plot(paste0(figurePrefix, "_plt_QC.pdf"), plt_cd_totalQC, base_width = 17, base_height = 9)

RPFcells <- pull(RPFmeta %>% 
                 filter(cds_frac >= 0.85, cell_total > 4500) %>%
                 filter(!(sort_population %in% c("NTC"))) %>%
                 filter(!(batch %in% c("Sync"))) %>%
                 filter(!(sort_population %in% c("KO14", "KO13", "KO15", "KO17",
                                       "9h_post", "5h_post", "2h_post", "20h_post", "10h_post", "6h_post"))) %>%
                 filter(!(CB %in% s_low)) %>%
                 select(CB))

RPFreads <- RPFreads %>%
  filter(CB %in% RPFcells)
RPFcounts <- RPFcounts[,RPFcells]
RPFfacs <- RPFfacs %>%
  filter(CB %in% RPFcells)
RPFmeta <- RPFmeta[RPFcells,]

# RPF slopes
lr_utr5 = 50
lr_cds = 100
lr_utr3 = 50
sc_pscaled <- scaleDF(RPFreads, 
                      group = "single", site = "psite", 
                      utr5 = lr_utr5, cds = lr_cds, utr3 = lr_utr3) %>%
mutate(CB = as.character(CB)) %>%
separate(CB, c("sample",NA), remove = FALSE, sep="_") %>%
mutate(fraction = n/tot)

RPFcoverage <- as.data.frame(sc_pscaled %>%
                             select(CB, pos, fraction) %>%
                             spread(pos, fraction))

rownames(RPFcoverage) <- RPFcoverage$CB
RPFcoverage <- RPFcoverage[, -which(colnames(RPFcoverage) %in% "CB")]

CDS <- RPFcoverage[rownames(RPFmeta),(lr_utr5+1):(lr_utr5+lr_cds)]
fitslope <- function(x) {
  ll <- lm(seq(0,1,length.out = 100) ~ x)
  return(ll$coefficients[2])
}

RPFmeta$rpf_slope <- apply(CDS,1,fitslope)

## Tours
tours <- readRDS('/hpc/hub_oudenaarden/mvanins/RPF/analysis/scRiboSeq/figures/f3/trajectory/RPF_IS_trajectory_bigger_all_raw_tours_5e+05iter_30pts_1noise.RDS') %>%
  group_by(CB) %>%
  summarize(median_path = median(path_distance),
            mean_path = mean(path_distance),
            n_tours = n()) %>% 
  ungroup() %>%
  mutate(path = case_when(median_path == 0 ~ mean_path,
                          TRUE ~ median_path)) %>%
  select(CB, n_tours, path) %>%
  as.data.frame()
rownames(tours) <- tours$CB

RPFmeta <- cbind(RPFmeta, 
                 tours[rownames(RPFmeta),c("n_tours","path")])
RPFmeta <- RPFmeta[!is.na(RPFmeta$CB),]


RPFmeta$log_green <- log10(RPFmeta$green)
RPFmeta$log_red <- log10(RPFmeta$red)

P <- CreateSeuratObject(counts = RPFcounts,
                        project = "RPF_Fucci",
                        assay = "RNA",
                        min.cells = 10,
                        min.features = 10,
                        names.field = 1,
                        names.delim = "-",
                        meta.data = RPFmeta)
P[["percent_mt"]] <- PercentageFeatureSet(P, pattern = "-MT-")

## Batch correction/Integration

P.list <- SplitObject(object = P, split.by = "batch")
for(i in 1:length(x=P.list)){
  P.list[[i]] <- NormalizeData(object = P.list[[i]], verbose = FALSE)
  #   P.list[[i]] <- FindVariableFeatures(object = P.list[[i]], selection.method = "vst", nfeatures = 1000, verbose = FALSE)
  P.list[[i]] <- FindVariableFeatures(object = P.list[[i]], selection.method = "mvp", verbose = FALSE)
}

integration_dims <- 1:30
P.anchors <- FindIntegrationAnchors(object.list = P.list[unique(RPFmeta$batch)], dims = integration_dims, k.filter = 120)

PI <- IntegrateData(anchorset = P.anchors, dims=integration_dims)

DefaultAssay(object = PI) <- "integrated"

PI <- ScaleData(object = PI, verbose=FALSE)
PI <- RunPCA(object=PI, verbose=FALSE)
PI <- JackStraw(PI, num.replicate = 100)
PI <- ScoreJackStraw(PI, dims = 1:20)

ndims <- 13
PI <- RunUMAP(PI, dims=1:ndims, verbose=FALSE, n.neighbors = 50, n.epochs = 1000, return.model = TRUE)
PI <- FindNeighbors(PI, dims=1:ndims, k.param=15)
PI <- FindClusters(PI,  n.start=50, n.iter=50)

RPFmeta$seurat_clusters <- PI$seurat_clusters[rownames(RPFmeta)]
RPFmeta$percent_mt <- PI$percent_mt[rownames(RPFmeta)]
RPFmeta$UMAP_1 <- as.data.frame(Embeddings(PI, reduction="umap"))[rownames(RPFmeta),"UMAP_1"]
RPFmeta$UMAP_2 <- as.data.frame(Embeddings(PI, reduction="umap"))[rownames(RPFmeta),"UMAP_2"]

# collapse sort populations
RPFmeta$population <- pull(RPFmeta %>%
                           mutate(population = case_when(sort_population == "G0" ~ "G0",
                                                         sort_population == "Interphase" ~ "Interphase",
                                                         sort_population == "Mit" ~ "Mitotic",
                                                         sort_population == "WT" ~ "Interphase",
                                                         sort_population == "cycling_pre" ~ "Interphase",
                                                         sort_population == "cycling_post" ~ "Interphase",
                                                         TRUE ~ "unknown") ) %>%
                           select(population))
RPFmeta$population_order <- pull(RPFmeta %>%
                                 mutate(population_order = case_when(population == "G0" ~ 0,
                                                                     population == "Mitotic" ~ 1,
                                                                     population == "Interphase" ~ 2,
                                                                     TRUE ~ -1) ) %>%
                                 select(population_order) ) 

## export integrated info for GP pseudotime
# remove G0 cells
export_cells <- pull(RPFmeta %>%
                     #filter(seurat_clusters != 6) %>%
                     select(CB))

## z-score fluorescence only:
UF <- t(RPFmeta[export_cells, c("log_red", "log_green")])
UF_mean_correction <- apply(UF, 1, mean)
UF_sd_correction <- apply(sweep(UF, 1, UF_mean_correction, FUN = "-"), 1, sd)

UF <- sweep(UF, 1, UF_mean_correction, FUN = "-")
UF <- sweep(UF, 1, UF_sd_correction, FUN = "/")

UZ <- rbind(as.data.frame(PI@assays$integrated@scale.data[,export_cells]),
            UF)


PUZ <- rbind(t(Embeddings(PI, reduction = 'pca'))[,export_cells],
             t(RPFmeta[export_cells, c("log_red", "log_green")]))

write.csv(UZ, file = paste0(figurePrefix, "_RPF_seurat_integrated.csv"))
write.csv(RPFmeta[colnames(UZ),], file = paste0(figurePrefix, "_RPF_seurat_integrated_meta.csv"))
write.csv((PUZ), file=paste0(figurePrefix, "_RPF_seurat_integrated_pca.csv"))

rm(UF, UZ, PUZ)


## read pseudotimes 
# pst_dir <- "/hpc/hub_oudenaarden/mvanins/RPF/analysis/CS2integration/pausing_fig/prediction_nolowS/"
pst_dir <- "pseudotimes_g0/"

prediction_mean <- as.data.frame(fread(paste0(pst_dir, "prediction_withfacs_posterior_mean.csv.gz")))
prediction_var <- as.data.frame(fread(paste0(pst_dir, "prediction_withfacs_posterior_var.csv.gz")))
pseudotimes <- fread(paste0(pst_dir,"pseudotimes_withfacs.csv.gz"))
cellnames <- fread(paste0(pst_dir, "cellNames.csv.gz"))

cellnames <- as.data.frame(cellnames[!is.na(cellnames$V1),"V2"])$V2

RPFmeta2 <- RPFmeta[cellnames,]
RPFmeta2$pseudotime <- as.data.frame(pseudotimes)$V1

RPFmeta$pseudotime <- NA
RPFmeta[cellnames, "pseudotime"] <- as.data.frame(pseudotimes)$V1


min_pst <- sort(RPFmeta2$pseudotime)[4]
max_pst <- sort(RPFmeta2$pseudotime)[nrow(RPFmeta2)-6]

RPFmeta$pseudotime <- pull(RPFmeta %>%
                           mutate(pst_bounded = case_when(pseudotime < min_pst ~ min_pst,
                                                          pseudotime > max_pst ~ max_pst,
                                                          TRUE ~ pseudotime)) %>%
                           mutate(pst_scaled = (pst_bounded - min_pst)/(max_pst-min_pst)) %>%
                           select(pst_scaled))

RPFmeta2$pseudotime <- pull(RPFmeta2 %>%
                           mutate(pst_bounded = case_when(pseudotime < min_pst ~ min_pst,
                                                          pseudotime > max_pst ~ max_pst,
                                                          TRUE ~ pseudotime)) %>%
                           mutate(pst_scaled = (pst_bounded - min_pst)/(max_pst-min_pst)) %>%
                           select(pst_scaled))
min_pst <- 0
max_pst <- 1

## pseudotime checks
# plt_pst_tsp <- ggplot(RPFmeta2, aes(x=path, y=pseudotime, colour = seurat_clusters))+
#   geom_point()+
#   coord_equal()
# save_plot(paste0(figurePrefix, "_pseudotime_tsp_comparison.pdf"), plt_pst_tsp, base_width = 6, base_height=6)
# plt_umap_pst <- ggplot(RPFmeta2, aes(x=UMAP_1, y=UMAP_2, colour = pseudotime))+
#   geom_point()+
#   coord_equal()+
#   scale_colour_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), na.value = "grey40")
# plt_umap_path <- ggplot(RPFmeta2, aes(x=UMAP_1, y=UMAP_2, colour = path))+
#   geom_point()+
#   coord_equal()+
#   scale_colour_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), na.value = "grey40")
# save_plot(paste0(figurePrefix, "_pseudotime_umap.pdf"), plot_grid(plt_umap_path, plt_umap_pst, nrow=1), base_width = 16, base_height=6)
# 
# plt_pst_g <- ggplot(RPFmeta2, aes(x=pseudotime, y=log_green))+
#   geom_point()
# plt_pst_r <- ggplot(RPFmeta2, aes(x=pseudotime, y=log_red))+
#   geom_point()
# save_plot(paste0(figurePrefix, "pst_rgfacs.pdf"), plot_grid(plt_pst_r, plt_pst_g, ncol=1), base_width=5, base_height=5)


## Fluorescence prediction reprojection

F_mean <- prediction_mean[, names(UF_mean_correction)]
F_mean <- sweep(F_mean, 2, UF_sd_correction, FUN = "*")
F_mean <- sweep(F_mean, 2, UF_mean_correction, FUN = "+")
rownames(F_mean) <- 1:nrow(F_mean)
F_mean <- as.data.frame(F_mean)
F_mean$pseudotime <- seq(from = min_pst, 
                      to = max_pst,
                      length.out = nrow(F_mean))

F_var <- prediction_var[, c("log_green", "log_red")]
rownames(F_var) <- 1:nrow(F_var)

# plt_FACS_pst <- ggplot(RPFmeta, aes(x=log_green, y=log_red, colour = seurat_clusters))+
#   geom_point()+
#   geom_point(data=F_mean[seq(from=1, to=nrow(F_mean), length.out = 100), ], colour = "black", alpha = 0.5)+
#   geom_path(data=F_mean[seq(from=1, to=nrow(F_mean), length.out = 100), ], colour = "black", arrow=arrow(length=unit(0.5, "cm"), type="closed",), alpha = 0.5)+
#   coord_equal()
# save_plot(paste0(figurePrefix, "_plt_facs_pstpred.pdf"), plt_FACS_pst)
# 
# 
# plt_FACS_pst <- ggplot(F_mean[seq(from = 1, to = nrow(F_mean), length.out = 100), ], aes(x=log_green, y=log_red, colour = pseudotime))+
#   geom_point(data = RPFmeta %>% select(log_green, log_red), colour = "grey60")+
#   geom_point()+
#   scale_colour_distiller(palette = "Spectral")+
#   #scale_colour_gradientn(colours = rev(brewer.pal(9,"YlGnBu")), na.value = "grey40")+
#   coord_equal()
# save_plot(paste0(figurePrefix, "_plt_facs_pstpred.pdf"), plt_FACS_pst)
# 
# plt_pfacs_green <- ggplot(F_mean, aes(x=pseudotime, y=log_green))+
#   geom_point()
# plt_pfacs_red <- ggplot(F_mean, aes(x=pseudotime, y=log_red))+
#   geom_point()
# save_plot(paste0(figurePrefix, "_plt_predicted_facs.pdf"), plot_grid(plt_pfacs_green, plt_pfacs_red, ncol=1))


## project pseudotime path in umap

# PMR <- P_mean[!rownames(P_mean) %in% c("log_green", "log_red"),]
PMR <- t(prediction_mean[,colnames(prediction_mean) != "V1"])
colnames(PMR) <- 1:ncol(PMR)
PMR <- PMR[!rownames(PMR) %in% c("log_green", "log_red"),]

RMR <- t(PMR) %*% PI@reductions$pca@feature.loadings %*% diag(PI@reductions$pca@stdev*sqrt(max(1, nrow(PI@reductions$pca@cell.embeddings)-1)))
u_t <- as.data.frame(umap_transform(RMR[,1:ndims],  model = PI@reductions$umap@misc$model))
colnames(u_t) <- c("UMAP_1", "UMAP_2")

u_t$pseudotime <- seq(from = min_pst,
                      to = max_pst,
                      length.out = nrow(u_t))

# u_t <- u_t[seq(from=1, to=nrow(u_t), length.out=500),]

# plt_umap_pst <- ggplot(RPFmeta, aes(x=UMAP_1, y=UMAP_2, colour = seurat_clusters))+
#   geom_point()+
#   geom_point(data = as.data.frame(u_t), colour = "black", alpha = 0.6)+ 
#   geom_path(data = as.data.frame(u_t), colour = "black", arrow=arrow(length=unit(0.5, "cm"), type="closed"), alpha = 0.6)+
#   coord_equal()
# save_plot(paste0(figurePrefix, "_plt_umap_pstpred.pdf"), plt_umap_pst)
# 
# plt_umap_pst <- ggplot(u_t, aes(x=UMAP_1, y=UMAP_2, colour = pseudotime))+
#   geom_point(data = RPFmeta %>% select(UMAP_1, UMAP_2), colour = "grey40")+
#   geom_point()+
#   scale_colour_distiller(palette = "Spectral")+
#   #scale_colour_gradientn(colours = rev(brewer.pal(9,"YlGnBu")), na.value = "grey40")+
#   coord_equal()
# save_plot(paste0(figurePrefix, "_plt_umap_pstpred.pdf"), plt_umap_pst)


#population_colours <- c("#F2AD00", "#5BBCD6", "#00A08A")
#population_colours <- c("#ECCBAE", "#5BBCD6", "#00A08A")
#population_colours <- c("#4daf4a", "#377eb8","#984ea3") # G0, Interphase, Mit
population_colours <- c("#f781bf", "#549ed1","#984ea3") # G0, Interphase, Mit

## reorder clusters
RPFmeta$clusters <- factor(pull(RPFmeta %>%
                                mutate(clusters = case_when(seurat_clusters == 5 ~ 0,
                                                            seurat_clusters == 2 ~ 1,
                                                            seurat_clusters == 0 ~ 2,
                                                            seurat_clusters == 3 ~ 3,
                                                            seurat_clusters == 4 ~ 4,
                                                            seurat_clusters == 7 ~ 5,
                                                            seurat_clusters == 1 ~ 6,
                                                            seurat_clusters == 6 ~ 7,
                                                            TRUE ~ -1)) %>%
                                select(clusters)))

### FACS plots
set.seed(112)
plt_facs_clusters <- ggplot(RPFmeta[sample(nrow(RPFmeta), replace=FALSE),], aes(x=log_green, y=log_red, colour = clusters))+
  geom_point()+
  xlab("log10 mAG-GMNN")+
  ylab("log10 mKO2-CDT1")+
  coord_equal() #+
#   theme(axis.line = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
plt_facs_population <- ggplot(RPFmeta %>% arrange(-population_order), aes(x=log_green, y=log_red, colour = population))+
  geom_point()+
  coord_equal()+
  xlab("log10 mAG-GMNN")+
  ylab("log10 mKO2-CDT1")+
  scale_colour_manual(values = population_colours) #+
#   theme(axis.line = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
plt_facs_pseudotime_proj <- ggplot(F_mean, aes(x=log_green, y=log_red, colour = pseudotime))+
  geom_point(data = RPFmeta %>% select(log_green, log_red), colour = "grey50")+
  geom_point()+
  xlab("log10 mAG-GMNN")+
  ylab("log10 mKO2-CDT1")+
  scale_colour_gradientn(colours = rev(brewer.pal(11,"Spectral")), na.value = "grey40")+
  coord_equal() #+
#   theme(axis.line = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
#save_plot(paste0(figurePrefix, "_FACS_panels.pdf"), plot_grid(plt_facs_clusters, plt_facs_population, plt_facs_pseudotime_proj, nrow=1), base_width = 10, base_height = 3)

plt_umap_clusters <- ggplot(RPFmeta[sample(nrow(RPFmeta), replace=FALSE),], aes(x=UMAP_1, y=UMAP_2, colour = clusters))+
  geom_point()+
  coord_equal() +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_population <- ggplot(RPFmeta %>% arrange(-population_order), aes(x=UMAP_1, y=UMAP_2, colour = population))+
  geom_point()+
  scale_colour_manual(values = population_colours)+
  coord_equal()+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_pseudotime <- ggplot(RPFmeta2 %>% arrange(pseudotime), aes(x=UMAP_1, y=UMAP_2, colour = pseudotime))+
  geom_point()+
  coord_equal()+
  scale_colour_gradientn(colours = rev(brewer.pal(11,"Spectral")), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_pseudotime_path <- ggplot(u_t, aes(x=UMAP_1, y=UMAP_2, colour = pseudotime))+
  geom_point(data = RPFmeta %>% select(UMAP_1, UMAP_2), colour = "grey50", size=1)+
  geom_point()+
  scale_colour_gradientn(colours = rev(brewer.pal(11,"Spectral")), na.value = "grey40")+
  coord_equal()+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
#save_plot(paste0(figurePrefix, "_UMAP_panels.pdf"), plot_grid(plt_umap_clusters, plt_umap_population, plt_umap_pseudotime_path, nrow=1), base_width = 10, base_height = 3)
plt_umap_red <- ggplot(RPFmeta[sample(nrow(RPFmeta), replace=FALSE),], aes(x=UMAP_1, y=UMAP_2, colour = log_red))+
  geom_point()+
  coord_equal()+
  scale_colour_gradientn(colours = brewer.pal(8,"Reds"), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_green <- ggplot(RPFmeta[sample(nrow(RPFmeta), replace=FALSE),], aes(x=UMAP_1, y=UMAP_2, colour = log_green))+
  geom_point()+
  coord_equal()+
  scale_colour_gradientn(colours = brewer.pal(8,"Greens"), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
#save_plot(paste0(figurePrefix, "_UMAP_redgreen.pdf"), plot_grid(plt_umap_green, plt_umap_red, ncol=1), base_width = 10, base_height = 3)

save_plot(paste0(figurePrefix, "_reduction_panels.pdf"), plot_grid(plt_facs_population, plt_facs_clusters, plt_facs_pseudotime_proj,
                                                                   plt_umap_population, plt_umap_clusters, plt_umap_pseudotime_path,
                                                                   plt_umap_green, plt_umap_red, ncol = 3, align="hv"), base_width = 10, base_height = 10)

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
  inner_join(RPFmeta %>% select(CB, red, green, seurat_clusters, sort_population, plate, pseudotime), by = "CB") %>%
  inner_join(total_baseline, by = "codon_display") %>%
  inner_join(site_baseline, by = c("codon_display" = "codon_display", "site" = "site")) %>%
  mutate(total_fc = cell_total_fraction/(baseline_total_fraction/nrow(site_order)),
         site_fc = cell_site_fraction/baseline_site_fraction) %>%
  inner_join(site_order, by = "site") %>%
  mutate(site = fct_reorder(site, site_order))



## Pausing statistics
pausing_site <- site_change %>%
  filter(nchar(codon_display) == 7) %>%
  mutate(codon_display_site = paste(codon_display, site, sep = "_")) %>%
  # select(CB, codon_display_site, site_fc) %>%
  # pivot_wider(names_from = codon_display_site, values_from = site_fc) %>%
  select(CB, codon_display_site, cell_site_fraction) %>%
  pivot_wider(names_from = codon_display_site, values_from = cell_site_fraction) %>%
  as.data.frame()
rownames(pausing_site) <- pausing_site$CB
pausing_site <- pausing_site[, colnames(pausing_site) != "CB"]
pausing_site <- t(pausing_site)

PIp <- PI
PIp[['pausing_site']] <- CreateAssayObject(data = pausing_site)
pausing_markers <- FindAllMarkers(PIp, assay = "pausing_site", logfc.threshold = 0, only.pos = FALSE, pseudocount.use = 0) %>%
  separate(gene, c("codon", "three", "site"), remove = FALSE, sep = "-") %>%
  mutate(codon_display = paste(codon, three, sep = "_")) %>%
  mutate(cluster_redef = case_when(cluster == 5 ~ 0,
                              cluster == 2 ~ 1,
                              cluster == 0 ~ 2,
                              cluster == 3 ~ 3,
                              cluster == 4 ~ 4,
                              cluster == 7 ~ 5,
                              cluster == 1 ~ 6,
                              cluster == 6 ~ 7,
                              TRUE ~ -1)) %>%
  mutate(avg_logFC = log2(exp(avg_logFC)))

label_codons <- c("GAA_Glu", "GAG_Glu", "AUA_Ile", "CGA_Arg" )
plt_significant_codons <- pausing_markers %>%
  mutate(display = case_when(codon_display %in% label_codons ~ paste(paste(codon_display, site, sep = "-"), cluster_redef, sep = "."),
                             TRUE ~ "")) %>%
  mutate(display = case_when(p_val_adj < 1E-5 & abs(avg_logFC) > .25 ~ display,
                             TRUE ~ "")) %>%
  ggplot(aes(x=avg_logFC, y=-log10(p_val_adj), label = display))+
  geom_point()+
  geom_text_repel(segment.color = "black", box.padding = .6, size = 3, seed = 112, max.iter = 10000)+
  geom_vline(xintercept = .25, colour = "grey50", linetype = "dashed")+
  geom_vline(xintercept = -.25, colour = "grey50", linetype = "dashed")+
  geom_hline(yintercept = -log10(1E-5), colour = "grey50", linetype = "dashed")+
  xlim(-.75,.75)+
  xlab("average log2 FC")+
  ylab("-log10 adjusted p-value")
save_plot(paste0(figurePrefix, "_plt_sigsite.pdf"), plt_significant_codons, base_width = 6, base_height = 6)


sigcod_plt <- function(label_codons) {
  plt <- pausing_markers %>%
    mutate(display = case_when(codon_display %in% label_codons ~ paste(paste(codon_display, site, sep = "-"), cluster_redef, sep = "."),
                               TRUE ~ "")) %>%
    mutate(display = case_when(p_val_adj < 1E-5 & abs(avg_logFC) > .25 ~ display,
                               TRUE ~ "")) %>%
    ggplot(aes(x=avg_logFC, y=-log10(p_val_adj), label = display))+
    geom_point()+
    geom_text_repel(segment.color = "black", box.padding = .6, size = 3, seed = 112, max.iter = 10000)+
    geom_vline(xintercept = .25, colour = "grey50", linetype = "dashed")+
    geom_vline(xintercept = -.25, colour = "grey50", linetype = "dashed")+
    geom_hline(yintercept = -log10(1E-5), colour = "grey50", linetype = "dashed")+
    xlim(-.75,.75)+
    xlab("average log2 FC")+
    ylab("-log10 adjusted p-value")
  return(plt)
}

plt_sigcod_gaa <- sigcod_plt("GAA_Glu")
plt_sigcod_gag <- sigcod_plt("GAG_Glu")
plt_sigcod_aua <- sigcod_plt("AUA_Ile")
plt_sigcod_cga <- sigcod_plt("CGA_Arg")
plt_sigcod_ugc <- sigcod_plt("UGC_Cys")
plt_sigcod_cgc <- sigcod_plt("CGC_Arg")
plt_sigcod_cgu <- sigcod_plt("CGU_Arg")

save_plot(paste0(figurePrefix, "_plt_sigsites_by_codon.pdf"), plot_grid(plt_sigcod_ugc, plt_sigcod_cgc, plt_sigcod_cgu, plt_sigcod_gaa, plt_sigcod_gag, plt_sigcod_aua, plt_sigcod_cga, nrow = 1), base_width = 16, base_height = 3)
save_plot(paste0(figurePrefix, "_plt_sigsite.pdf"), plt_significant_codons, base_width = 6, base_height = 6)



sig_codons <- pull(pausing_markers %>%
                   filter(site %in% c("e", "p", "a")) %>%
                  filter(abs(avg_logFC) > .25) %>%
                  filter(p_val_adj < 1E-5) %>%
                  select(codon_display) %>%
                  distinct())

# sig_asite_codons <- pull(pausing_markers %>%
#   filter(abs(avg_logFC) > .25) %>%
#   filter(p_val_adj < 1E-5) %>%
#   filter(codon_display %in% sig_codons) %>%
#   select(codon_display, site) %>%
#   distinct() %>%
#   group_by(codon_display) %>%
#   summarize(n=n()) %>%
#   ungroup() %>%
#   filter(n < 3) %>%
#   select(codon_display))
## the above approach is noisy and relies on the threshold cutoffs,
## also, AAA_Lys is a bit of an outlier in terms of A-site occupancies

sig_asite_codons <- c("AUA_Ile", "CGA_Arg", "GAA_Glu", "GAG_Glu")

clusters_for_asite_codons <- pausing_markers %>%
  filter(abs(avg_logFC) > .25) %>%
  filter(p_val_adj < 1E-5) %>%
  filter(codon_display %in% sig_asite_codons) %>%
  group_by(cluster_redef) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  arrange(n)


## Pausing heatmaps
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


## selected codons pausing heatmap

selected_codons <- c(
                     "CAG_Gln",
                     "UGC_Cys", 
                     "CGC_Arg", "CGU_Arg", 
                     "GAA_Glu", "GAG_Glu", "AUA_Ile", "CGA_Arg")

selectedHM <- site_change %>%
  codons_heatmap(codons = selected_codons, feature = site_fc)
selectedHM <- t(selectedHM)
selectedHM <- selectedHM[, pull(RPFmeta %>% filter(!is.na(pseudotime)) %>% arrange(pseudotime) %>% select(CB))]

pdf(paste0(figurePrefix, "_selected_codon_stalling_sitefc_l2fc_pseudotime.pdf"), height = 8, width = 10)
#site_labels <- HeatmapAnnotation(site = anno_text(rep(c("-2","-1","E","P","A","+1","+2"), times = length(selected_codons)), rot=0, gp = gpar(fontsize = 5)), which="row")
# site_labels <- HeatmapAnnotation(site = anno_mark(at = c(grep("_a$", rownames(selectedHM)),
#                                                          grep("_p$", rownames(selectedHM)),
#                                                          grep("_e$", rownames(selectedHM))),
#                                                   labels = c(rep("A", length(grep("_a$", rownames(selectedHM)))),
#                                                              rep("P", length(grep("_p$", rownames(selectedHM)))),
#                                                              rep("E", length(grep("_e$", rownames(selectedHM))))) ) , which = "row")
site_labels <- HeatmapAnnotation(site = anno_mark(at = grep("^GAA_Glu", rownames(selectedHM)),
                                                  labels = c("-2", "-1", "E", "P", "A", "+1", "+2")), which = "row")
cell_annotation <- HeatmapAnnotation(population = as.character(RPFmeta[colnames(selectedHM), "population"]),
                                     clusters = as.character(RPFmeta[colnames(selectedHM), "clusters"]),
                                     pseudotime = RPFmeta[colnames(selectedHM), "pseudotime"],
                                     mKO2 = RPFmeta[colnames(selectedHM), "log_red"],
                                     mAG = RPFmeta[colnames(selectedHM), "log_green"],
                                     col = list(clusters = c("0" = gg_color_hue(8)[1],
                                                             "1" = gg_color_hue(8)[2],
                                                             "2" = gg_color_hue(8)[3],
                                                             "3" = gg_color_hue(8)[4],
                                                             "4" = gg_color_hue(8)[5],
                                                             "5" = gg_color_hue(8)[6],
                                                             "6" = gg_color_hue(8)[7],
                                                             "7" = gg_color_hue(8)[8]),
                                                population = c("G0" = population_colours[1],
                                                               "Interphase" = population_colours[2],
                                                               "Mitotic" = population_colours[3]),
                                                mKO2 = colorRamp2(seq(0, max(RPFmeta$log_red), length.out = 8), brewer.pal(8, 'Reds')),
                                                mAG = colorRamp2(seq(0, max(RPFmeta$log_red), length.out = 8), brewer.pal(8, 'Greens')),
                                                pseudotime = colorRamp2(seq(min(RPFmeta$pseudotime), max(RPFmeta$pseudotime), length.out = 11), rev(brewer.pal(11, 'Spectral')))),
                                     annotation_name_side = "left")
Heatmap(log2(as.matrix(selectedHM)),
        right_annotation = site_labels,
        top_annotation = cell_annotation,
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
        # row_split = factor(RPFmeta[rownames(all_hm), "seurat_clusters"], levels = levels(RPFmeta$seurat_clusters), ordered = TRUE),
        row_split = factor(rep(selected_codons, each = length(levels(site_change$site))), levels = selected_codons, ordered = TRUE),
        row_title_rot = 0,
        border = TRUE,
        use_raster = FALSE,
        raster_quality = 10)
dev.off()


## selected codons pausing boxplots
# n_pst <- 20
# plt_selected_box_bin <- site_change %>%
#   filter(codon_display %in% selected_codons) %>%
#   #filter(site %in% c("e", "p", "a")) %>%
#   filter(site %in% c("a")) %>%
#   mutate(pst_bin = round(n_pst*((pseudotime - min_pst)/(max_pst-min_pst)))) %>%
#   mutate(pst_bin = case_when(pst_bin < 0 ~ 0,
#                              pst_bin > n_pst ~ n_pst,
#                              TRUE ~ pst_bin)) %>%
#   ggplot(aes(x=as.factor(pst_bin), y=log2(site_fc), colour = site))+
#   geom_boxplot()+
#   facet_wrap(~codon_display, ncol=2)
# save_plot(paste0(figurePrefix, "_selected_codons_box_bin.pdf"), plt_selected_box_bin, base_width = 8, base_height = 10)
# 
# plt_selected_box_cluster <- site_change %>%
#   filter(codon_display %in% selected_codons) %>%
#   inner_join(RPFmeta %>% select(CB, clusters), by = "CB") %>%
#   #filter(site %in% c("e", "p", "a")) %>%
#   filter(site %in% c("a")) %>%
#   mutate(pst_bin = round(n_pst*((pseudotime - min_pst)/(max_pst-min_pst)))) %>%
#   mutate(pst_bin = case_when(pst_bin < 0 ~ 0,
#                              pst_bin > n_pst ~ n_pst,
#                              TRUE ~ pst_bin)) %>%
#   ggplot(aes(x=as.factor(clusters), y=log2(site_fc), colour = site))+
#   geom_boxplot()+
#   facet_wrap(~codon_display, ncol=2)
# save_plot(paste0(figurePrefix, "_selected_codons_box_clusters.pdf"), plt_selected_box_cluster, base_width = 8, base_height = 10)
# 
# plt_selected_pst_sca <- site_change %>%
#   filter(codon_display %in% selected_codons) %>%
#   filter(site %in% c("e", "p", "a")) %>%
#   # filter(site %in% c("a")) %>%
#   ggplot(aes(x=pseudotime, y=log2(site_fc), colour = site))+
#   geom_point()+
#   facet_wrap(~codon_display, ncol=2)
# save_plot(paste0(figurePrefix, "_selected_codons_sca_pst.pdf"), plt_selected_pst_sca, base_width = 8, base_height = 10)


## site pausing codons
# highlight_codons <- c("GAA_Glu", "AUA_Ile", "CGA_Arg", "GAG_Glu")
highlight_codons <- c("GAA_Glu", "AUA_Ile")

plt_umap_highlight <- site_change %>%
  filter(codon_display %in% highlight_codons) %>%
  filter(site %in% c("a")) %>%
  inner_join(RPFmeta %>% select(CB, UMAP_1, UMAP_2), by = "CB") %>%
  mutate(arrange_value = case_when(codon_display == "CGA_Arg" ~ -site_fc,
                                   TRUE ~ site_fc)) %>%
  arrange(arrange_value) %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, colour = log2(site_fc)))+
  geom_point(size=1)+
  coord_equal()+
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu")), na.value = "grey40", limits = c(-2,2), oob=squish)+
  facet_wrap(~codon_display)+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank())
save_plot(paste0(figurePrefix, "_highlight_codons.pdf"), plt_umap_highlight, base_width = 10, base_height = 6)

plt_aua_mit <- site_change %>%
  filter(codon_display == "AUA_Ile") %>%
  filter(site %in% c("a")) %>%
  inner_join(RPFmeta %>% select(clusters, CB, log_green), by = "CB") %>%
  filter(clusters == 6) %>%
  ggplot(aes(x=log_green, y=log2(site_fc)))+
  geom_point()
save_plot(paste0(figurePrefix, "_highlight_aua.pdf"), plt_aua_mit, base_width = 5, base_height = 5)




## all codons pausing heatmap

all_codons <- unique(site_change$codon_display)
all_codons <- all_codons[order(gsub(".*_(.*?)$","\\1",all_codons,perl=TRUE))]
all_codons <- all_codons[nchar(all_codons)==7]


# average fc across all sites to cluster codons
allHMa <- site_change %>%
  filter(codon_display %in% all_codons) %>%
  filter(site %in% c("em2", "em1", "e", "p", "a", "ap1", "ap2")) %>%
  group_by(CB, codon_display) %>%
  summarize(mean_site_fc = mean(site_fc),
            mean_total_fc = mean(total_fc)) %>%
  ungroup() %>%
  select(CB, codon_display, mean_site_fc) %>%
  spread(codon_display, mean_site_fc, fill = 0) %>%
  as.data.frame()
rownames(allHMa) <- allHMa$CB
allHMa <- t(allHMa[,all_codons])

allHMa_hc <- hclust(as.dist(1-cosine_similarity(as.matrix(log2(allHMa)))), method = "ward.D2")

all_codons <- all_codons[allHMa_hc$order]

allHM <- site_change %>%
  codons_heatmap(codons = all_codons, feature = site_fc)
allHM <- t(allHM)

allHM <- allHM[, pull(RPFmeta %>% filter(!is.na(pseudotime)) %>% arrange(pseudotime) %>% select(CB))]

all_codons_sig <- all_codons
all_codons_sig[all_codons_sig %in% sig_codons] <- paste0("*", all_codons_sig[all_codons_sig %in% sig_codons]) 

pdf(paste0(figurePrefix, "_all_codon_stalling_sitefc_l2fc_pseudotime.pdf"), height = 10, width = 10)
site_labels <- HeatmapAnnotation(site = anno_mark(at = grep("^GAA_Glu", rownames(allHM)),
                                                  labels = c("-2", "-1", "E", "P", "A", "+1", "+2")), which = "row")
#site_labels <- HeatmapAnnotation(site = anno_text(rep(c("-2","-1","E","P","A","+1","+2"), times = length(all_codons)), rot=0, gp = gpar(fontsize = 2)), which="row")
cell_annotation <- HeatmapAnnotation(population = as.character(RPFmeta[colnames(allHM), "population"]),
                                     clusters = as.character(RPFmeta[colnames(allHM), "clusters"]),
                                     pseudotime = RPFmeta[colnames(allHM), "pseudotime"],
                                     mKO2 = RPFmeta[colnames(allHM), "log_red"],
                                     mAG = RPFmeta[colnames(allHM), "log_green"],
                                     col = list(clusters = c("0" = gg_color_hue(8)[1],
                                                             "1" = gg_color_hue(8)[2],
                                                             "2" = gg_color_hue(8)[3],
                                                             "3" = gg_color_hue(8)[4],
                                                             "4" = gg_color_hue(8)[5],
                                                             "5" = gg_color_hue(8)[6],
                                                             "6" = gg_color_hue(8)[7],
                                                             "7" = gg_color_hue(8)[8]),
                                                population = c("G0" = population_colours[1],
                                                               "Interphase" = population_colours[2],
                                                               "Mitotic" = population_colours[3]),
                                                mKO2 = colorRamp2(seq(0, max(RPFmeta$log_red), length.out = 8), brewer.pal(8, 'Reds')),
                                                mAG = colorRamp2(seq(0, max(RPFmeta$log_red), length.out = 8), brewer.pal(8, 'Greens')),
                                                pseudotime = colorRamp2(seq(min(RPFmeta$pseudotime), max(RPFmeta$pseudotime), length.out = 11), rev(brewer.pal(11, 'Spectral')))),
                                     annotation_name_side = "left")
Heatmap(log2(as.matrix(allHM)),
        right_annotation = site_labels,
        top_annotation = cell_annotation,
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
        # row_split = factor(RPFmeta[rownames(all_hm), "seurat_clusters"], levels = levels(RPFmeta$seurat_clusters), ordered = TRUE),
        row_split = factor(rep(all_codons_sig, each = length(levels(site_change$site))), levels = all_codons_sig, ordered = TRUE),
        row_title_rot = 0,
        border = FALSE,
        use_raster = TRUE,
        raster_device = "CairoPNG",
        raster_quality = 15)
dev.off()


## Expression heatmap

DefaultAssay(PI) <- "RNA"
markers <- FindAllMarkers(PI, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

topmarkers <- markers %>% 
  filter(p_val_adj < 0.01) %>% 
  filter(!grepl("-AC[0-9][0-9]+", gene)) %>% 
  filter(!grepl("-AL[0-9][0-9]+", gene)) %>% 
  group_by(cluster) %>% 
  top_n(n = 5000, wt = avg_logFC)



outmarkers <- topmarkers %>%
  ungroup() %>%
  mutate(out_cluster = case_when(cluster == 5 ~ 0,
                              cluster == 2 ~ 1,
                              cluster == 0 ~ 2,
                              cluster == 3 ~ 3,
                              cluster == 4 ~ 4,
                              cluster == 7 ~ 5,
                              cluster == 1 ~ 6,
                              cluster == 6 ~ 7,
                              TRUE ~ -1)) %>%
  select(-cluster)
write.csv(outmarkers, file = paste0(figurePrefix, "_outmarkers.csv"))


markerHM <- as.data.frame(PI@assays$RNA[ pull(topmarkers %>% ungroup %>% select(gene) %>% distinct()), pull(RPFmeta2 %>% arrange(pseudotime) %>% select(CB)) ])

markerHM <- sweep(markerHM, 1, apply(markerHM,1,mean), FUN = "-")
markerHM <- sweep(markerHM, 1, apply(markerHM,1,sd), FUN = "/")

log2_fc_RPFpC = log2(RPFmeta2[colnames(markerHM), "cds"]/mean(RPFmeta2[colnames(markerHM), "cds"]))
pdf(paste0(figurePrefix, "_cc_marker_hm_pstwithfacs.pdf"), height = 9, width = 9)
cell_annotation <- HeatmapAnnotation(population = as.character(RPFmeta[colnames(markerHM), "population"]),
                                     clusters = as.character(RPFmeta[colnames(markerHM), "clusters"]),
                                     pseudotime = RPFmeta[colnames(markerHM), "pseudotime"],
                                     mKO2 = RPFmeta[colnames(markerHM), "log_red"],
                                     mAG = RPFmeta[colnames(markerHM), "log_green"],
                                     #log2_fc_RPFpC = anno_lines(log2(RPFmeta2[colnames(markerHM), "cds"]/mean(RPFmeta2[colnames(markerHM), "cds"])), height=unit(2, "cm"), smooth=TRUE, pch=20, gp = gpar(col="#b30000"), pt_gp = gpar(col = "grey80")),
                                     #log2_fc_RPFpC = log2_fc_RPFpC, 
                                     #cds = (log2(RPFmeta2[colnames(markerHM), "cds"]/mean(RPFmeta2[colnames(markerHM), "cds"]))),
                                     col = list(clusters = c("0" = gg_color_hue(8)[1],
                                                             "1" = gg_color_hue(8)[2],
                                                             "2" = gg_color_hue(8)[3],
                                                             "3" = gg_color_hue(8)[4],
                                                             "4" = gg_color_hue(8)[5],
                                                             "5" = gg_color_hue(8)[6],
                                                             "6" = gg_color_hue(8)[7],
                                                             "7" = gg_color_hue(8)[8]),
                                                population = c("G0" = population_colours[1],
                                                               "Interphase" = population_colours[2],
                                                               "Mitotic" = population_colours[3]),
                                                mKO2 = colorRamp2(seq(0, max(RPFmeta$log_red), length.out = 8), brewer.pal(8, 'Reds')),
                                                mAG = colorRamp2(seq(0, max(RPFmeta$log_red), length.out = 8), brewer.pal(8, 'Greens')),
                                                pseudotime = colorRamp2(seq(min(RPFmeta$pseudotime), max(RPFmeta$pseudotime), length.out = 11), rev(brewer.pal(11, 'Spectral'))),
                                                log2_fc_RPFpC = colorRamp2(seq(min(log2_fc_RPFpC), max(log2_fc_RPFpC), length.out = 11), rev(brewer.pal(11, 'PuOr')))),
                                     annotation_name_side = "right")
gene_annotation <- rowAnnotation(foo = anno_mark(at = c(
                                                        grep("-FBXO5$", rownames(markerHM)),
                                                        grep("-POLE$", rownames(markerHM)),
                                                        grep("-CCND1$", rownames(markerHM)),
                                                        grep("-MCM2$", rownames(markerHM)),
                                                        grep("-MCM5$", rownames(markerHM)),
                                                        grep("-H3C2$", rownames(markerHM)),
                                                        grep("-H1-4$", rownames(markerHM)),
                                                        grep("-MKI67$", rownames(markerHM)),
                                                        grep("-CENPE$", rownames(markerHM)),
                                                        grep("-UBE2C$", rownames(markerHM)),
                                                        grep("-CDK1$", rownames(markerHM)),
                                                        grep("-CCNB1$", rownames(markerHM)),
                                                        grep("-KIF11$", rownames(markerHM)),
                                                        grep("-TPX2$", rownames(markerHM)),
                                                        grep("-CCNA2$", rownames(markerHM)),
                                                        grep("-TOP2A$", rownames(markerHM)),
                                                        grep("-CCNE2$", rownames(markerHM)),
                                                        grep("-UNG$", rownames(markerHM)),
                                                        grep("-KPNA2$", rownames(markerHM)),
                                                        grep("-TUBA1C$", rownames(markerHM)),
                                                        grep("-BRCA1$", rownames(markerHM)),
                                                        grep("-CDKN1A$", rownames(markerHM))
                                                        ),
                                                 labels = c("FBXO5",
                                                            "POLE",
                                                            "CCND1",
                                                            "MCM2",
                                                            "MCM5",
                                                            "H3C2",
                                                            "H1-4",
                                                            "MKI67",
                                                            "CENPE",
                                                            "UBE2C",
                                                            "CDK1",
                                                            "CCNB1",
                                                            "KIF11",
                                                            "TPX2",
                                                            "CCNA2",
                                                            "TOP2A",
                                                            "CCNE2",
                                                            "UNG",
                                                            "KPNA2",
                                                            "TUBA1C",
                                                            "BRCA1",
                                                            "CDKN1A") ) )
cc_annotation <- rowAnnotation(cc = anno_mark(at = c(
                                                     grep("-PLK2$", rownames(markerHM)),
                                                     grep("-CDKN1A$", rownames(markerHM)),
                                                     grep("-PLK1$", rownames(markerHM)),
                                                     grep("-CCNB1$", rownames(markerHM)),
                                                     grep("-AURKA$", rownames(markerHM)),
                                                     grep("-CDC20$", rownames(markerHM)),
                                                     grep("-CCNA2$", rownames(markerHM)),
                                                     grep("-BUB1$", rownames(markerHM)),
                                                     grep("-CENPE$", rownames(markerHM)),
                                                     grep("-CCND1$", rownames(markerHM)),
                                                     grep("-PLK2$", rownames(markerHM)),
                                                     grep("-CDC6$", rownames(markerHM)),
                                                     grep("-CHAF1A$", rownames(markerHM)),
                                                     grep("-PCNA$", rownames(markerHM)),
                                                     grep("-POLE$", rownames(markerHM)),
                                                     grep("-FEN1$", rownames(markerHM)),
                                                     grep("-PRIM1$", rownames(markerHM)),
                                                     grep("-CHAF1B$", rownames(markerHM)),
                                                     grep("-CDC45$", rownames(markerHM)),
                                                     grep("-POLA1$", rownames(markerHM)),
                                                     grep("-GMNN$", rownames(markerHM)),
                                                     grep("-CDC25A$", rownames(markerHM)),
                                                     grep("-H3C2$", rownames(markerHM)),
                                                     grep("-H1-4$", rownames(markerHM)),
                                                     grep("-MKI67$", rownames(markerHM)),
                                                     grep("-CDKN2B$", rownames(markerHM)),
                                                     grep("-DNMT3A$", rownames(markerHM)),
                                                     grep("-RRM2B$", rownames(markerHM)),
                                                     grep("-CDKN2A$", rownames(markerHM)),
                                                     grep("-SGO1$", rownames(markerHM)),
                                                     grep("-FBXO5$", rownames(markerHM)),
                                                     grep("-CDC25B$", rownames(markerHM)),
                                                     grep("-CDC25C$", rownames(markerHM))
                                                     ),
                                              labels = c(
                                                     "PLK2",
                                                     "CDKN1A",
                                                     "PLK1",
                                                     "CCNB1",
                                                     "AURKA",
                                                     "CDC20",
                                                     "CCNA2",
                                                     "BUB1",
                                                     "CENPE",
                                                     "CCND1",
                                                     "PLK2",
                                                     "CDC6",
                                                     "CHAF1A",
                                                     "PCNA",
                                                     "POLE",
                                                     "FEN1",
                                                     "PRIM1",
                                                     "CHAF1B",
                                                     "CDC45",
                                                     "POLA1",
                                                     "GMNN",
                                                     "CDC25A",
                                                     "H3C2",
                                                     "H1-4",
                                                     "MKI67",
                                                     "CDKN2B",
                                                     "DNMT3A",
                                                     "RRM2B",
                                                     "CDKN2A",
                                                     "SGO1",
                                                     "FBXO5",
                                                     "CDC25B",
                                                     "CDC25C" 
                                                     ) ) )
Heatmap(as.matrix(markerHM),
        top_annotation = cell_annotation,
        right_annotation = cc_annotation,
        #col = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))), # log fold change
        #col = colorRamp2(seq(0, 5, length = 11), rev(brewer.pal(11, "RdBu"))),
        na_col = "grey0",
        name = "z-score expression",
        cluster_rows = TRUE,
        show_row_dend = FALSE,
        #clustering_distance_rows = function(x) as.dist(1-cosine_similarity(x)),
        clustering_distance_rows = "pearson",
        clustering_method_rows = "ward.D2",
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        #row_split = factor(treatment[rownames(selected_hm), "treatment"], levels = c("Rich","Arg3h", "Arg6h", "Leu3h", "Leu6h"), ordered = TRUE),
        #column_split = factor(rep(display_codons, each = length(levels(site_change$site))), levels = display_codons, ordered = TRUE),
        column_title_rot = 90,
        border = TRUE,
        use_raster = FALSE)
dev.off()

## example hm

asite_pausing <- site_change %>%
  filter(codon_display %in% highlight_codons) %>%
  filter(site %in% c("a")) %>%
  select(CB, codon_display, site_fc) %>%
  pivot_wider(names_from = codon_display, values_from = site_fc) %>%
  as.data.frame()
rownames(asite_pausing) <- asite_pausing$CB



display_gene <- "MYL6"

EXGcodons <- codons %>%
  filter(transcript_id == pull(annot %>% filter(gene_name == display_gene, set == "canonical") %>% select(transcript_id)))

EXG <- EXGcodons %>%
  inner_join(codon_reads %>% filter(site %in% c("a")), by = c("transcript_id" = "transcript_id", "codon_position" = "site_position")) %>%
  #inner_join(codon_reads %>% filter(site %in% c("p")) , by = c("transcript_id" = "transcript_id", "codon_position" = "site_position")) %>%
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
EXG_hm <- EXG_hm[pull(RPFmeta %>% filter(CB %in% rownames(EXG_hm)) %>% arrange(-pseudotime) %>% select(CB)), ]

aa_highlight <- rep("other", nrow(EXGcodons))
aa_highlight[grep("GAA", EXGcodons$codon)] <- "GAA"
#aa_highlight[grep("ATA", EXGcodons$codon)] <- "AUA"
#aa_highlight[grep("CGA", EXGcodons$codon)] <- "CGA"

pdf(paste0(figurePrefix, "_pausing_hm_", display_gene, ".pdf"), height = 3.5, width = 9)
#aa_annotation <- HeatmapAnnotation(aa = anno_text(EXGcodons$single, gp = gpar(fontsize = 2)))
aa_show <- HeatmapAnnotation(#an_text = anno_text(EXGcodons$single, gp = gpar(fontsize=1)),
                             aa = aa_highlight,
                             gaa_label = anno_mark(at = c(6,91),
                                                   labels = c("E6", "E91")),
                             col = list(aa = c("other" = "white",
                                               "GAA" = "red",
                                               "AUA" = "blue",
                                               "CGA" = "green")),
                             annotation_name_side = "left")
#aa_annotation <- HeatmapAnnotation(aa = as.character(rep(EXGcodons$sin
exprn <- as.numeric(PI@assays$RNA[ pull(annot %>% 
                                        filter(gene_name == display_gene) %>% 
                                        mutate(gene_out = paste(gene_id, gene_name, sep = '-')) %>% 
                                        select(gene_out) %>% 
                                        distinct()), rownames(EXG_hm) ] )
cell_annotation <- HeatmapAnnotation(exprn = exprn,
                                     GAA = log2(asite_pausing[rownames(EXG_hm), "GAA_Glu"]),
                                     #population = as.character(RPFmeta[rownames(EXG_hm), "population"]),
                                     clusters = as.character(RPFmeta[rownames(EXG_hm), "clusters"]),
                                     #pseudotime = RPFmeta[rownames(EXG_hm), "pseudotime"],
                                     #mKO2 = RPFmeta[rownames(EXG_hm), "log_red"],
                                     #mAG = RPFmeta[rownames(EXG_hm), "log_green"],
                                     col = list(clusters = c("0" = gg_color_hue(8)[1],
                                                             "1" = gg_color_hue(8)[2],
                                                             "2" = gg_color_hue(8)[3],
                                                             "3" = gg_color_hue(8)[4],
                                                             "4" = gg_color_hue(8)[5],
                                                             "5" = gg_color_hue(8)[6],
                                                             "6" = gg_color_hue(8)[7],
                                                             "7" = gg_color_hue(8)[8]),
                                                population = c("G0" = population_colours[1],
                                                               "Interphase" = population_colours[2],
                                                               "Mitotic" = population_colours[3]),
                                                mKO2 = colorRamp2(seq(0, max(RPFmeta$log_red), length.out = 8), brewer.pal(8, 'Reds')),
                                                mAG = colorRamp2(seq(0, max(RPFmeta$log_red), length.out = 8), brewer.pal(8, 'Greens')),
                                                GAA = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "RdBu"))),
                                                exprn = colorRamp2(seq(min(exprn), max(exprn), length.out = 6), (brewer.pal(6, "Blues"))),
                                                pseudotime = colorRamp2(seq(min(RPFmeta$pseudotime), max(RPFmeta$pseudotime), length.out = 11), rev(brewer.pal(11, 'Spectral')))),
                                     which = "row",
                                     annotation_name_side = "top")
length_annotation <- HeatmapAnnotation(pos = anno_mark(at = c(seq(25,ncol(EXG_hm), by = 50)),
                                                       labels = as.character(c(seq(25,ncol(EXG_hm), by = 50))), 
                                                       which = "column", side = "bottom", labels_rot = 0, link_height = unit(3, "mm")))
Heatmap(as.matrix(EXG_hm),
        #top_annotation = aa_annotation,
        top_annotation = aa_show,
        right_annotation = cell_annotation,
        bottom_annotation = length_annotation,
        #col = colorRamp2(seq(0,max(EXG_hm), length = 11), rev(brewer.pal(11, "RdYlBu"))),
        #col = colorRamp2(seq(0,max(EXG_hm), length = 9), brewer.pal(9, "Greys")),
        col = colorRamp2(seq(0,.5, length = 9), brewer.pal(9, "Greys")),
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
        #row_split = factor(argSortSplit[rownames(EXG_hm), "sortSplit"], levels = c("first", "second", "third", "fourth", "rest"), ordered = TRUE),
        row_title = NULL,
        column_title_rot = 0,
        border = TRUE,
        use_raster = FALSE)
dev.off()



## 
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

cluster_change <- codon_reads %>%
  filter(site == "a") %>%
  #filter(gene_id %in% pull(gene_codon_background %>% filter(gene_total > 2500) %>% select(gene_id) %>% distinct())) %>%
  filter(gene_id %in% pull(gene_codon_background %>% filter(gene_mean > 1) %>% select(gene_id) %>% distinct())) %>%
  inner_join(RPFmeta %>% group_by(clusters) %>% mutate(cells_per_cluster = n()) %>% ungroup() %>% select(CB, clusters, cells_per_cluster), by = "CB") %>%
  group_by(gene_id, gene_name, codon, clusters, cells_per_cluster) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(gene_name, clusters) %>%
  mutate(cluster_gene_frac = n/sum(n),
         cluster_gene_total = sum(n)) %>%
  ungroup() %>% 
  mutate(cluster_gene_mean = cluster_gene_total / cells_per_cluster) %>%
  inner_join(gene_codon_background %>% select(gene_id, codon, codon_frac), by = c("gene_id" = "gene_id", "codon" = "codon")) %>%
  mutate(lfc_change = log2(cluster_gene_frac/codon_frac))


shared_genes <- pull(cluster_change %>%
  group_by(gene_id, clusters) %>%
  summarize(asites_per_cluster = sum(n)) %>%
  ungroup() %>%
  complete(gene_id, clusters) %>%
  mutate(asites_per_cluster = replace_na(asites_per_cluster,0)) %>%
  inner_join(RPFmeta %>% group_by(clusters) %>% mutate(cells_per_cluster = n()) %>% ungroup() %>% select(clusters, cells_per_cluster) %>% distinct(), by = "clusters") %>%
  mutate(mean_asites = asites_per_cluster / cells_per_cluster) %>%
  filter(mean_asites > 2.5) %>%
  group_by(gene_id) %>%
  summarize(detected_in_clusters = n()) %>%
  ungroup() %>%
  filter(detected_in_clusters > 3) %>%
  select(gene_id) %>%
  distinct())



plt_cluster_change_gaa_facet <- cluster_change %>%
  filter(codon == "GAA") %>%
  filter(gene_id %in% shared_genes) %>%
  #mutate(gene_display = case_when(gene_name == "MYL6" ~ "MYL6",
  mutate(gene_display = case_when(gene_name == "MYL6" ~ "MYL6",
                                  TRUE ~ "")) %>%
  ggplot(aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = clusters, label = gene_display ))+
  geom_point(size=1)+
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~clusters, nrow=2)+
  geom_text_repel(segment.color = "black", box.padding = 0.4)+
  xlab("log10 mean cluster expression")+ylab("log2 f.c. GAA A-site FOC")+
  theme(strip.background = element_blank())
save_plot(paste0(figurePrefix, "_plt_cluster_change_GAA.pdf"), plt_cluster_change_gaa_facet, base_width = 6, base_height =3.25)


xlims <- log10(c(min(pull(cluster_change %>% filter(gene_id %in% shared_genes) %>% filter(codon %in% c("GAA","CAG", "GAG", "ATA", "CGA")) %>% select(cluster_gene_mean))),
                 max(pull(cluster_change %>% filter(gene_id %in% shared_genes) %>% filter(codon %in% c("GAA","CAG", "GAG", "ATA", "CGA")) %>% select(cluster_gene_mean)))))
ylims <- c(min(pull(cluster_change %>% filter(gene_id %in% shared_genes) %>% filter(codon %in% c("GAA","CAG", "GAG", "ATA", "CGA")) %>% mutate(l2fc_frac_change = log2(cluster_gene_frac/codon_frac)) %>% select(l2fc_frac_change))),
           max(pull(cluster_change %>% filter(gene_id %in% shared_genes) %>% filter(codon %in% c("GAA","CAG", "GAG", "ATA", "CGA")) %>% mutate(l2fc_frac_change = log2(cluster_gene_frac/codon_frac)) %>% select(l2fc_frac_change))))

plt_cluster_change_cag_facet <- cluster_change %>%
  filter(codon == "CAG") %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name == "MYL6" ~ "MYL6",
                                  TRUE ~ "")) %>%
  ggplot(aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = clusters, label = gene_display ))+
  geom_point()+
  coord_cartesian(xlim = xlims, ylim = ylims)+ 
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~clusters, nrow=1)+
  geom_text_repel(segment.color = "black", box.padding = 0.4)+
  xlab("log10 mean cluster expression")+ylab("log2 f.c. CAG A-site FOC")+
  theme(strip.background = element_blank())
plt_cluster_change_gag_facet <- cluster_change %>%
  filter(codon == "GAG") %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name == "MYL6" ~ "MYL6",
                                  TRUE ~ "")) %>%
  ggplot(aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = clusters, label = gene_display ))+
  geom_point()+
  coord_cartesian(xlim = xlims, ylim = ylims)+ 
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~clusters, nrow=1)+
  geom_text_repel(segment.color = "black", box.padding = 0.4)+
  xlab("log10 mean cluster expression")+ylab("log2 f.c. GAG A-site FOC")+
  theme(strip.background = element_blank())
plt_cluster_change_gaa_facet <- cluster_change %>%
  filter(codon == "GAA") %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name == "MYL6" ~ "MYL6",
                                  TRUE ~ "")) %>%
  ggplot(aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = clusters, label = gene_display ))+
  geom_point()+
  coord_cartesian(xlim = xlims, ylim = ylims)+ 
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~clusters, nrow=1)+
  geom_text_repel(segment.color = "black", box.padding = 0.4)+
  xlab("log10 mean cluster expression")+ylab("log2 f.c. GAA A-site FOC")+
  theme(strip.background = element_blank())
plt_cluster_change_aua_facet <- cluster_change %>%
  filter(codon == "ATA") %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name == "MYL6" ~ "MYL6",
                                  TRUE ~ "")) %>%
  ggplot(aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = clusters, label = gene_display ))+
  geom_point()+
  coord_cartesian(xlim = xlims, ylim = ylims)+ 
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~clusters, nrow=1)+
  geom_text_repel(segment.color = "black", box.padding = 0.4)+
  xlab("log10 mean cluster expression")+ylab("log2 f.c. ATA A-site FOC")+
  theme(strip.background = element_blank())
plt_cluster_change_cga_facet <- cluster_change %>%
  filter(codon == "CGA") %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name == "MYL6" ~ "MYL6",
                                  TRUE ~ "")) %>%
  ggplot(aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = clusters, label = gene_display ))+
  geom_point()+
  coord_cartesian(xlim = xlims, ylim = ylims)+ 
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~clusters, nrow=1)+
  geom_text_repel(segment.color = "black", box.padding = 0.4)+
  xlab("log10 mean cluster expression")+ylab("log2 f.c. CGA A-site FOC")+
  theme(strip.background = element_blank())
save_plot(paste0(figurePrefix, "_plt_cluster_change_rest.pdf"), plot_grid(plt_cluster_change_cag_facet, plt_cluster_change_gaa_facet, plt_cluster_change_gag_facet, plt_cluster_change_aua_facet, plt_cluster_change_cga_facet, ncol = 1), base_width = 14, base_height = 2.5*5)


plt_cluster_box <- cluster_change %>%
  #filter(codon %in% c("GAA")) %>%
  filter(codon %in% c("GAA","CAG", "GAG", "ATA", "CGA")) %>%
  filter(gene_id %in% shared_genes) %>%
  #filter(codon %in% c("GAA", "ATA", "CGA", "GAG")) %>%
  mutate(gene_display = case_when(gene_name == "MYL6" ~ "MYL6",
                                  TRUE ~ "")) %>%
  ggplot(aes(x=clusters, y=log2(cluster_gene_frac/codon_frac), colour = clusters ))+
  geom_boxplot()+
  facet_wrap(~codon)#+
  #geom_text_repel(segment.color = "black", box.padding = 0.4)
save_plot(paste0(figurePrefix, "_plt_cluster_box.pdf"), plt_cluster_box, base_width = 5, base_height =5)




## GAA gene-wise pausing stats
## Pausing statistics


GAA_a_genes <- codon_reads %>%
  filter(site == "a") %>%
  filter(gene_id %in% shared_genes) %>%
  group_by(gene_id, gene_name, codon, CB) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(CB, gene_name) %>%
  mutate(cell_gene_total = sum(n)) %>%
  ungroup() %>%
  mutate(cell_gene_codon_frac = n/cell_gene_total) %>%
  filter(codon == "GAA") %>%
  mutate(gene_display = paste(gene_id, gene_name, sep = "-")) %>%
  select(gene_display, CB, cell_gene_codon_frac) %>%
  complete(gene_display, CB) %>%
  mutate(cell_gene_codon_frac = replace_na(cell_gene_codon_frac, 0)) %>%
  pivot_wider(names_from = CB, values_from = cell_gene_codon_frac) %>%
  as.data.frame()
rownames(GAA_a_genes) <- GAA_a_genes$gene_display
GAA_a_genes <- GAA_a_genes[, colnames(PI@assays$RNA)]

PIp <- PI
PIp[['GAA_apausing']] <- CreateAssayObject(data = as.matrix(GAA_a_genes))
gaa_a_pausing_markers <- FindAllMarkers(PIp, assay = "GAA_apausing", logfc.threshold = 0, only.pos = FALSE, pseudocount.use = 0, min.cells.feature = 10, min.cells.group =0) %>%
  mutate(cluster_redef = case_when(cluster == 5 ~ 0,
                              cluster == 2 ~ 1,
                              cluster == 0 ~ 2,
                              cluster == 3 ~ 3,
                              cluster == 4 ~ 4,
                              cluster == 7 ~ 5,
                              cluster == 1 ~ 6,
                              cluster == 6 ~ 7,
                              TRUE ~ -1)) %>%
  mutate(avg_logFC = log2(exp(avg_logFC)))


cl6_sig_genes <- gaa_a_pausing_markers %>%
  filter(p_val_adj < 0.01) %>%
  filter(avg_logFC > 0) %>%
  filter(cluster_redef == 6) %>%
  select(gene) %>%
  distinct()
pct_sig <- nrow(cl6_sig_genes)/sum(apply(GAA_a_genes > 0, 1, sum) > 10)



## UGC stats

ugc_stats <- site_change %>%
  filter(codon_display == "UGC_Cys") %>%
  inner_join(RPFmeta %>% select(CB, clusters), by="CB") %>%
  mutate(grouping_cluster = case_when(clusters %in% c(2,7) ~ "2_7",
                                      TRUE ~ "rest")) %>%
  group_by(grouping_cluster) %>%
  summarize(mean_site_fraction = mean(100*cell_site_fraction),
            sd_site_fraction = sd(100*cell_site_fraction)) %>%
  ungroup()

ugc_stats <- codon_reads %>%
  group_by(CB, codon) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(CB) %>%
  mutate(cell_total = sum(n)) %>%
  ungroup() %>%
  mutate(cell_codon_fraction = n/cell_total) %>%
  inner_join(RPFmeta %>% select(CB, clusters), by = "CB") %>%
  filter(codon == "TGC") %>%
  filter(clusters %in% c(2,7,6)) %>%
  mutate(grouping_cluster = case_when(clusters %in% c(2,7) ~ "2_7",
                                      TRUE ~ "rest")) %>%
  group_by(grouping_cluster) %>%
  summarize(mean_ugc = mean(100*cell_codon_fraction),
            sd_ugc = sd(100*cell_codon_fraction)) %>%
  ungroup()


## GAA a-site stats

gaa_stats <- codon_reads %>%
  filter(site == "a") %>%
  group_by(CB, codon) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(CB) %>%
  mutate(cell_total = sum(n)) %>%
  ungroup() %>%
  mutate(cell_codon_fraction = n/cell_total) %>%
  inner_join(RPFmeta %>% select(CB, clusters), by = "CB") %>%
  filter(codon == "GAA") %>%
  mutate(grouping_cluster = case_when(clusters == 6 ~ 6,
                                      TRUE ~ -1)) %>%
  group_by(grouping_cluster) %>%
  summarize(mean_gaa = mean(100*cell_codon_fraction),
            sd_gaa = sd(100*cell_codon_fraction)) %>%
  ungroup()


## MYL6 GAA a-site stats for text
myl6_gaa_stats <- codon_reads %>%
  filter(site == "a",
         gene_name == "MYL6") %>%
  group_by(CB, codon) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  group_by(CB) %>%
  mutate(cell_total = sum(n)) %>%
  ungroup() %>%
  mutate(cell_codon_fraction = n/cell_total) %>%
  inner_join(RPFmeta %>% select(CB, clusters), by = "CB") %>%
  filter(codon == "GAA") %>%
  mutate(grouping_cluster = case_when(clusters == 6 ~ 6,
                                      TRUE ~ -1)) %>%
  group_by(grouping_cluster) %>%
  summarize(mean_gaa = mean(100*cell_codon_fraction),
            sd_gaa = sd(100*cell_codon_fraction)) %>%
  ungroup()

## GAA a-stie average increase
cl6_increase <- cluster_change %>%
  filter(codon == "GAA") %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(grouping_cluster = case_when(clusters == 6 ~ 6,
                                      TRUE ~ -1)) %>%
  group_by(grouping_cluster) %>%
  summarize(mean_gaa = mean(cluster_gene_frac/codon_frac),
            sd_gaa = sd(cluster_gene_frac/codon_frac)) %>%
  ungroup()


### pst-prediction-genes
selected_genes <- c("ENSG00000286522.1-H3C2", "ENSG00000100297.16-MCM5", "ENSG00000124762.13-CDKN1A", "ENSG00000138778.13-CENPE", "ENSG00000148773.14-MKI67", "ENSG00000131747.15-TOP2A", "ENSG00000112029.10-FBXO5")
selexp <- as.data.frame(PI@assays$integrated[ selected_genes , pull(RPFmeta2 %>% arrange(pseudotime) %>% select(CB)) ])
selexp$gene <- rownames(selexp)
selexp <- selexp %>% 
  pivot_longer(!gene, names_to = "CB", values_to = "expression") %>%
  inner_join(RPFmeta %>% select(CB, pseudotime, clusters), by = "CB")


Pmeta <- data.frame(PB = paste0("V", rep(1:nrow(prediction_mean))),
                    pseudotime = seq(from = min_pst, to = max_pst, length.out = nrow(prediction_mean)))
selpred <- as.data.frame(t(prediction_mean[,selected_genes]))
selpred$gene <- rownames(selpred)
selpred <- selpred %>%
  pivot_longer(!gene, names_to = "PB", values_to = "expression") %>%
  inner_join(Pmeta %>% select(PB, pseudotime), by = "PB")



plt_selsca <- ggplot(selexp, aes(x=pseudotime, y=expression, colour = gene))+
  geom_point()+
  geom_line(data = selpred)+
  facet_wrap(~gene, ncol=1, scales = "free_y")
save_plot(paste0(figurePrefix, "_selsca.pdf"), plt_selsca, base_width = 10, base_height = 10)


### codon usage

cell_counts <- RPFreads %>%
  filter(gene_type == "protein_coding",
         set == "canonical") %>%
  group_by(CB, transcript_id) %>%
  summarize(n_reads = n()) %>%
  ungroup()

# cell_codon_counts <- do.call('rbind', lapply(unique(cell_counts$CB), function(cb) {
#                                          codon_counts <- cell_counts %>%
#                                            filter(CB == cb) %>%
#                                            inner_join(codons %>% filter(nchar(codon) == 3) %>% select(transcript_id, codon), by = "transcript_id") %>%
#                                            group_by(CB, codon) %>%
#                                            summarize(one = n(),
#                                                      sep = sum(n_reads), .groups = "drop_last") %>%
#                                            ungroup()
#                                   } ) )
# saveRDS(cell_codon_counts, paste0(figurePrefix, "_cell_codon_counts.RDS"))

cell_counts_if <- RPFreads %>%
  filter(gene_type == "protein_coding",
         set == "canonical") %>%
  mutate(frameP = (psite - cds_start) %% 3) %>%
  filter(frameP == 0) %>%
  group_by(CB, transcript_id) %>%
  summarize(n_reads = n()) %>%
  ungroup()

codons_filtered <- codons %>%
  filter(nchar(codon) == 3) %>%
  filter(!(three %in% c("Ter"))) %>%
  mutate(codon_display = paste(gsub("T", "U", codon), three, sep = "_")) %>%
  select(transcript_id, codon_display)

# cell_codon_counts_if <- do.call('rbind', lapply(unique(cell_counts_if$CB), function(cb) {
#                                             codon_counts <- cell_counts_if %>%
#                                               filter(CB == cb) %>%
#                                               inner_join(codons_filtered, by = "transcript_id") %>%
#                                               group_by(CB, codon_display) %>%
#                                               summarize(one = n(),
#                                                         sep = sum(n_reads), .groups = "drop_last") %>%
#                                               ungroup()
#                                   } ) )
# saveRDS(cell_codon_counts_if, paste0(figurePrefix, "_cell_codon_counts_if.RDS"))
cell_codon_counts_if <- readRDS(paste0(figurePrefix, "_cell_codon_counts_if.RDS"))

cell_codon_site <- codon_reads %>%
  filter(nchar(codon) == 3) %>%
  filter(!(three %in% c("Ter"))) %>%
  mutate(codon_display = paste(gsub("T", "U", codon), three, sep = "_")) %>%
  group_by(CB, site, codon_display) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(CB, site) %>%
  mutate(cell_site_fraction = n/sum(n)) %>%
  ungroup() %>%
  group_by(site) %>%
  mutate(site_total = sum(n)) %>%
  ungroup() %>%
  group_by(codon_display, site) %>%
  mutate(baseline_site_fraction = sum(n)/site_total) %>%
  ungroup() %>%
  select(-site_total)



cell_codon_usage <- cell_codon_counts_if %>%
  pivot_longer(cols = c('one', 'sep'), names_to = "site", values_to = "n") %>%
  group_by(CB, site) %>%
  mutate(cell_site_fraction = n/sum(n)) %>%
  ungroup() %>%
  group_by(site) %>%
  mutate(site_total = sum(n)) %>%
  ungroup() %>%
  group_by(codon_display, site) %>%
  mutate(baseline_site_fraction = sum(n)/site_total) %>%
  ungroup() %>%
  select(-site_total) %>%
  # rbind(site_change %>% select(CB, codon_display, site, n, cell_site_fraction, baseline_site_fraction)) %>%
  rbind(cell_codon_site %>% select(CB, codon_display, site, n, cell_site_fraction, baseline_site_fraction)) %>%
  inner_join(RPFmeta %>% select(CB, pseudotime, clusters), by = "CB") %>%
  mutate(site_order = case_when(site == "em2" ~ 1,
                                site == "em1" ~ 2,
                                site == "e" ~ 3,
                                site == "p" ~ 4,
                                site == "a" ~ 5,
                                site == "ap1" ~ 6,
                                site == "ap2" ~ 7,
                                site == "one" ~ 8,
                                site == "sep" ~ 9,
                                TRUE ~ 100)) %>%
  mutate(site = fct_reorder(site, site_order)) 



#plt_codon_frequency <- ggplot(cell_codon_usage %>% filter(codon_display %in% c("UGC_Cys")) , aes(x=pseudotime, y=log2(cell_site_fraction/baseline_site_fraction), colour = site))+



pltlist <- list()
for(i in 1:length(selected_codons)){
  pltlist[[i]] <- ggplot(cell_codon_usage %>% filter(codon_display %in% selected_codons[i]) , aes(x=pseudotime, y=cell_site_fraction, colour = clusters))+
  geom_point(size=0.5)+
  geom_smooth(colour = "black", se = TRUE, method = "loess")+
  facet_wrap(~codon_display+site, ncol = 9)+
  theme(legend.position = "none")
}
save_plot(paste0(figurePrefix, "_plt_codon_frequency_usage.pdf"), plot_grid(plotlist = pltlist, nrow = length(selected_codons), align = "hv"), base_width = 13, base_height = 2.0*length(selected_codons))

pltlist <- list()
for(i in 1:length(selected_codons)){
  pltlist[[i]] <- ggplot(cell_codon_usage %>% filter(codon_display %in% selected_codons[i]) , aes(x=pseudotime, y=log2(cell_site_fraction/baseline_site_fraction), colour = clusters))+
  geom_point(size=0.5)+
  geom_smooth(colour = "black", se = TRUE, method = "loess")+
  facet_wrap(~codon_display+site, ncol = 9)+
  theme(legend.position = "none")
}
save_plot(paste0(figurePrefix, "_plt_codon_frequency_usage_lfc.pdf"), plot_grid(plotlist = pltlist, nrow = length(selected_codons)), base_width = 17, base_height = 3*length(selected_codons))


sanity <- cell_codon_usage %>%
  group_by(CB, site) %>%
  mutate(cell_site_total = sum(n)) %>%
  ungroup()



## GEO output
# KO/SR2 whitelists
KO_wells <- c(columnFeatures(17:24, LETTERS[1:16]), #WT
              columnFeatures(22:24, LETTERS[15:16])) #NTC
SR2_wells <- c(columnFeatures(19:21, LETTERS[1:16]), # cycling_post
               columnFeatures(22:24, LETTERS[1:14]), # cycling_pre
               columnFeatures(22:24, LETTERS[15:16])) #NTC
KO_whitelist <- barcodes %>%
  filter(well %in% KO_wells)
SR2_whitelist <- barcodes %>%
  filter(well %in% SR2_wells)

KO_blacklist <- barcodes %>%
  filter(!(well %in% KO_wells)) %>%
  mutate(category = "black") %>%
  select(CB, category) 

SR2_blacklist <- barcodes %>%
  filter(!(well %in% SR2_wells)) %>%
  mutate(category = "black") %>%
  select(CB, category) 

write.table(KO_blacklist, file = paste0(figurePrefix, "_KO_blacklist.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")
write.table(SR2_blacklist, file = paste0(figurePrefix, "_SR2_blacklist.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# write.csv(KO_whitelist, file = paste0(figurePrefix, "_KO_whitelist.csv"), row.names = FALSE, quote = FALSE)
# write.csv(SR2_whitelist, file = paste0(figurePrefix, "_SR2_whitelist.csv"), row.names = FALSE, quote = FALSE)
write.csv(RPFmeta, file = paste0(figurePrefix, "_RPE-FUCCI_metadata.csv"))
write.csv(RPFcounts, file = paste0(figurePrefix, "_RPE-FUCCI_RPFcounts.csv"))


## Plate layouts
First_G0 <- barcodes
First_G0$fraction <- "Interphase"
First_G0[First_G0$well %in% columnFeatures(1:8, LETTERS[1:16]), "fraction"] = "G0"
First_G0[First_G0$well %in% columnFeatures(21:24, LETTERS[15:16]), "fraction"] = "NTC"
write.csv(First_G0, file = paste0(figurePrefix, "_RPE_G0_layout.csv"), row.names = FALSE, quote = FALSE)

First_Mit <- barcodes
First_Mit$fraction <- "Interphase"
First_Mit[First_Mit$well %in% columnFeatures(1:8, LETTERS[1:16]), "fraction"] = "Mitotic_shake_off"
First_Mit[First_Mit$well %in% columnFeatures(21:24, LETTERS[15:16]), "fraction"] = "NTC"
write.csv(First_Mit, file = paste0(figurePrefix, "_RPE_Mit_layout.csv"), row.names = FALSE, quote = FALSE)


KO <- barcodes
KO$fraction <- "Interphase"
KO[ KO$well %in% columnFeatures(1:4, LETTERS[1:16]) , "fraction"] = "KO17"
KO[ KO$well %in% columnFeatures(5:8, LETTERS[1:16]) , "fraction"] = "KO15"
KO[ KO$well %in% columnFeatures(9:12, LETTERS[1:16]) , "fraction"] = "KO14"
KO[ KO$well %in% columnFeatures(13:16, LETTERS[1:16]) , "fraction"] = "KO13"
KO[ KO$well %in% columnFeatures(17:24, LETTERS[1:16]) , "fraction"] = "Interphase"
KO[ KO$well %in% columnFeatures(22:24, LETTERS[15:16]) , "fraction"] = "NTC"
KO <- KO[KO$fraction %in% c("Interphase", "NTC"),]
write.csv(KO, file = paste0(figurePrefix, "_RPE_KO_layout.csv"), row.names = FALSE, quote = FALSE)

SR2 <- barcodes
SR2[ SR2$well %in% c(columnFeatures(1:4, LETTERS[1:16]), columnFeatures(5, LETTERS[1:8])) , "fraction"] = "2h_post"
SR2[ SR2$well %in% c(columnFeatures(6:9, LETTERS[1:16]), columnFeatures(5, LETTERS[9:16])) , "fraction"] = "6h_post"
SR2[ SR2$well %in% c(columnFeatures(10:13, LETTERS[1:16]), columnFeatures(14, LETTERS[1:8])) , "fraction"] = "10h_post"
SR2[ SR2$well %in% c(columnFeatures(15:18, LETTERS[1:16]), columnFeatures(14, LETTERS[9:16])) , "fraction"] = "20h_post"
SR2[ SR2$well %in% columnFeatures(19:21, LETTERS[1:16]) , "fraction"] = "Interphase"
SR2[ SR2$well %in% columnFeatures(22:24, LETTERS[1:14]) , "fraction"] = "Interphase"
SR2[ SR2$well %in% columnFeatures(22:24, LETTERS[15:16]) , "fraction"] = "NTC"
SR2 <- SR2[SR2$fraction %in% c("Interphase", "NTC"),]
write.csv(SR2, file = paste0(figurePrefix, "_RPE_SR2_layout.csv"), row.names = FALSE, quote = FALSE)



