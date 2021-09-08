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
library(dendsort)
theme_set(theme_cowplot())

source('plotFunctions.R')

data_dir <- "../scRiboSeq/data"

figurePrefix <- "replicate_rpe/F3"

# barcodes
barcodes <- read.table('/hpc/hub_oudenaarden/mvanins/local/lib/barcodes/miRv6_barcodemap.tsv',
                       sep="\t",
                       header = FALSE,
                       stringsAsFactors = FALSE,
                       row.names = NULL,
                       col.names = c("CB","number","well"))
rownames(barcodes) <- barcodes$CB

# annot
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



RPFmeta <- readRDS(file.path(data_dir, "Fucci/compiled/meta.rds"))
RPFcounts <- readRDS(file.path(data_dir, "Fucci/compiled/counts.rds"))
RPFbiotypes <- readRDS(file.path(data_dir, "Fucci/compiled/biotypes.rds"))
RPFreads <- readRDS(file.path(data_dir, "Fucci/compiled/reads.rds"))




ENSG_to_name <- data.frame(count_names = rownames(RPFcounts), stringsAsFactors = FALSE) %>%
  inner_join(annot, by = c("count_names" = "gene_id")) %>% 
  mutate(gene_out = paste(count_names, gene_name, sep = "-")) %>%
  select(count_names, gene_out) %>%
  distinct() %>%
  column_to_rownames(var = "count_names")

rownames(RPFcounts) <- ENSG_to_name[rownames(RPFcounts),"gene_out"]
RPFcounts <- RPFcounts[!grepl("-PCDH", rownames(RPFcounts)),]


## QC check
plt_cd_totalQC <- ggplot(RPFmeta %>% filter(cell_total > 100), aes(x=log10(cell_total), y=cds_frac))+
  geom_point()+
  geom_vline(xintercept = log10(4500), colour = 'red', linetype = 'dashed')+
  geom_vline(xintercept = log10(50000), colour = 'red', linetype = 'dashed')+
  geom_hline(yintercept = .87, colour = 'red', linetype = 'dashed')+
  facet_wrap(~plate, ncol = 7)
save_plot(paste0(figurePrefix, "_plt_QC.pdf"), plt_cd_totalQC, base_width = 17, base_height = 9)

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
tours <- readRDS('trajectory/RPF_replicate_trajectory_20210413_all_raw_tours_5e+05iter_30pts_1noise.RDS') %>%
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


plt_jackstraw <- JackStrawPlot(PI,dims=1:20)
plt_elbow <- ElbowPlot(PI, ndims=20)
save_plot(paste0(figurePrefix,"_qc_plt_dimscores.pdf"), plot_grid(plt_elbow, plt_jackstraw, nrow=1), base_width = 15, base_height = 5)


ndims <- 16
PI <- RunUMAP(PI, dims=1:ndims, verbose=FALSE, n.neighbors = 30, n.epochs = 500, return.model = TRUE, seed.use = 112)
PI <- FindNeighbors(PI, dims=1:ndims, k.param=15)
PI <- FindClusters(PI,  n.start=50, n.iter=50)

RPFmeta$seurat_clusters <- PI$seurat_clusters[rownames(RPFmeta)]
RPFmeta$percent_mt <- PI$percent_mt[rownames(RPFmeta)]
RPFmeta$UMAP_1 <- as.data.frame(Embeddings(PI, reduction="umap"))[rownames(RPFmeta),"UMAP_1"]
RPFmeta$UMAP_2 <- as.data.frame(Embeddings(PI, reduction="umap"))[rownames(RPFmeta),"UMAP_2"]

saveRDS(RPFmeta, file.path(data_dir,"Fucci/compiled/clustered_meta.rds"))


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


plt_umap_clusters <- ggplot(RPFmeta[sample(nrow(RPFmeta), replace=FALSE),], aes(x=UMAP_1, y=UMAP_2, colour = seurat_clusters))+
  geom_point(size=1)+
  coord_equal() +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_population <- ggplot(RPFmeta, aes(x=UMAP_1, y=UMAP_2, colour = sort_population))+
  geom_point(size=1)+
  #scale_colour_manual(values = population_colours)+
  coord_equal()+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_batch <- ggplot(RPFmeta, aes(x=UMAP_1, y=UMAP_2, colour = batch))+
  geom_point(size=1)+
  #scale_colour_manual(values = population_colours)+
  coord_equal()+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_red <- ggplot(RPFmeta %>% arrange(log_red), aes(x=UMAP_1, y=UMAP_2, colour = log_red))+
  geom_point(size=1)+
  coord_equal()+
  scale_colour_gradientn(colours = brewer.pal(8,"Reds"), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_green <- ggplot(RPFmeta %>% arrange(log_green), aes(x=UMAP_1, y=UMAP_2, colour = log_green))+
  geom_point(size=1)+
  coord_equal()+
  scale_colour_gradientn(colours = brewer.pal(8,"Greens"), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_dapi <- ggplot(RPFmeta %>% arrange(X..405..460.50), aes(x=UMAP_1, y=UMAP_2, colour = X..405..460.50))+
  geom_point(size=1)+
  coord_equal()+
  scale_colour_gradientn(colours = brewer.pal(8,"Blues"), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_fsc <- ggplot(RPFmeta %>% arrange(FSC), aes(x=UMAP_1, y=UMAP_2, colour = FSC))+
  geom_point(size=1)+
  coord_equal()+
  scale_colour_gradientn(colours = brewer.pal(8,"Purples"), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_ssc <- ggplot(RPFmeta %>% arrange(SSC), aes(x=UMAP_1, y=UMAP_2, colour = SSC))+
  geom_point(size=1)+
  coord_equal()+
  scale_colour_gradientn(colours = brewer.pal(8,"Oranges"), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_rpf_slope <- ggplot(RPFmeta %>% arrange(rpf_slope), aes(x=UMAP_1, y=UMAP_2, colour = rpf_slope))+
  geom_point(size=1)+
  coord_equal()+
  scale_colour_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_cds <- ggplot(RPFmeta %>% arrange(cds), aes(x=UMAP_1, y=UMAP_2, colour = log10(cds)))+
  geom_point(size=1)+
  coord_equal()+
  scale_colour_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_cds_frac <- ggplot(RPFmeta %>% arrange(cds_frac), aes(x=UMAP_1, y=UMAP_2, colour = cds_frac))+
  geom_point(size=1)+
  coord_equal()+
  scale_colour_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

save_plot(paste0(figurePrefix, "_umap_plots.pdf"), plot_grid(plt_umap_clusters, plt_umap_population, plt_umap_batch,
                                                             plt_umap_red, plt_umap_green,
                                                             plt_umap_dapi, plt_umap_fsc, plt_umap_ssc,
                                                             plt_umap_rpf_slope, plt_umap_cds, plt_umap_cds_frac, nrow=4),
          base_width = 15, base_height = 15)


plt_umap_plate <- ggplot(RPFmeta[sample(nrow(RPFmeta), replace=FALSE),], aes(x=UMAP_1, y=UMAP_2, colour = plate))+
  geom_point(data = RPFmeta %>% select(UMAP_1, UMAP_2), colour = "grey50")+
  geom_point(size=1)+
  coord_equal() +
  facet_wrap(~plate)+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
save_plot(paste0(figurePrefix, "_umap_plate.pdf"), 
          plt_umap_plate,
          base_width = 15, base_height = 15)


## export integrated info for GP pseudotime
export_cells <- pull(RPFmeta %>%
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


## read pseudotimes from GrandPrix, run separately in pyhon notebook
# pst_dir <- "/hpc/hub_oudenaarden/mvanins/RPF/analysis/CS2integration/pausing_fig/prediction_nolowS/"
pst_dir <- "pseudotimes_replicateU/"

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
plt_pst_tsp <- ggplot(RPFmeta2, aes(x=path, y=pseudotime, colour = seurat_clusters))+
  geom_point()
save_plot(paste0(figurePrefix, "_pseudotime_tsp_comparison.pdf"), plt_pst_tsp, base_width = 6, base_height=6)
plt_umap_pst <- ggplot(RPFmeta2, aes(x=UMAP_1, y=UMAP_2, colour = pseudotime))+
  geom_point()+
  coord_equal()+
  scale_colour_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), na.value = "grey40")
plt_umap_path <- ggplot(RPFmeta2, aes(x=UMAP_1, y=UMAP_2, colour = path))+
  geom_point()+
  coord_equal()+
  scale_colour_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), na.value = "grey40")
save_plot(paste0(figurePrefix, "_pseudotime_umap.pdf"), plot_grid(plt_umap_path, plt_umap_pst, nrow=1), base_width = 16, base_height=6)

plt_pst_g <- ggplot(RPFmeta2, aes(x=pseudotime, y=log_green))+
  geom_point()
plt_pst_r <- ggplot(RPFmeta2, aes(x=pseudotime, y=log_red))+
  geom_point()
save_plot(paste0(figurePrefix, "pst_rgfacs.pdf"), plot_grid(plt_pst_r, plt_pst_g, ncol=1), base_width=5, base_height=5)


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

population_colours <- c("#f781bf", "#549ed1","#984ea3") # G0, Interphase, Mit

## reorder clusters
RPFmeta$clusters <- factor(pull(RPFmeta %>%
                                mutate(clusters = case_when(seurat_clusters == 4 ~ 0,
                                                            seurat_clusters == 6 ~ 1,
                                                            seurat_clusters == 2 ~ 7,
                                                            seurat_clusters == 1 ~ 2,
                                                            seurat_clusters == 3 ~ 3,
                                                            seurat_clusters == 5 ~ 4,
                                                            seurat_clusters == 7 ~ 5,
                                                            seurat_clusters == 0 ~ 6,
                                                            TRUE ~ -1)) %>%
                                select(clusters)))

### FACS plots
set.seed(112)
plt_facs_clusters <- ggplot(RPFmeta[sample(nrow(RPFmeta), replace=FALSE),], aes(x=log_green, y=log_red, colour = clusters))+
  geom_point(size = 0.5)+
  xlab("log10 mAG-GMNN")+
  ylab("log10 mKO2-CDT1")+
  coord_equal() #+
#   theme(axis.line = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
plt_facs_population <- ggplot(RPFmeta %>% arrange(-population_order), aes(x=log_green, y=log_red, colour = population))+
  geom_point(size = 0.5)+
  coord_equal()+
  xlab("log10 mAG-GMNN")+
  ylab("log10 mKO2-CDT1")+
  scale_colour_manual(values = population_colours) #+
#   theme(axis.line = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank())
plt_facs_pseudotime_proj <- ggplot(F_mean, aes(x=log_green, y=log_red, colour = pseudotime))+
  geom_point(size = 0.5, data = RPFmeta %>% select(log_green, log_red), colour = "grey50")+
  geom_point(size = 0.5)+
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
  geom_point(size = 0.5)+
  coord_equal() +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_population <- ggplot(RPFmeta %>% arrange(-population_order), aes(x=UMAP_1, y=UMAP_2, colour = population))+
  geom_point(size = 0.5)+
  scale_colour_manual(values = population_colours)+
  coord_equal()+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_pseudotime <- ggplot(RPFmeta2 %>% arrange(pseudotime), aes(x=UMAP_1, y=UMAP_2, colour = pseudotime))+
  geom_point(size = 0.5)+
  coord_equal()+
  scale_colour_gradientn(colours = rev(brewer.pal(11,"Spectral")), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_pseudotime_path <- ggplot(u_t, aes(x=UMAP_1, y=UMAP_2, colour = pseudotime))+
  geom_point(size = 0.5, data = RPFmeta %>% select(UMAP_1, UMAP_2), colour = "grey50", size=1)+
  geom_point(size = 0.5)+
  scale_colour_gradientn(colours = rev(brewer.pal(11,"Spectral")), na.value = "grey40")+
  coord_equal()+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_red <- ggplot(RPFmeta[sample(nrow(RPFmeta), replace=FALSE),], aes(x=UMAP_1, y=UMAP_2, colour = log_red))+
  geom_point(size = 0.5)+
  coord_equal()+
  scale_colour_gradientn(colours = brewer.pal(8,"Reds"), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
plt_umap_green <- ggplot(RPFmeta[sample(nrow(RPFmeta), replace=FALSE),], aes(x=UMAP_1, y=UMAP_2, colour = log_green))+
  geom_point(size = 0.5)+
  coord_equal()+
  scale_colour_gradientn(colours = brewer.pal(8,"Greens"), na.value = "grey40")+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

save_plot(paste0(figurePrefix, "_reduction_panels.pdf"), plot_grid(plt_facs_population, plt_facs_clusters, plt_facs_pseudotime_proj,
                                                                   plt_umap_population, plt_umap_clusters, plt_umap_pseudotime_path,
                                                                   plt_umap_green, plt_umap_red, ncol = 3 ), base_width = 10, base_height = 10)





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
  inner_join(RPFmeta %>% select(CB, red, green, clusters, seurat_clusters, sort_population, plate, pseudotime), by = "CB") %>%
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
pausing_markers <- FindAllMarkers(PIp, assay = "pausing_site", logfc.threshold = log(2^.25), only.pos = FALSE, pseudocount.use = 0) %>%
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
                  #filter(abs(avg_logFC) > .25) %>%
                  filter(p_val_adj < 1E-5) %>%
                  #filter(p_val_adj < 0.01) %>%
                  select(codon_display) %>%
                  distinct())


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


selected_to_quant <- site_change %>%
  filter(site == "a") %>%
  select(CB, site, site_fc) %>%
  mutate(l2fc = log2(site_fc))
quantile(abs(selected_to_quant$l2fc), probs = 1-0.01)

pausing_limit <- 1.0

selected_to_winzorize <- site_change %>%
  filter(codon_display %in% selected_codons) %>%
  filter(site == "a") %>%
  ggplot(aes(x=log2(site_fc)))+
  geom_histogram(bins = 50)+
  facet_wrap(~codon_display, scales = "free_y")+
  geom_vline(xintercept = -2, linetype = 'dashed', colour = 'red')+
  geom_vline(xintercept = 2, linetype = 'dashed', colour = 'red')
save_plot(paste0(figurePrefix, "_selected_codon_stalling_sitefc_l2fc_hist.pdf"), selected_to_winzorize, base_width = 12, base_height = 12)

pdf(paste0(figurePrefix, "_selected_codon_stalling_sitefc_l2fc_pseudotime.pdf"), height = 8, width = 10)
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
                                                mAG = colorRamp2(seq(0, max(RPFmeta$log_green), length.out = 8), brewer.pal(8, 'Greens')),
                                                pseudotime = colorRamp2(seq(min(RPFmeta$pseudotime), max(RPFmeta$pseudotime), length.out = 11), rev(brewer.pal(11, 'Spectral')))),
                                     annotation_name_side = "left")
Heatmap(log2(as.matrix(selectedHM)),
        right_annotation = site_labels,
        top_annotation = cell_annotation,
        col = colorRamp2(seq(-pausing_limit, pausing_limit, length = 11), rev(brewer.pal(11, "RdBu"))), # log fold change
        na_col = "grey0",
        name = "log2 fold change",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_split = factor(rep(selected_codons, each = length(levels(site_change$site))), levels = selected_codons, ordered = TRUE),
        row_title_rot = 0,
        border = TRUE,
        use_raster = TRUE,
        raster_device = "CairoPNG",
        raster_quality = 20)
dev.off()



## CGC CGU pausing heatmap with histone expression for Reviewer 2 fig R2.7
# make above heatmap with CGC and CGU, concatenate histone expression

selected_codons <- c("CGC_Arg", "CGU_Arg", "GAA_Glu")
selectedHM <- site_change %>%
  codons_heatmap(codons = selected_codons, feature = site_fc)
selectedHM <- t(selectedHM)
selectedHM <- selectedHM[, pull(RPFmeta %>% filter(!is.na(pseudotime)) %>% arrange(pseudotime) %>% select(CB))]

histone_genes <- read.table('histones_group-384.csv', sep = ',', header = TRUE, stringsAsFactors = FALSE) %>%
  inner_join(canonical_pc %>% select(gene_name, gene_id), by = c("Approved.symbol" = "gene_name")) %>%
  mutate(gene_out = gsub("_", "-", paste(gene_id, Approved.symbol, sep = "-")))

histoneHM <- as.data.frame(PI@assays$RNA[rownames(PI@assays$RNA) %in% pull(histone_genes %>% select(gene_out)), pull(RPFmeta2 %>% arrange(pseudotime) %>% select(CB)) ])
histoneHM <- sweep(histoneHM, 1, apply(histoneHM,1,mean), FUN = "-")
histoneHM <- sweep(histoneHM, 1, apply(histoneHM,1,sd), FUN = "/")


pdf(paste0(figurePrefix, "_selected_codon_stalling_sitefc_l2fc_pseudotime_withHistones.pdf"), height = 5, width = 10)
site_labels <- HeatmapAnnotation(site = anno_mark(at = grep("^CGU_Arg", rownames(selectedHM)),
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
                                                mAG = colorRamp2(seq(0, max(RPFmeta$log_green), length.out = 8), brewer.pal(8, 'Greens')),
                                                pseudotime = colorRamp2(seq(min(RPFmeta$pseudotime), max(RPFmeta$pseudotime), length.out = 11), rev(brewer.pal(11, 'Spectral')))),
                                     annotation_name_side = "left")
pausing_HM <- Heatmap(log2(as.matrix(selectedHM)),
                      right_annotation = site_labels,
                      top_annotation = cell_annotation,
                      col = colorRamp2(seq(-pausing_limit, pausing_limit, length = 11), rev(brewer.pal(11, "RdBu"))), # log fold change
                      na_col = "grey0",
                      name = "log2 fold change",
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      cluster_row_slices = FALSE,
                      cluster_column_slices = FALSE,
                      show_column_names = FALSE,
                      show_row_names = FALSE,
                      row_split = factor(rep(selected_codons, each = length(levels(site_change$site))), levels = selected_codons, ordered = TRUE),
                      row_title_rot = 0,
                      border = TRUE,
                      use_raster = TRUE,
                      raster_device = "CairoPNG",
                      raster_quality = 15)
gene_annotation <- rowAnnotation(foo = anno_mark(at = c(grep("-H3C2$", rownames(histoneHM)),
                                                        grep("-H2AZ1$", rownames(histoneHM)),
                                                        grep("-CENPA$", rownames(histoneHM)),
                                                        grep("-H1-4$", rownames(histoneHM))),
                                                 labels = c("H3C2",
                                                            "H2AZ1",
                                                            "CENPA",
                                                            "H1-4") ))
histone_HM <- Heatmap(as.matrix(histoneHM),
                      right_annotation = gene_annotation,
                      na_col = "grey0",
                      name = "z-score expression",
                      cluster_rows = TRUE,
                      show_row_dend = FALSE,
                      clustering_distance_rows = "pearson",
                      clustering_method_rows = "ward.D2",
                      cluster_columns = FALSE,
                      cluster_row_slices = FALSE,
                      cluster_column_slices = FALSE,
                      show_column_names = FALSE,
                      show_row_names = FALSE,
                      row_title_rot = 0,
                      border = TRUE,
                      use_raster = TRUE,
                      raster_device = "CairoPNG",
                      raster_quality = 15,
                      height = unit(1.5, "in"))
heat_list <- pausing_HM %v% histone_HM
draw(heat_list)
dev.off()


## for the text, what is the fraction of RPF reads that align to histone genes?
RPFlong <- as.data.frame(t(RPFcounts))
RPFlong$CB <- rownames(RPFlong)
RPFlong <- RPFlong %>%
  pivot_longer(cols = starts_with("ENSG"), names_to = "gene_id", values_to = "counts") 

histone_fraction <- RPFlong %>%
  mutate(gene_group = case_when(gene_id %in% histone_genes$gene_out ~ "histone",
                                TRUE ~ "rest")) %>%
  group_by(CB, gene_group) %>%
  summarize(cell_counts = sum(counts)) %>%
  ungroup() %>%
  group_by(CB) %>% 
  mutate(cell_frac = cell_counts/sum(cell_counts)) %>%
  ungroup() %>%
  inner_join(RPFmeta %>% select(CB, clusters), by = "CB")

plt_histone_frac <- histone_fraction %>%
  filter(gene_group == "histone") %>%
  ggplot(aes(x=clusters, y=cell_frac))+
  geom_boxplot()
save_plot(paste0(figurePrefix, "_cluster_histone_frac.pdf"), plt_histone_frac, base_width = 3, base_height = 5)

histone_fraction %>%
  filter(gene_group == "histone") %>%
  group_by(clusters) %>%
  summarize(mean_histone = mean(cell_frac),
            sd_histone = sd(cell_frac),
            n_cells = n(),
            se_histone = sd(cell_frac)/sqrt(n())) %>%
  ungroup()

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
  stat_summary(fun.data = function(x){ return( c(y=0.1, label = length(x)) ) }, geom = "text", fun = max )+
  scale_fill_manual(values = brewer.pal(3, "Greys")[c(2,3)])+
  scale_colour_manual(values = brewer.pal(3, "Greys")[c(2,3)])+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
save_plot(paste0(figurePrefix, "_histone_codon_enrichment.pdf"), plt_codon_enrichment, base_width = 6, base_height = 5)


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
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu")), na.value = "grey40", limits = c(-pausing_limit,pausing_limit), oob=squish)+
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
  inner_join(RPFmeta %>% select(CB, log_green), by = "CB") %>%
  filter(clusters == 6) %>%
  ggplot(aes(x=log_green, y=log2(site_fc)))+
  geom_point()
save_plot(paste0(figurePrefix, "_highlight_aua.pdf"), plt_aua_mit, base_width = 5, base_height = 5)



plt_umap_highlight_batch <- site_change %>%
  filter(codon_display %in% c("GAA_Glu")) %>%
  filter(site %in% c("a")) %>%
  inner_join(RPFmeta %>% select(CB, UMAP_1, UMAP_2, batch), by = "CB") %>%
  mutate(arrange_value = case_when(codon_display == "CGA_Arg" ~ -site_fc,
                                   TRUE ~ site_fc)) %>%
  arrange(arrange_value) %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, colour = log2(site_fc)))+
  geom_point(size=.5)+
  coord_equal()+
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu")), na.value = "grey40", limits = c(-pausing_limit,pausing_limit), oob=squish)+
  facet_wrap(~codon_display+batch)+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank())
save_plot(paste0(figurePrefix, "_highlight_codons_batch.pdf"), plt_umap_highlight_batch, base_width = 10, base_height = 8)



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

allHMa_hc <- dendsort(hclust(as.dist(1-cosine_similarity(as.matrix(log2(allHMa)))), method = "ward.D2"), type = "min")

all_codons <- all_codons[allHMa_hc$order]

allHM <- site_change %>%
  codons_heatmap(codons = all_codons, feature = site_fc)
allHM <- t(allHM)

allHM <- allHM[, pull(RPFmeta %>% filter(!is.na(pseudotime)) %>% arrange(pseudotime) %>% select(CB))]

all_codons_sig <- all_codons
all_codons_sig[all_codons_sig %in% sig_codons] <- paste0("*", all_codons_sig[all_codons_sig %in% sig_codons]) 

pdf(paste0(figurePrefix, "_all_codon_stalling_sitefc_l2fc_pseudotime.pdf"), height = 10, width = 9)
site_labels <- HeatmapAnnotation(site = anno_mark(at = grep("^GAA_Glu", rownames(allHM)),
                                                  labels = c("-2", "-1", "E", "P", "A", "+1", "+2")), which = "row")
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
                                                mAG = colorRamp2(seq(0, max(RPFmeta$log_green), length.out = 8), brewer.pal(8, 'Greens')),
                                                pseudotime = colorRamp2(seq(min(RPFmeta$pseudotime), max(RPFmeta$pseudotime), length.out = 11), rev(brewer.pal(11, 'Spectral')))),
                                     annotation_name_side = "left")
Heatmap(log2(as.matrix(allHM)),
        right_annotation = site_labels,
        top_annotation = cell_annotation,
        col = colorRamp2(seq(-pausing_limit, pausing_limit, length = 11), rev(brewer.pal(11, "RdBu"))), # log fold change
        na_col = "grey0",
        name = "log2 fold change",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_split = factor(rep(all_codons_sig, each = length(levels(site_change$site))), levels = all_codons_sig, ordered = TRUE),
        row_title_rot = 0,
        border = FALSE,
        use_raster = TRUE,
        raster_device = "CairoPNG",
        raster_quality = 20)
dev.off()


## Expression heatmap

DefaultAssay(PI) <- "RNA"
markers <- FindAllMarkers(PI, min.pct = 0.25, logfc.threshold = 0.25)

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

row_dend <- dendsort(hclust(as.dist(1-cosine_similarity(as.matrix(markerHM))), method = "ward.D2"), type = "min")
log2_fc_RPFpC = log2(RPFmeta2[colnames(markerHM), "cds"]/mean(RPFmeta2[colnames(markerHM), "cds"]))
pdf(paste0(figurePrefix, "_cc_marker_hm_pstwithfacs.pdf"), height = 8, width = 9)
cell_annotation <- HeatmapAnnotation(population = as.character(RPFmeta[colnames(markerHM), "population"]),
                                     clusters = as.character(RPFmeta[colnames(markerHM), "clusters"]),
                                     pseudotime = RPFmeta[colnames(markerHM), "pseudotime"],
                                     mKO2 = RPFmeta[colnames(markerHM), "log_red"],
                                     mAG = RPFmeta[colnames(markerHM), "log_green"],
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
                                                mAG = colorRamp2(seq(0, max(RPFmeta$log_green), length.out = 8), brewer.pal(8, 'Greens')),
                                                pseudotime = colorRamp2(seq(min(RPFmeta$pseudotime), max(RPFmeta$pseudotime), length.out = 11), rev(brewer.pal(11, 'Spectral'))),
                                                log2_fc_RPFpC = colorRamp2(seq(min(log2_fc_RPFpC), max(log2_fc_RPFpC), length.out = 11), rev(brewer.pal(11, 'PuOr')))),
                                     annotation_name_side = "right")
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
                                                     "SGO1",
                                                     "FBXO5",
                                                     "CDC25B",
                                                     "CDC25C" 
                                                     ) ) )
Heatmap(as.matrix(markerHM),
        top_annotation = cell_annotation,
        right_annotation = cc_annotation,
        na_col = "grey0",
        name = "z-score expression",
        cluster_rows = row_dend,
        show_row_dend = FALSE,
        clustering_distance_rows = function(x) as.dist(1-cosine_similarity(x)),
        clustering_method_rows = "ward.D2",
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        column_title_rot = 90,
        border = TRUE,
        use_raster = TRUE,
        raster_device = "CairoPNG",
        raster_quality = 15)
dev.off()




## cell_removal

plt_total_pst <- ggplot(RPFmeta, aes(x=pseudotime, y=log10(cds), colour = seurat_clusters))+
  geom_point(size=1)
save_plot(paste0(figurePrefix, "_total_pst.pdf"), plt_total_pst, base_width = 7, base_height = 5)

plt_total_pst_f <- ggplot(RPFmeta, aes(x=pseudotime, y=log10(cds), colour = seurat_clusters))+
  geom_point(size=1, data = RPFmeta %>% select(cds, pseudotime), colour = "grey70")+
  geom_point(size=1)+
  facet_wrap(~seurat_clusters+batch)
save_plot(paste0(figurePrefix, "_total_pst.pdf"), plt_total_pst_f, base_width = 15, base_height = 15)

plt_cds_cluster <- ggplot(RPFmeta, aes(x=log10(cds), fill = seurat_clusters))+
  geom_histogram(data = RPFmeta %>% select(cds), fill= "grey70", bins=50)+
  geom_histogram(bins = 50)+
  facet_wrap(~seurat_clusters, scales = "free_y")
save_plot(paste0(figurePrefix, "_total_hist_cluster.pdf"), plt_cds_cluster, base_width = 8, base_height = 8)


## example hm

asite_pausing <- site_change %>%
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

pdf(paste0(figurePrefix, "_pausing_hm_", display_gene, ".pdf"), height = 3.5, width = 9)
aa_show <- HeatmapAnnotation(aa = aa_highlight,
                             gaa_label = anno_mark(at = c(6,91),
                                                   labels = c("E6", "E91")),
                             col = list(aa = c("other" = "white",
                                               "GAA" = "red",
                                               "AUA" = "blue",
                                               "CGA" = "green")),
                             annotation_name_side = "left")
exprn <- as.numeric(PI@assays$RNA[ pull(annot %>% 
                                        filter(gene_name == display_gene) %>% 
                                        mutate(gene_out = paste(gene_id, gene_name, sep = '-')) %>% 
                                        select(gene_out) %>% 
                                        distinct()), rownames(EXG_hm) ] )
cell_annotation <- HeatmapAnnotation(exprn = exprn,
                                     GAA = log2(asite_pausing[rownames(EXG_hm), "GAA_Glu"]),
                                     clusters = as.character(RPFmeta[rownames(EXG_hm), "clusters"]),
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
                                                mAG = colorRamp2(seq(0, max(RPFmeta$log_green), length.out = 8), brewer.pal(8, 'Greens')),
                                                GAA = colorRamp2(seq(-pausing_limit, pausing_limit, length = 11), rev(brewer.pal(11, "RdBu"))),
                                                exprn = colorRamp2(seq(min(exprn), max(exprn), length.out = 6), (brewer.pal(6, "Blues"))),
                                                pseudotime = colorRamp2(seq(min(RPFmeta$pseudotime), max(RPFmeta$pseudotime), length.out = 11), rev(brewer.pal(11, 'Spectral')))),
                                     which = "row",
                                     annotation_name_side = "top")
length_annotation <- HeatmapAnnotation(pos = anno_mark(at = c(seq(25,ncol(EXG_hm), by = 50)),
                                                       labels = as.character(c(seq(25,ncol(EXG_hm), by = 50))), 
                                                       which = "column", side = "bottom", labels_rot = 0, link_height = unit(3, "mm")))
Heatmap(as.matrix(EXG_hm),
        top_annotation = aa_show,
        right_annotation = cell_annotation,
        bottom_annotation = length_annotation,
        col = colorRamp2(seq(0.025,.20, length = 9), brewer.pal(9, "Greys")),
        na_col = "red",
        name = "fraction_of_reads p cell",
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        clustering_distance_rows = "pearson", 
        clustering_method_rows = "ward.D2",
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        row_title = NULL,
        column_title_rot = 0,
        border = TRUE,
        use_raster = TRUE,
        raster_device = "CairoPNG",
        raster_quality = 20)
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


xlims <- log10(c(min(pull(cluster_change %>% filter(gene_id %in% shared_genes) %>% filter(codon %in% c("GAA", "GAG", "ATA", "CGA")) %>% select(cluster_gene_mean))),
                 max(pull(cluster_change %>% filter(gene_id %in% shared_genes) %>% filter(codon %in% c("GAA", "GAG", "ATA", "CGA")) %>% select(cluster_gene_mean)))))
ylims <- c(min(pull(cluster_change %>% filter(gene_id %in% shared_genes) %>% filter(codon %in% c("GAA", "GAG", "ATA", "CGA")) %>% mutate(l2fc_frac_change = log2(cluster_gene_frac/codon_frac)) %>% select(l2fc_frac_change))),
           max(pull(cluster_change %>% filter(gene_id %in% shared_genes) %>% filter(codon %in% c("GAA", "GAG", "ATA", "CGA")) %>% mutate(l2fc_frac_change = log2(cluster_gene_frac/codon_frac)) %>% select(l2fc_frac_change))))

plt_cluster_change_cag_facet <- cluster_change %>%
  filter(codon == "CAG") %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name == "MYL6" ~ "MYL6",
                                  TRUE ~ "")) %>%
  ggplot(aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = clusters, label = gene_display ))+
  geom_point(size=0.5)+
  coord_cartesian(xlim = xlims, ylim = ylims)+ 
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~clusters, nrow=1)+
  geom_text_repel(segment.color = "black", box.padding = 0.4, min.segment.length = 0.0)+
  xlab("log10 mean cluster expression")+ylab("log2 f.c. CAG A-site FOC")+
  theme(strip.background = element_blank(),
        legend.position = "none")
plt_cluster_change_gag_facet <- cluster_change %>%
  filter(codon == "GAG") %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name == "MYL6" ~ "MYL6",
                                  TRUE ~ "")) %>%
  ggplot(aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = clusters, label = gene_display ))+
  geom_point(size=0.5)+
  coord_cartesian(xlim = xlims, ylim = ylims)+ 
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~clusters, nrow=1)+
  geom_text_repel(segment.color = "black", box.padding = 0.4, min.segment.length = 0.0)+
  xlab("log10 mean cluster expression")+ylab("log2 f.c. GAG A-site FOC")+
  theme(strip.background = element_blank(),
        legend.position = "none")
plt_cluster_change_gaa_facet <- cluster_change %>%
  filter(codon == "GAA") %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name == "MYL6" ~ "MYL6",
                                  TRUE ~ "")) %>%
  ggplot(aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = clusters, label = gene_display ))+
  geom_point(size=0.5)+
  coord_cartesian(xlim = xlims, ylim = ylims)+ 
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~clusters, nrow=1)+
  geom_text_repel(segment.color = "black", box.padding = 0.4, min.segment.length = 0.0)+
  xlab("log10 mean cluster expression")+ylab("log2 f.c. GAA A-site FOC")+
  theme(strip.background = element_blank(),
        legend.position = "none")
plt_cluster_change_aua_facet <- cluster_change %>%
  filter(codon == "ATA") %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name == "MYL6" ~ "MYL6",
                                  TRUE ~ "")) %>%
  ggplot(aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = clusters, label = gene_display ))+
  geom_point(size=0.5)+
  coord_cartesian(xlim = xlims, ylim = ylims)+ 
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~clusters, nrow=1)+
  geom_text_repel(segment.color = "black", box.padding = 0.4, min.segment.length = 0.0)+
  xlab("log10 mean cluster expression")+ylab("log2 f.c. ATA A-site FOC")+
  theme(strip.background = element_blank(),
        legend.position = "none")
plt_cluster_change_cga_facet <- cluster_change %>%
  filter(codon == "CGA") %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(gene_display = case_when(gene_name == "MYL6" ~ "MYL6",
                                  TRUE ~ "")) %>%
  ggplot(aes(x=log10(cluster_gene_mean), y=log2(cluster_gene_frac/codon_frac), colour = clusters, label = gene_display ))+
  geom_point(size=0.5)+
  coord_cartesian(xlim = xlims, ylim = ylims)+ 
  geom_hline(yintercept = 0, colour = 'grey50', linetype = "dashed")+
  facet_wrap(~clusters, nrow=1)+
  geom_text_repel(segment.color = "black", box.padding = 0.4, min.segment.length = 0.0)+
  xlab("log10 mean cluster expression")+ylab("log2 f.c. CGA A-site FOC")+
  theme(strip.background = element_blank(),
        legend.position = "none")
save_plot(paste0(figurePrefix, "_plt_cluster_change_rest.pdf"), plot_grid(plt_cluster_change_gaa_facet, plt_cluster_change_gag_facet, plt_cluster_change_aua_facet, plt_cluster_change_cga_facet, ncol = 1), base_width = 6, base_height = 2*4)



### MYL6 example numbers for R2 question 3
# For eample, in the cluster representing mitotic cells (purple, clusters == 6), there are 
cluster_change %>%
  filter(gene_name == "MYL6",
         clusters == 6) %>%
  summarize(n_cluster = sum(n))
# in-frame RPFs measured for MYL6, 
cluster_change %>%
  filter(gene_name == "MYL6",
         clusters == 6,
         codon == "GAA") %>%
  summarize(n_cluster_GAA = sum(n))
# of which contain a GAA in their A-site. 
# Similarly, in the entire dataset, 
gene_codon_background %>%
  filter(gene_name == "MYL6",
         codon == "GAA") %>%
  summarize(n_GAA = sum(n))
# of the
gene_codon_background %>%
  filter(gene_name == "MYL6") %>%
  summarize(n_total = sum(n))



plt_cluster_box <- cluster_change %>%
  #filter(codon %in% c("GAA")) %>%
  filter(codon %in% c("GAA","GAG", "ATA", "CGA")) %>%
  mutate(display_order = case_when(codon == "GAA" ~ 1,
                                   codon == "GAG" ~ 2,
                                   codon == "ATA" ~ 3,
                                   codon == "CGA" ~ 4)) %>%
  mutate(codon = fct_reorder(codon, display_order)) %>%
  filter(gene_id %in% shared_genes) %>%
  #filter(codon %in% c("GAA", "ATA", "CGA", "GAG")) %>%
  mutate(gene_display = case_when(gene_name == "MYL6" ~ "MYL6",
                                  TRUE ~ "")) %>%
  ggplot(aes(x=clusters, y=log2(cluster_gene_frac/codon_frac), colour = clusters ))+
  geom_boxplot()+
  geom_hline(yintercept = 0, colour = 'grey50', linetype = 'dashed')+
  facet_wrap(~codon, nrow=1)+
  theme(legend.position = "none",
        strip.background = element_blank())
  #geom_text_repel(segment.color = "black", box.padding = 0.4)
save_plot(paste0(figurePrefix, "_plt_cluster_box.pdf"), plt_cluster_box, base_width = 7, base_height =3)


cluster_diff_stats <- do.call('bind_rows', lapply(levels(cluster_change$clusters), function(cl){
                                          cdstats <- test_clusters %>%
                                            mutate(group = case_when(clusters == cl ~ cl,
                                                                     clusters != cl ~ "rest")) %>%
                                            group_by(codon) %>%
                                            summarize(p_val = wilcox.test(lf2_change ~ group)$p.value, .groups = "drop_last"
                                                      ) %>%
                                            mutate(cluster = cl) 
                                  } ) ) %>% 
  mutate(p_val_adj = p.adjust(p_val, method = "bonferroni"))
                                            

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

getMode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

PIp <- PI
PIp[['GAA_apausing']] <- CreateAssayObject(data = as.matrix(GAA_a_genes))
gaa_a_pausing_markers <- FindAllMarkers(PIp, assay = "GAA_apausing", logfc.threshold = 0, only.pos = FALSE, pseudocount.use = 0, min.cells.feature = 10, min.cells.group =0) %>%
  inner_join(RPFmeta %>%
               group_by(seurat_clusters) %>%
               summarize(cluster_redef = getMode(clusters)) %>%
               ungroup(),
             by = c("cluster" = "seurat_clusters")) %>%
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
  #inner_join(RPFmeta %>% select(CB, clusters), by="CB") %>%
  mutate(grouping_cluster = case_when(clusters %in% c(2,7) ~ "2_7",
                                      TRUE ~ "rest")) %>%
  group_by(grouping_cluster) %>%
  summarize(mean_site_fraction = mean(100*cell_site_fraction),
            sd_site_fraction = sd(100*cell_site_fraction),
            n = n(),
            sem_site_fraction = sd(100*cell_site_fraction)/sqrt(n())) %>%
  ungroup()

ugc_stats_test <- site_change %>%
  filter(codon_display == "UGC_Cys") %>%
  #inner_join(RPFmeta %>% select(CB, clusters), by="CB") %>%
  mutate(grouping_cluster = case_when(clusters %in% c(2,7) ~ "2_7",
                                      TRUE ~ "rest")) %>%
  select(cell_site_fraction, grouping_cluster)
wilcox.test(cell_site_fraction ~ grouping_cluster, data = ugc_stats_test)

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
            sd_ugc = sd(100*cell_codon_fraction),
            sem_ugc = sd(100*cell_codon_fraction)/sqrt(n())) %>%
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
            sd_gaa = sd(100*cell_codon_fraction),
            sem_gaa = sd(100*cell_codon_fraction)/sqrt(n())) %>%
  ungroup()

gaa_stats_test <- codon_reads %>%
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
  mutate(grouping_cluster = case_when(clusters == 6 ~ "6",
                                      TRUE ~ "rest")) %>%
  wilcox.test(cell_codon_fraction ~ grouping_cluster, data = .)


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
            sd_gaa = sd(100*cell_codon_fraction),
            sem_gaa = sd(100*cell_codon_fraction)/sqrt(n())) %>%
  ungroup()

## GAA a-stie average increase
cl6_increase <- cluster_change %>%
  filter(codon == "GAA") %>%
  filter(gene_id %in% shared_genes) %>%
  mutate(grouping_cluster = case_when(clusters == 6 ~ 6,
                                      TRUE ~ -1)) %>%
  group_by(grouping_cluster) %>%
  summarize(mean_gaa = mean(cluster_gene_frac/codon_frac),
            sd_gaa = sd(cluster_gene_frac/codon_frac),
            sem_gaa = sd(cluster_gene_frac/codon_frac)/sqrt(n())) %>%
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

# to filter fastq
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

MGI <- barcodes
MGI[ MGI$well %in% columnFeatures(1:15, LETTERS[1:16]), "sort_population"] = "Interphase"
MGI[ MGI$well %in% columnFeatures(16:21, LETTERS[1:16]), "sort_population"] = "Mit"
MGI[ MGI$well %in% columnFeatures(22:24, LETTERS[1:14]), "sort_population"] = "G0"
MGI[ MGI$well %in% columnFeatures(22:24, LETTERS[15:16]), "sort_population"] = "NTC"
write.csv(MGI, file = paste0(figurePrefix, "_RPE_MGI_layout.csv"), row.names = FALSE, quote = FALSE)
