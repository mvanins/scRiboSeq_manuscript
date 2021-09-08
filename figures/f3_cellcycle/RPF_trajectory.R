#!/usr/bin/env Rscript
# Michael VanInsberghe 2020-06-29
# Defines trajectory through cell cycle using TSP

library(dplyr)
library(tidyr)
library(forcats)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(RColorBrewer)
theme_set(theme_cowplot())

library(TSP)
library(doParallel)
library(foreach)
library(doRNG)

library(RANN)

source('plotFunctions.R')
source('TSPtrajectory.R')
figurePrefix <- "RPF_replicate_trajectory_20210413"

data_dir <- "/hpc/hub_oudenaarden/mvanins/RPF/analysis/scRiboSeq/data"

# barcodes
barcodes <- read.table('/hpc/hub_oudenaarden/mvanins/local/lib/barcodes/miRv6_barcodemap.tsv',
                       sep="\t",
                       header = FALSE,
                       stringsAsFactors = FALSE,
                       row.names = NULL,
                       col.names = c("CB","number","well"))
rownames(barcodes) <- barcodes$CB


RPFmeta <- readRDS(file.path(data_dir, "Fucci/compiled/meta.rds"))

RPFpts <- RPFmeta %>%
  mutate(cgreen = X..488..530.40,
         cred = X..561..585.29) %>%
  mutate(rgreen = X.488..530.40,
         rred = X.561..585.29) %>%
  mutate(cgreen = log10(cgreen),
         cred = log10(cred)) %>%
  mutate(rgreen = log10(rgreen),
         rred = log10(rred)) %>%
  select(CB, cred, cgreen, rred, rgreen, batch) %>%
  as.data.frame()
rownames(RPFpts) <- RPFpts$CB


# how do the two FACS experiments overlap?
# FACS
# Green: 488 530/40 X..488..530.40
# Red: 561 585/29 X..561..585.29

pts <- rbind(RPFpts %>% mutate(type = "RPF"))
rownames(pts) <- pts$CB

plt_facs_comparison_comp <- ggplot(pts, aes(y=cred, x=cgreen, colour = batch)) +
  geom_point(data = pts %>% select(cred, cgreen), colour = "grey80")+
  geom_point(size=1)+
  facet_wrap(~batch, nrow=1)+
  coord_equal()
plt_facs_comparison_raw <- ggplot(pts, aes(y=rred, x=rgreen, colour = batch)) +
  geom_point(data=pts %>% select(rred, rgreen), colour = "grey80")+
  geom_point(size=1)+
  facet_wrap(~batch, nrow=1)+
  coord_equal()

save_plot(paste0(figurePrefix, "_facs_comparison.pdf"), plot_grid(plt_facs_comparison_raw, plt_facs_comparison_comp,ncol=1), base_width = 20, base_height = 8)

nn_dist <- function(tbl, x, y){
  if(dplyr::is_grouped_df(tbl)){
    indices <- group_indices(tbl)
    tbl <- do(tbl, nn_dist(., x=x, y=y))
    tbl <- tbl[order(order(indices)), , drop = FALSE]
    return(tbl)
  }

  tbl$nn_dist <- nn2(tbl %>% select(.data[[x]], .data[[y]]), k=2)$nn.dists[,2]

  return(tbl)
}

nn_downsample <- function(tbl, x=x, y=y, rounds = 10, nn_final = 100, final = 50, noise_fraction = 1.0, seed = NULL){
  if(dplyr::is_grouped_df(tbl)){
    tblds <- do(tbl, nn_downsample(., x=x, y=y, rounds = rounds, nn_final = nn_final, final = final, noise_fraction = noise_fraction, seed = seed))
    return(tblds)
  }

  downsample_fraction = (nn_final/nrow(tbl))^(1/rounds)
  if(!is.null(seed)){
    set.seed(seed)
  }

  tblds <- tbl
  for(i in 1:rounds){
    nndists <- nn2(tblds %>% select(.data[[x]], .data[[y]]), k=2)$nn.dists[,2]
    tblds$nn_dist <- nndists + runif(nrow(tblds))*(noise_fraction*mean(nndists))
    tblds <- tblds %>%
      slice_sample(n=round(downsample_fraction * nrow(tblds)), weight_by = nn_dist)
  }

  return(tblds %>% slice_sample(n=final) %>% select(-nn_dist))
}


# downsample, single subtractive offset
reference_dataset = "RPF_First"
offsets <- data.frame()
set.seed(112)
for(i in 1:500){
  offset_samp <- RPFpts %>%
    filter(CB %in% pull(RPFmeta %>% filter(sort_population != "Mit") %>% select(CB))) %>%
    filter(cgreen > 0.2) %>%
    group_by(batch) %>%
    nn_downsample(x="rgreen", y="rred", rounds = 5, nn_final = 100, final = 30, noise_fraction = 1.5) %>%
    summarize(mean_cred = mean(cred),
              mean_cgreen = mean(cgreen),
              mean_rred = mean(rred),
              mean_rgreen = mean(rgreen), .groups = "drop_last") %>%
    ungroup() %>%
    mutate(off_cred = mean_cred - nth(x=mean_cred, n=which(batch == reference_dataset)),
           off_cgreen = mean_cgreen - nth(x=mean_cgreen, n=which(batch == reference_dataset)),
           off_rred = mean_rred - nth(x=mean_rred, n=which(batch == reference_dataset)),
           off_rgreen = mean_rgreen - nth(x=mean_rgreen, n=which(batch == reference_dataset))) %>%
    select(batch, starts_with("off_")) %>%
    pivot_longer(!batch, names_to = "channel", values_to = "offset") %>%
    mutate(index = i)

    offsets <- rbind(offsets, offset_samp)
}

FACSoffsets <- offsets %>%
  group_by(batch, channel) %>%
  summarize(mean_offset = mean(offset)) %>%
  ungroup() %>%
  pivot_wider(names_from = channel, values_from = mean_offset)

RPFptsC <- RPFpts %>%
  filter(CB %in% pull(RPFmeta %>% filter(sort_population != "Mit") %>% select(CB))) %>%
  inner_join(FACSoffsets, by = "batch") %>%
  mutate(cor_cred = cred - off_cred,
         cor_cgreen = cgreen - off_cgreen,
         cor_rred = rred - off_rred,
         cor_rgreen = rgreen - off_rgreen) %>%
  select(-starts_with("off_"))

plt_facs_comparison_comp <- ggplot(RPFptsC, aes(y=cor_cred, x=cor_cgreen, colour = batch)) +
  geom_point(data = RPFptsC %>% select(cor_cred, cor_cgreen), colour = "grey80")+
  geom_point(size=1)+
  facet_wrap(~batch, nrow=1)+
  coord_equal()
plt_facs_comparison_raw <- ggplot(RPFptsC, aes(y=cor_rred, x=cor_rgreen, colour = batch)) +
  geom_point(data=RPFptsC %>% select(cor_rred, cor_rgreen), colour = "grey80")+
  geom_point(size=1)+
  facet_wrap(~batch, nrow=1)+
  coord_equal()

save_plot(paste0(figurePrefix, "_facs_comparison_shifted.pdf"), plot_grid(plt_facs_comparison_raw, plt_facs_comparison_comp,ncol=1), base_width = 20, base_height = 8)


pts <- RPFpts %>%
  inner_join(FACSoffsets, by = "batch") %>%
  mutate(cor_cred = cred - off_cred,
         cor_cgreen = cgreen - off_cgreen,
         cor_rred = rred - off_rred,
         cor_rgreen = rgreen - off_rgreen) %>%
  select(CB, batch, starts_with("cor_")) %>%
  rename_at(.vars = vars(starts_with("cor_")), .funs = list(~sub("cor_", "", .)))



# fit separate trajectories for each

## raw start list
RPF_start_list <- pull(pts %>%
  filter((rred > 1.39 & rred < 1.80) & rgreen < 1.25 & rgreen > 0.5) %>%
  select(CB))

plt_start <- ggplot(pts %>% mutate(type = "RPF") %>% mutate(start = case_when(CB %in% c(RPF_start_list) ~ "start", TRUE ~ "not")), aes(x=rgreen, y=rred, fill=type, colour = start))+
  geom_point(shape=21)+
  geom_vline(xintercept = c(.5, 1.25), colour = 'grey50', linetype = 'dashed')+
  geom_hline(yintercept = c(1.39, 1.80), colour = 'grey50', linetype = 'dashed')+
  facet_wrap(~type)+
  coord_equal()
save_plot(paste0(figurePrefix, "_start.pdf"), plt_start, base_width = 10, base_height = 6)


### solve tours

cores <- 24*2
registerDoParallel(cores=cores)

niter <- 500000
npoints <- 30

# raw tours
# all
set.seed(112)
tours <- foreach(a=idivix(niter, chunks = cores), .combine = 'rbind') %dorng% {
  do.call('rbind', lapply(seq(a$i, length.out = a$m), function(i) {
    solveTour(RPFpts,
              iteration = i,
              x="rgreen", y="rred", 
              start_list = RPF_start_list,
              downsample_rounds = 5, 
              downsample_nn = 10*npoints,
              downsample_final = npoints,
              downsample_noise_fraction = 1)
    }
  ))
}

saveRDS(tours, paste0(figurePrefix,"_all_raw_tours_", niter, "iter_", npoints, "pts", "_1noise.RDS"))
