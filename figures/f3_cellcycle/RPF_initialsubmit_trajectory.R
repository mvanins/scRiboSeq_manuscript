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
#library(Matrix)
#library(igraph)

source('/hpc/hub_oudenaarden/mvanins/local/avopipelines/RPF/bin/scRibo/plotFunctions.R')
source('TSPtrajectory.R')
# source('Swanderlust.R')
figurePrefix <- "RPF_IS_trajectory_bigger"

# barcodes
barcodes <- read.table('/hpc/hub_oudenaarden/mvanins/local/lib/barcodes/miRv6_barcodemap.tsv',
                       sep="\t",
                       header = FALSE,
                       stringsAsFactors = FALSE,
                       row.names = NULL,
                       col.names = c("CB","number","well"))
rownames(barcodes) <- barcodes$CB

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
                         counts = c("../../../data/Fucci/counts/RPFv4-Fucci-G011",
                                    "../../../data/Fucci/counts/RPFv4-Fucci-G012",
                                    "../../../data/Fucci/counts/RPFv4-Fucci-Mit01",
                                    "../../../data/Fucci/counts/RPFv4-Fucci-Mit02",
                                    "../../../data/Fucci/counts/RPFv4-Fucci-Mit03",
                                    "../../../data/Fucci/counts/RPFv4-Fucci-Mit04",
                                    "../../../data/Fucci/counts/RPFv4-Fucci-Mit06",
                                    "../../../data/Fucci/counts/RPFv4-Fucci-Mit07",
                                    "../../../data/Fucci/counts/RPFv4-Fucci-Mit08",
                                    "../../../data/Fucci/counts/RPFv4-CDKN1AKO01",
                                    "../../../data/Fucci/counts/RPFv4-CDKN1AKO02",
                                    "../../../data/Fucci/counts/RPFv4-CDKN1AKO03",
                                    "../../../data/Fucci/counts/RPFv4-CDKN1AKO04",
                                    "../../../data/Fucci/counts/RPFv4-SR201",
                                    "../../../data/Fucci/counts/RPFv4-SR202",
                                    "../../../data/Fucci/counts/RPFv4-SR203",
                                    "../../../data/Fucci/counts/RPFv4-SR204",
                                    "../../../data/Fucci/counts/RPFv4-SR205",
                                    "../../../data/Fucci/counts/RPFv4-SR206"),
                         facs = c("../../../data/Fucci/facs/G011.csv",
                                  "../../../data/Fucci/facs/G012.csv",
                                  "../../../data/Fucci/facs/Mit01.csv",
                                  "../../../data/Fucci/facs/Mit02.csv",
                                  "../../../data/Fucci/facs/Mit03.csv",
                                  "../../../data/Fucci/facs/Mit04.csv",
                                  "../../../data/Fucci/facs/Mit06.csv",
                                  "../../../data/Fucci/facs/Mit07.csv",
                                  "../../../data/Fucci/facs/Mit08.csv",
                                  "../../../data/Fucci/facs/KO01.csv",
                                  "../../../data/Fucci/facs/KO02.csv",
                                  "../../../data/Fucci/facs/KO03.csv",
                                  "../../../data/Fucci/facs/KO04.csv",
                                  "../../../data/Fucci/facs/SR201.csv",
                                  "../../../data/Fucci/facs/SR202.csv",
                                  "../../../data/Fucci/facs/SR203.csv",
                                  "../../../data/Fucci/facs/SR204.csv",
                                  "../../../data/Fucci/facs/SR205.csv",
                                  "../../../data/Fucci/facs/SR206.csv"),
                         stringsAsFactors = FALSE)

RPFfacs <- readFACS(RPFsamples, barcodes = barcodes) %>%
  mutate(cgreen = X..488..530.40,
         cred = X..561..585.29) %>%
  mutate(rgreen = X.488..530.40,
         rred = X.561..585.29) #%>%
  #filter(sample != "G0")

RPFcounts <- readMergeCounts(RPFsamples)

## Assemble treatment information
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

RPFpts <- RPFfacs %>%
  inner_join(RPFmeta, by="CB") %>%
  #filter(sort_population != "G0") %>%
  #filter(cred > 1.5, cgreen > 1.5) %>%
  mutate(cgreen = log10(cgreen),
         cred = log10(cred)) %>%
  mutate(rgreen = log10(rgreen),
         rred = log10(rred)) %>%
  select(CB, cred, cgreen, rred, rgreen, batch) %>%
  as.data.frame()
rownames(RPFpts) <- RPFpts$CB

RPFcells <- pull(RPFmeta %>% 
                 filter(cds_frac >= 0.85, cell_total > 4500) %>%
                 filter(!(sort_population %in% c("NTC"))) %>%
                 filter(!(batch %in% c("Sync"))) %>%
                 filter(!(sort_population %in% c("KO14", "KO13", "KO15", "KO17",
                                       "9h_post", "5h_post", "2h_post", "20h_post", "10h_post", "6h_post"))) %>%
                 #filter(batch %in% c("First")) %>%
                 #filter(!(batch %in% c("NB", "Sync"))) %>%
                 select(CB))

RPFpts <- RPFpts[RPFcells,]

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

## simple offsets
# FACSoffsets <- RPFpts %>%
#   filter(CB %in% pull(RPFmeta %>% filter(sort_population != "Mit") %>% select(CB))) %>%
#   filter(cgreen > 0) %>%
#   group_by(batch) %>%
#   summarize(mean_cred = mean(cred),
#             mean_cgreen = mean(cgreen),
#             mean_rred = mean(rred),
#             mean_rgreen = mean(rgreen)) %>%
#   ungroup() %>%
#   mutate(off_cred = mean_cred - nth(x=mean_cred, n=which(batch == "RPF_First")),
#          off_cgreen = mean_cgreen - nth(x=mean_cgreen, n=which(batch == "RPF_First")),
#          off_rred = mean_rred - nth(x=mean_rred, n=which(batch == "RPF_First")),
#          off_rgreen = mean_rgreen - nth(x=mean_rgreen, n=which(batch == "RPF_First"))) %>%
#   select(batch, starts_with("off_"))

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

## Compensated start list
# RPF_start_list <- pull(RPFpts %>%
#   filter((cred > 1.32 & cred < 1.80) & cgreen < 1.25) %>%
#   select(CB))
# plt_start <- ggplot(pts %>% mutate(type = "RPF") %>% mutate(start = case_when(CB %in% c(RPF_start_list) ~ "start", TRUE ~ "not")), aes(x=cgreen, y=cred, fill=type, colour = start))+
#   geom_point(shape=21)+
#   facet_wrap(~type)+
#   coord_equal()
# save_plot(paste0(figurePrefix, "_start.pdf"), plt_start, base_width = 10, base_height = 6)

## raw start list
RPF_start_list <- pull(pts %>%
  filter((rred > 1.32 & rred < 1.80) & rgreen < 1.25 & rgreen > 0.6) %>%
  select(CB))

plt_start <- ggplot(pts %>% mutate(type = "RPF") %>% mutate(start = case_when(CB %in% c(RPF_start_list) ~ "start", TRUE ~ "not")), aes(x=rgreen, y=rred, fill=type, colour = start))+
  geom_point(shape=21)+
  facet_wrap(~type)+
  coord_equal()
save_plot(paste0(figurePrefix, "_start.pdf"), plt_start, base_width = 10, base_height = 6)


### solve tours

cores <- 48*2
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



tours_summary <- tours %>%
  inner_join(RPFpts, by = "CB") %>%
  group_by(CB, cred, cgreen) %>%
  summarize(median_path = median(path_distance),
            mean_path = mean(path_distance),
            n_tours = n()) %>%
  ungroup() 

plt_tours_median <- ggplot(tours_summary %>% arrange(median_path), aes(x=cgreen,y=cred,colour=median_path))+
  geom_point(size=1)+
  scale_colour_distiller(palette = "Spectral")+
  coord_equal()
plt_tours_mean <- ggplot(tours_summary %>% arrange(mean_path), aes(x=cgreen,y=cred,colour=mean_path))+
  geom_point(size=1)+
  scale_colour_distiller(palette = "Spectral")+
  coord_equal()
plt_tours_n <- ggplot(tours_summary %>% arrange(n_tours), aes(x=cgreen,y=cred,colour=log10(n_tours)))+
  geom_point(size=1)+
  scale_colour_distiller(palette = "Spectral")+
  coord_equal()
plt_cdf <- ggplot(tours_summary, aes(x=median_path))+
  stat_ecdf(geom = "step")

save_plot(paste0(figurePrefix,"_all_compensated_tours_", niter, "iter_", npoints, "pts", "_1noise_summary.pdf"), plot_grid(plt_tours_mean, plt_tours_median, plt_tours_n, plt_cdf, nrow=2), base_width = 15, base_height = 15)



dynamics <- rbind(RPFtours_summary %>% mutate(type = "RPF"),
                  CS2tours_summary %>% mutate(type = "CS2")) %>%
  select(CB, type, median_path, cgreen, cred) %>%
  pivot_longer(cols = c("cgreen", "cred"), names_to = "channel", values_to = "fluorescence") %>%
  mutate(channel = gsub("^c", "", channel))


plt_dynamics <- ggplot(dynamics, aes(x=median_path, y=fluorescence, colour = channel))+
  geom_point()+
  scale_colour_manual(values = c("#4daf4a", "#e41a1c")) +
  facet_wrap(~type, ncol = 1)
#   facet_wrap(~channel)
save_plot(paste0(figurePrefix, "_fluo_dynamics.pdf"), plt_dynamics, base_width =8, base_height = 8)

plt_methcomp <- dynamics %>%
  group_by(type) %>%
  mutate(max_path = max(median_path)) %>%
  ungroup() %>%
  ggplot(aes(x=median_path/max_path, y=fluorescence, colour = channel))+
  geom_point()+
  scale_colour_manual(values = c("#4daf4a", "#e41a1c")) +
  facet_wrap(~type, ncol = 1)
save_plot(paste0(figurePrefix, "_fluo_dynamics.pdf"), plt_methcomp, base_width =8, base_height = 8)
