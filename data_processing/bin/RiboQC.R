#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(feather))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggseqlogo))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(forcats))
theme_set(theme_cowplot())
source("/hpc/hub_oudenaarden/mvanins/local/avopipelines/RPF/bin/scRibo/plotFunctions.R")

option_list = list(
                   make_option(c("-r", "--reads"), type = "character", default = NULL, help = "feather file containing read information dataframe", metavar="character"),
                   make_option(c("-a", "--annotations"), type = "character", default = NULL, help = "feather file containing annotation information", metavar="character"),
                   make_option(c("-p", "--prefix"), type = "character", default = NULL, help = "prefix for output filenames", metavar = "character"),
                   make_option(c("-w", "--whitelist"), type = "character", default = NULL, help = "whitelist for barcodes to create plate layout plots. Assumed to be in increasing well order A1-A24, B1-B24, etc.", metavar = "character") )

opt_parser = OptionParser(option_list = option_list)
args = parse_args(opt_parser)


## for testing
# args <- data.frame(reads = "/hpc/hub_oudenaarden/mvanins/RPF/data/CS2_Fucci/RPF_update/scRibo/CDKN1AKO05-CS2_Aligned.toTranscriptome.out.dedup.sorted.csv.gz",
#                    annotations = "/hpc/hub_oudenaarden/mvanins/local/lib/reference/gencode/hsa_33/RPF_reference/annotations/gencode.v33.chr_patch_hapl_scaff.annotation.annotations.feather",
#                    whitelist = "/hpc/hub_oudenaarden/mvanins/local/avopipelines/RPF/barcodes/CS2.whitelist",
#                    stringsAsFactors = FALSE)

if(is.null(args$reads) | is.null(args$annotations)){
  print_help(opt_parser)
  stop("Please supply reads and annotations", call.= FALSE)
}

if(is.null(args$prefix)){
  figurePrefix = strsplit(tools::file_path_sans_ext(basename(args$reads)),"_")[[1]][1]
}else{
  figurePrefix = args$prefix
}

print(figurePrefix)

if(grepl("feather",args$annotations)){
  annot <- read_feather(args$annotations)
}else{
  annot <- fread(args$annotations)
}

if(grepl("feather",args$reads)){
  reads <- read_feather(args$reads)
}else{
  reads <- fread(args$reads)
}

if(!is.null(args$whitelist)){
  barcodes <- read.table(args$whitelist, header = FALSE, col.names = c("CB"), stringsAsFactors = FALSE)
  if(nrow(barcodes) == 96){
    barcodes$well <- paste0(rep(LETTERS[1:8],each=12),rep(1:12,8))
    barcodes$row <- rep(LETTERS[1:8], each=12)
    barcodes$column <- rep(1:12,8)
  }else if(nrow(barcodes) == 384){
    barcodes$well <- paste0(rep(LETTERS[1:16], each=24), rep(1:24,16))
    barcodes$row <- rep(LETTERS[1:16], each=24)
    barcodes$column <- rep(1:24,16)
  }
}

canonical_pc <- annot %>%
  filter(set=="canonical", transcript_type == "protein_coding") %>%
  filter(chr != "chrM") %>%
  filter(transcript_id != "ENST00000437139.7")

pc_reads <- canonical_pc %>%
  inner_join(reads, by = "transcript_id") %>%
  filter(length > 20)

lr_cds <- 200
lr_utr3 <- 100
lr_utr5 <- 100


pdf(paste0(figurePrefix,"_RiboQC.pdf"), height = 8, width = 8)

### bulk 


plt_bulk_metawiggle <- wiggleDF(pc_reads, group = "bulk", site = "cut5") %>%
  ggplot(aes(x=csd, y=position_total))+
  geom_bar(stat="identity")+
  facet_wrap(~name,scales="free_x")+
  xlab("distance from codon")+ylab("number of 5' cuts")+
  ggtitle(paste0(figurePrefix," bulk wiggle metaplot"))

print(plt_bulk_metawiggle)

## transcript wiggle metaplot

start_df_b <- pc_reads %>%
  mutate(csd5 = cut5-cds_start) %>%
  filter(csd5 > -50, csd5 < 100) %>%
  group_by(transcript_id,csd5) %>%
  summarise(n=n()) %>%
  group_by(transcript_id) %>%
  mutate(tot = sum(n,na.rm=TRUE)) %>%
  ungroup() %>%
  complete(csd5,transcript_id) %>%
  mutate(set="start")

stop_df_b <- pc_reads %>%
  mutate(csd5 = cut5- (cds_end+3)) %>%
  filter(csd5 > -100, csd5 < 50) %>%
  group_by(transcript_id,csd5) %>%
  summarise(n=n()) %>%
  group_by(transcript_id) %>%
  mutate(tot = sum(n,na.rm=TRUE)) %>%
  ungroup() %>%
  complete(csd5,transcript_id) %>%
  mutate(set="stop")

bulk_gene_wiggle <- bind_rows(start_df_b,stop_df_b) %>%
  group_by(transcript_id) %>%
  mutate(bothtot = sum(n,na.rm=TRUE)) %>%
  filter(bothtot>100) %>%
  ungroup() %>%
  mutate(transcript_id = fct_reorder(transcript_id, bothtot))

if(length(unique(bulk_gene_wiggle$set)) > 1){
  plt_bulk_gene_wiggle <- ggplot(bulk_gene_wiggle,aes(x=csd5,y=transcript_id,fill=(n/tot)))+
    geom_tile()+
    scale_fill_gradientn(colours=rev(brewer.pal(11,"RdYlBu")),na.value="grey40")+
    scale_y_discrete(expand=c(0,0))+
    scale_x_continuous(expand=c(0,0))+
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border=element_blank(),
          axis.line = element_blank())+
facet_wrap(~set,nrow=1,strip.position = "top",scales="free_x")+
xlab("distance from codon")+
ggtitle(paste0(figurePrefix," bulk transcript metaplot"))

  print(plt_bulk_gene_wiggle)
}
## scaled CDS
# scaledf_b <- pc_reads %>%
#   mutate( site = cut5) %>%
#   mutate( pos = case_when( site < l_utr5 ~ round( (site/l_utr5)*lr_utr5 ),
#                            ((site >= l_utr5) & (site < (l_utr5+l_cds) )) ~ round( ((site-l_utr5)/l_cds)*lr_cds ) + lr_utr5, 
#                            (site >= (l_utr5 + l_cds)) ~ round( ((site - l_utr5 - l_cds)/l_utr3)*lr_utr3 ) + lr_utr5 + lr_cds ) ) %>%
#   filter(pos > 0) %>%
#   group_by(transcript_id,pos) %>%
#   summarise(n=n()) %>%
#   ungroup() %>%
#   complete(pos = full_seq(c(1,lr_utr5 + lr_cds + lr_utr3),1),transcript_id) %>%
#   mutate(n=replace_na(n,0)) %>%
#   mutate(tot=sum(n)) %>%
#   filter(tot> 100) %>%
#   mutate(transcript_id = fct_reorder(transcript_id,tot)) 
# 
# hm_scaled_gene <- ggplot(scaledf_b,aes(x=pos,y=transcript_id,fill=(n/tot)))+
#   geom_tile()+
#   scale_fill_gradientn(colours=rev(brewer.pal(11,"RdYlBu")),na.value="grey40")+
#   scale_y_discrete(expand=c(0,0))+
#   scale_x_continuous(expand=c(0,0),
#                      breaks=c(0,(lr_utr5/2),lr_utr5,(lr_utr5+lr_cds/2),(lr_utr5+lr_cds),(lr_utr5+lr_cds+(lr_utr3/2)),(lr_utr5+lr_cds+lr_utr3)),
#                      labels=c("","5' UTR","","CDS","","3' UTR",""))+
#   theme(axis.text.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.border=element_blank())+
#   xlab("distance along scaled metagene")+
#   ggtitle(paste0(figurePrefix," CDS per gene"))
# 
# print(hm_scaled_gene)

plt_bulk_coverage <- scaleDF(pc_reads, group = "bulk", site = "cut5", utr5 = lr_utr5, cds = lr_cds, utr3 = lr_utr3) %>%
  ggplot(aes(x=pos,y=n))+
  geom_line()+
  scale_x_continuous(breaks=c(0,(lr_utr5/2),lr_utr5,(lr_utr5+lr_cds/2),(lr_utr5+lr_cds),(lr_utr5+lr_cds+(lr_utr3/2)),(lr_utr5+lr_cds+lr_utr3)),
                     labels=c("","5' UTR","","CDS","","3' UTR",""))+
xlab("position along scaled metagene")+
ylab("number of 5' cuts")+
ggtitle(paste0(figurePrefix," CDS metagene"))

print(plt_bulk_coverage)



## icon

seqdf5 <- pc_reads %>%
  mutate(motif = paste0(base5m8,base5m7,base5m6,base5m5,base5m4,base5m3,base5m2,base5m1,base5,base5p1,base5p2,base5p3,base5p4,base5p5,base5p6,base5p7,base5p8)) %>%
  filter(nchar(motif) == 17) %>%
  select(motif)

seqdf3 <- pc_reads %>%
  mutate(motif = paste0(base3m8,base3m7,base3m6,base3m5,base3m4,base3m3,base3m2,base3m1,base3,base3p1,base3p2,base3p3,base3p4,base3p5,base3p6,base3p7,base3p8)) %>%
  filter(nchar(motif) == 17) %>%
  select(motif)

plt_icon5 <- ggplot() + 
  geom_logo(seqdf5$motif)+
  scale_x_continuous(breaks=1:17,labels = -8:8)+
  xlab("distance from 5' base")

plt_icon5_p <- ggplot() + 
  geom_logo(seqdf5$motif,method='prob')+
  scale_x_continuous(breaks=1:17,labels = -8:8)+
  xlab("distance from 5' base")

plt_icon3 <- ggplot() +
  geom_logo(seqdf3$motif)+
  scale_x_continuous(breaks=1:17,labels = -8:8)+
  xlab("distance from 3' base")

plt_icon3_p <- ggplot() +
  geom_logo(seqdf3$motif,method='prob')+
  scale_x_continuous(breaks=1:17,labels = -8:8)+
  xlab("distance from 3' base")

title <- ggdraw()+draw_label(paste0(figurePrefix," cut location icon"),x=0,hjust=0)+theme(plot.margin=margin(0,0,0,7))
plt_icon_figure <- plot_grid(plt_icon5,plt_icon5_p,plt_icon3,plt_icon3_p,ncol=2)+  ggtitle(figurePrefix)
plt_icon <- plot_grid(title,plt_icon_figure,ncol=1,rel_heights=c(0.1,1))

print(plt_icon)




## periodicity-length heatmap

plt_length_periodicity <- lengthMetaHeatmap(pc_reads)
print(plt_length_periodicity)


plt_length <- pc_reads %>%
  group_by(length) %>%
  summarize(n=n()) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=length,y=n))+
  geom_bar(stat="identity")+
  ylab("number of alignments")+
  ggtitle(paste0(figurePrefix," length distribution"))

print(plt_length)


plt_region <- pc_reads %>%
  mutate(site=cut5) %>%
  mutate(region = case_when( site < l_utr5 ~ "utr5",
                            ((site >= l_utr5) & (site < (l_utr5+l_cds) )) ~ "cds", 
                            (site >= (l_utr5 + l_cds)) ~ "utr3" )) %>%
  group_by(region) %>%
  summarize(n=n()) %>%
  mutate(total = sum(n)) %>%
  ggplot(aes(x=1,y=n/total,fill=region))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("grey80","grey60","grey30"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
ylab("fraction of alignments per region")

print(plt_region)


if ("CB" %in% colnames(reads)){

  ## single-cell metaheatmaps

  ### Single-cell
  ## SC metaheatmaps

  sc_wiggle <- wiggleDF(pc_reads, group = "single", site = "cut5") %>%
    filter(cell_total > 100)	

  if(length(unique(sc_wiggle$set)) > 1){
    plt_sc_wiggle <- ggplot(sc_wiggle, aes(x=csd, y=CB, fill=(n/cell_total)))+
      geom_tile()+
      scale_fill_gradientn(colours=rev(brewer.pal(11,"RdYlBu")),na.value="grey40")+
      scale_y_discrete(expand=c(0,0))+
      scale_x_continuous(expand=c(0,0))+
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border=element_blank(),
            axis.line = element_blank())+
      facet_wrap(~set,nrow=1,strip.position = "top",scales="free_x")+
      ggtitle(paste0(figurePrefix, " sc cut5"))

    print(plt_sc_wiggle)
  }


  # scaled CDS/UTR

  hm_sc_scaled <- scaleDF(pc_reads, group = "single", site = "cut5", utr5 = lr_utr5, cds = lr_cds, utr3 = lr_utr3) %>% 
    filter(tot>100) %>%
    mutate(CB = fct_reorder(CB,tot)) %>%
    ggplot(aes(x=pos, y=CB, fill=(n/tot)))+
    geom_tile()+
    scale_fill_gradientn(colours=rev(brewer.pal(11,"RdYlBu")),na.value="grey40")+
    scale_y_discrete(expand=c(0,0))+
    scale_x_continuous(expand=c(0,0),
                       breaks=c(0,(lr_utr5/2),lr_utr5,(lr_utr5+lr_cds/2),(lr_utr5+lr_cds),(lr_utr5+lr_cds+(lr_utr3/2)),(lr_utr5+lr_cds+lr_utr3)),
                       labels=c("","5' UTR","","CDS","","3' UTR",""))+
    coord_equal()+
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border=element_blank())+
    xlab("distance along scaled metagene")+
    ggtitle(paste0(figurePrefix," CDS metagene"))

    print(hm_sc_scaled)


  ## cut location
  
  plt_sc_region <- pc_reads %>%
    mutate(site=cut5) %>%
    mutate(region = case_when( site < l_utr5 ~ "utr5",
                              ((site >= l_utr5) & (site < (l_utr5+l_cds) )) ~ "cds", 
                              (site >= (l_utr5 + l_cds)) ~ "utr3" )) %>%
    group_by(CB,region) %>%
    summarize(n=n()) %>%
    mutate(cell_total = sum(n)) %>%
    ggplot(aes(x=region,y=n/cell_total))+
    geom_boxplot()+
    ylab("fraction of reads per region")
  
  print(plt_sc_region)
  
  sc_length <- pc_reads %>%
    group_by(CB,length) %>%
    summarize(n=n()) %>%
    mutate(cell_total = sum(n)) %>%
    ungroup() %>%
    filter(cell_total>100)
  
  plt_sc_length_tot <- ggplot(sc_length,aes(x=length,y=n,color=CB))+
    geom_line()+
    theme(legend.position="none")+
    ylab('number of alignments')
  
  plt_sc_length_norm <- ggplot(sc_length,aes(x=length,y=n/cell_total,color=CB))+
    geom_line()+
    theme(legend.position="none")+
    ylab('fraction of alignments')
  
  title <- ggdraw()+draw_label(paste0(figurePrefix," length distribution"),x=0,hjust=0)+theme(plot.margin=margin(0,0,0,7))
  plt_sc_length_figure <- plot_grid(plt_sc_length_tot,plt_sc_length_norm,ncol=2)
  plt_sc_length <- plot_grid(title,plt_sc_length_figure,ncol=1,rel_heights=c(0.1,1))
  #save_plot(paste0(figurePrefix,"_sc_length_distribution.pdf"),plt_sc_length,base_height=4,base_width=8)
  print(plt_sc_length)


  ## plate layouts
  if(!is.null(args$whitelist)){
    if("row" %in% colnames(barcodes)){
      makepltlyt <- TRUE
    }else{
      makepltlyt <- FALSE
    }
  }else{
    makepltlyt <- FALSE
  }

  if(makepltlyt){
    plt_pl_reads <- reads %>%
      group_by(CB) %>% 
      summarize(n_reads = n()) %>%
      ungroup() %>%
      full_join(barcodes, by = "CB") %>%
      mutate(column = fct_reorder(as.factor(column), rank(as.numeric(column))),
             row = fct_reorder(as.factor(row), -rank(row))) %>%
      ggplot(aes(x=column, y=row))+
      geom_tile(aes(fill=log10(n_reads)))+
      coord_equal()+
      scale_fill_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), na.value = "grey50")+
      scale_x_discrete(expand = c(0,0), position = "bottom")+
      scale_y_discrete(expand = c(0,0))+
      panel_border(colour="grey20")+
      theme(legend.position = "bottom",
            axis.text = element_text(size=8),
            axis.line = element_blank(),
            axis.title = element_blank())

    plt_pl_pc_reads <- pc_reads %>%
      group_by(CB) %>% 
      summarize(pc_reads = n()) %>%
      ungroup() %>%
      full_join(barcodes, by = "CB") %>%
      mutate(column = fct_reorder(as.factor(column), rank(as.numeric(column))),
             row = fct_reorder(as.factor(row), -rank(row))) %>%
      ggplot(aes(x=column, y=row))+
      geom_tile(aes(fill=log10(pc_reads)))+
      coord_equal()+
      scale_fill_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), na.value = "grey50")+
      scale_x_discrete(expand = c(0,0), position = "bottom")+
      scale_y_discrete(expand = c(0,0))+
      panel_border(colour="grey20")+
      theme(legend.position = "bottom",
            axis.text = element_text(size=8),
            axis.line = element_blank(),
            axis.title = element_blank())

    plt_pl_region_frac <- pc_reads %>%
      mutate(site = cut5) %>%
      mutate(region = case_when(site < l_utr5 ~ "utr5",
                                ((site >= l_utr5) & (site < (l_utr5+l_cds) )) ~ "cds", 
                                (site >= (l_utr5 + l_cds)) ~ "utr3" )) %>%
      group_by(CB,region) %>%
      summarize(n=n()) %>%
      mutate(cell_total = sum(n)) %>%
      ungroup() %>%
      mutate(region_fraction = n/cell_total) %>%
      full_join(barcodes, by = "CB") %>%
      mutate(column = fct_reorder(as.factor(column), rank(as.numeric(column))),
             row = fct_reorder(as.factor(row), -rank(row))) %>%
      ggplot(aes(x=column, y=row))+
      geom_tile(aes(fill=region_fraction))+
      coord_equal()+
      facet_wrap(~region, ncol=2)+
      scale_fill_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), na.value = "grey50")+
      scale_x_discrete(expand = c(0,0), position = "bottom")+
      scale_y_discrete(expand = c(0,0))+
      panel_border(colour="grey20")+
      theme(legend.position = "bottom",
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank())

    plt_pl_frame5_frac <- pc_reads %>%
      group_by(CB,frame5) %>%
      summarize(n=n()) %>%
      mutate(cell_total = sum(n)) %>%
      ungroup() %>%
      mutate(frame5_fraction = n/cell_total) %>%
      full_join(barcodes, by = "CB") %>%
      mutate(column = fct_reorder(as.factor(column), rank(as.numeric(column))),
             row = fct_reorder(as.factor(row), -rank(row))) %>%
      ggplot(aes(x=column, y=row))+
      geom_tile(aes(fill=frame5_fraction))+
      coord_equal()+
      facet_wrap(~frame5, ncol=2)+
      scale_fill_gradientn(colours = rev(brewer.pal(11,"RdYlBu")), na.value = "grey50")+
      scale_x_discrete(expand = c(0,0), position = "bottom")+
      scale_y_discrete(expand = c(0,0))+
      panel_border(colour="grey20")+
      theme(legend.position = "bottom",
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank())

    plt_plates <- plot_grid(plt_pl_reads, plt_pl_pc_reads, plt_pl_region_frac, plt_pl_frame5_frac, nrow=2)
    print(plt_plates)
  }
}

dev.off()
