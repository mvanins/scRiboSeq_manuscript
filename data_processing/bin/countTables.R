#!/usr/bin/env Rscript
# Michael VanInsberghe 2020-05-22
# Creates count tables from read files, outputting unique counts per gene or per canonical transcript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(feather))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(R.utils))
#source("/hpc/hub_oudenaarden/mvanins/local/avopipelines/RPF/bin/scRibo/plotFunctions.R")

option_list = list(
                   make_option(c("-r", "--reads"), type = "character", default = NULL, help = "feather file containing read information dataframe", metavar="character"),
                   make_option(c("-a", "--annotations"), type = "character", default = NULL, help = "feather file containing annotation information", metavar="character"),
                   make_option(c("-p", "--prefix"), type = "character", default = NULL, help = "prefix for output filenames", metavar = "character"),
                   make_option(c("-s", "--site"), type = "character", default = "cut5", help = "site to use for overlaps", metavar = "character"),
                   make_option(c("-e", "--expansion"), type = "integer", default = 0, help = "number of baes to expand the CDS into the 5' and 3' UTRs"),
                   make_option(c("--rna"), action="store_true", default=FALSE, help = "Input is RNA-seq, count all alignments"),
                   make_option(c("--canonical"), action="store_true", default=FALSE, help = "keep reads in canonical set only"),
                   make_option(c("--strand"), default = "Forward", help = "Strand of reads to keep. Default 'Forward'. Allowed values are 'Forward', 'Reverse', 'All'"),
                   make_option(c("--minlength"), default = 30, help = "Minimum read length to count. Default 30. Only used for RPF datasets"),
                   make_option(c("--maxlength"), default = 45, help = "Maximum read length to count. Default 46. Only used for RPF datasets") )


opt_parser = OptionParser(option_list = option_list)
args = parse_args(opt_parser)

if(is.null(args$reads) | is.null(args$annotations)){
  print_help(opt_parser)
  stop("Please supply reads and annotations")
}

if(is.null(args$prefix)){
  countPrefix = strsplit(tools::file_path_sans_ext(basename(args$reads)),"_")[[1]][1]
}else{
  countPrefix = args$prefix
}

if(!(args$strand %in% c("Forward", "Reverse", "All"))){
  print_help(opt_parser)
  stop(paste0("Error! strand: ", args$strand, " is not one of Forward, Reverse, All"))
}

if(args$strand == "Forward"){
  keepStrand = "+"
}else if(args$strand == "Reverse"){
  keepStrand = "-"
}else{
  keepStrand = c("+", "-")
}


dir.create(countPrefix, showWarnings = FALSE)
print(countPrefix)

# annotations
if(grepl("feather",args$annotations)){
  annot <- read_feather(args$annotations)
}else{
  annot <- fread(args$annotations, stringsAsFactors = FALSE, data.table = FALSE)
}

# reads
if(grepl("feather",args$reads)){
  reads <- read_feather(args$reads)
}else{
  reads <- fread(args$reads)
}

reads <- reads %>%
  select(-starts_with("base"),
         -starts_with("frame"))

if("CB" %in% colnames(reads)){
  # single-cell libraries, remove reads with unassigned CB
  reads <- reads %>%
    filter(CB != "")
}else{
  # is bulk, add CB for grouping below
  reads <- reads %>%
    mutate(CB = countPrefix)
}


if(!(args$site %in% colnames(reads))){
  print_help(opt_parser)
  stop(paste0("Error! site:", args$site, "not found in reads dataframe"))
}

args

# for testing!
# annot <- read_feather("gencode.v33.chr_patch_hapl_scaff.annotation.annotations.feather")
# reads <- fread('RPFv4-Fucci-Mit03_Aligned.toTranscriptome.out.dedup.sorted.csv.gz') %>%
# 		select(-starts_with("base"), -starts_with("frame")) %>%
# 		filter(CB != "")
# args = data.frame(site = "cut5",
# 				  expansion = 0,
# 				  canonical = FALSE,
# 				  rna = FALSE,
# 				  minlength = 30,
# 				  maxlength = 46,
# 				  stringsAsFactors = FALSE)
# /for testing

if(args$canonical){
  countGroup = "transcript_id"
  expanded_pc <- annot %>%
    filter(set == "canonical",
           transcript_type == "protein_coding",
           transcript_id != "ENST00000437139.7")
}else{
  countGroup = "gene_id"
  expanded_pc <- annot %>%
    filter(transcript_type == "protein_coding",
           transcript_id != "ENST00000437139.7")
}

expanded_pc <- expanded_pc %>%
  mutate(l_utr5 = l_utr5 - args$expansion,
         l_cds = l_cds + 2*args$expansion,
         l_utr3 = l_utr3 - args$expansion ) %>%
  mutate(l_utr5 = case_when( l_utr5 < 0 ~ 0,
                            l_utr5 >= 0 ~ as.numeric(l_utr5)),
         l_cds = case_when(l_cds < 0 ~ 0,
                           l_cds >= 0 ~ as.numeric(l_cds)),
         l_utr3 = case_when(l_utr3 < 0 ~ 0,
                            l_utr3 >= 0 ~ as.numeric(l_utr3)))

if(!args$rna){
  coding_counts <- expanded_pc %>%
    inner_join(reads, by = "transcript_id") %>%
    filter(length == read_length,
           length >= args$minlength, 
           length <= args$maxlength,
           read_strand %in% keepStrand) %>%
    filter( (.data[[args$site]] >= (l_utr5)) & (.data[[args$site]] < (l_utr5 + l_cds))) 

}else{
  coding_counts <- expanded_pc %>%
    inner_join(reads, by = "transcript_id") %>%
    filter(length >= args$minlength,
           length <= args$maxlength,
           read_strand %in% keepStrand) 
}


coding_counts <- coding_counts %>%
  group_by(CB,.data[[countGroup]], id) %>%
  distinct(id) %>%
  ungroup() %>%
  group_by(CB, .data[[countGroup]]) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  rename(count_group = all_of(countGroup))

if(!args$rna){
  cds_fraction <- expanded_pc %>%
    inner_join(reads, by = "transcript_id") %>%
    filter(length == read_length,
           length >= args$minlength,
           length <= args$maxlength,
           read_strand %in% keepStrand) %>%
    group_by(CB, .data[[countGroup]], id) %>%
    distinct(id, .keep_all = TRUE) %>%
    ungroup() %>%
    mutate( region = case_when(.data[[args$site]] < (l_utr5) ~ "utr5",
                               ((.data[[args$site]]  >= (l_utr5)) & (.data[[args$site]] < ((l_utr5)+l_cds ) )) ~ "cds", 
                               (.data[[args$site]] >= (l_utr5 + l_cds)) ~ "utr3" )) %>%
    group_by(CB,region) %>%
    summarize(n_region = n()) %>%
    mutate(region = factor(region, levels = c('utr5', 'cds', 'utr3')))  %>%
    complete(region, CB, fill = list(n_region = 0)) %>%
    ungroup() %>%
    pivot_wider(names_from = "region", values_from = "n_region") %>%
    mutate(cds = replace_na(cds, 0),
           utr3 = replace_na(utr3, 0),
           utr5 = replace_na(utr5, 0)) %>%
    mutate(cds_frac = cds/(cds + utr3 + utr5),
           cell_total = (cds + utr3 + utr5)) %>%
    pivot_longer(cols = c(cds, utr3, utr5, cds_frac, cell_total), names_to = "count_group", values_to = "n")

  coding_counts <- bind_rows(coding_counts,
                             cds_fraction)
}

coding_counts <- coding_counts %>%
  mutate_if(is.character, as.factor)

# sparse matrix
# rows genes
# columns cells
# from https://www.biostars.org/p/312933/
sparseCounts <- with(coding_counts, sparseMatrix(i = as.numeric(count_group),
                                                 j = as.numeric(CB),
                                                 x = n,
                                                 dimnames = list(levels(count_group),
                                                                 levels(CB))))

writeMM(obj = sparseCounts, file = file.path(countPrefix, "matrix.mtx"))
write(x = rownames(sparseCounts), file = file.path(countPrefix, "features.tsv"))
write(x = colnames(sparseCounts), file = file.path(countPrefix, "barcodes.tsv"))

gzip(file.path(countPrefix, "matrix.mtx"), overwrite = TRUE, remove = TRUE)
gzip(file.path(countPrefix, "features.tsv"), overwrite = TRUE, remove = TRUE)
gzip(file.path(countPrefix, "barcodes.tsv"), overwrite = TRUE, remove = TRUE)
