#!/usr/bin/env Rscript
# Michael VanInsberghe 2021-03-05
# counts biotypes per cell/library

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
                   make_option(c("--rna"), action="store_true", default=FALSE, help = "Input is RNA-seq, count all alignments"),
                   make_option(c("--canonical"), action="store_true", default=FALSE, help = "keep reads in canonical set only"),
                   make_option(c("--strand"), default = "Forward", help = "Strand of reads to keep. Default 'Forward'. Allowed values are 'Forward', 'Reverse', 'All'"),
                   make_option(c("--minlength"), default = 30, help = "Minimum read length to count. Default 30. Only used for RPF datasets"),
                   make_option(c("--maxlength"), default = 45, help = "Maximum read length to count. Default 45. Only used for RPF datasets") )


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

# for testing!
#annot <- fread("/hpc/hub_oudenaarden/mvanins/local/lib/reference/gencode/hsa_33/RPF_conda_reference/annotations/gencode.v33.chr_patch_hapl_scaff.annotation.annotations.csv.gz", stringsAsFactors = FALSE, data.table = FALSE)
#reads <- fread('/hpc/hub_oudenaarden/mvanins/RPF/data/scRPF_Fucci/RPF_update/scRibo/RPFv4-MGI01_Aligned.toTranscriptome.out.csv.gz') %>%
#  select(-starts_with("base"), -starts_with("frame")) %>%
#  filter(CB != "")
#args = data.frame(rna = FALSE,
#                  canonical = FALSE,
#                  minlength = 30,
#                  maxlength = 46,
#                  stringsAsFactors = FALSE)
#keepStrand = "+"
# /for testing



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

args

if(args$canonical){
  countGroup = "transcript_type"
  expanded_pc <- annot %>%
    filter(set == "canonical",
           gene_id != "ENSMUSG00000098178.1",
           transcript_id != "ENST00000437139.7")
}else{
  countGroup = "gene_type"
  expanded_pc <- annot %>%
    filter(gene_id != "ENSMUSG00000098178.1",
           transcript_id != "ENST00000437139.7")
}

if(!args$rna){
  biotype_counts <- expanded_pc %>%
    inner_join(reads, by = "transcript_id") %>%
    filter(length == read_length,
           length >= args$minlength, 
           length <= args$maxlength,
           read_strand %in% keepStrand)

}else{
  biotype_counts <- expanded_pc %>%
    inner_join(reads, by = "transcript_id") %>%
    filter(length >= args$minlength,
           length <= args$maxlength,
           read_strand %in% keepStrand) 
}

biotype_counts <- biotype_counts %>%
  group_by(CB,.data[[countGroup]], id) %>%
  distinct(id) %>%
  ungroup() %>%
  group_by(CB, .data[[countGroup]]) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  rename(count_group = all_of(countGroup)) %>%
  mutate_if(is.character, as.factor)


sparseCounts <- with(biotype_counts, sparseMatrix(i = as.numeric(count_group),
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
