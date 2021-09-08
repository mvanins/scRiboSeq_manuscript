# Michael VanInsberghe 2021-03-15
# Compile, merge, and filter Fucci cells

library(dplyr)
library(tidyr)
library(tibble)
library(data.table)
library(Matrix)
library(Matrix.utils)
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

source('plotFunctions.R')

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


fucciSamples <- data.frame(names=c("G011",
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
                                   "SR206",
                                   "MGI01",
                                   "MGI02",
                                   "MGI03",
                                   "MGI04",
                                   "MGI05",
                                   "MGI06",
                                   "MGI07",
                                   "MGI08"),
                           counts = c(file.path(data_dir, "Fucci/counts/RPFv4-Fucci-G011"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-Fucci-G012"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-Fucci-Mit01"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-Fucci-Mit02"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-Fucci-Mit03"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-Fucci-Mit04"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-Fucci-Mit06"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-Fucci-Mit07"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-Fucci-Mit08"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-CDKN1AKO01"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-CDKN1AKO02"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-CDKN1AKO03"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-CDKN1AKO04"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-SR201"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-SR202"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-SR203"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-SR204"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-SR205"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-SR206"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-MGI01"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-MGI02"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-MGI03"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-MGI04"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-MGI05"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-MGI06"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-MGI07"),
                                      file.path(data_dir, "Fucci/counts/RPFv4-MGI08")),
                           allcounts = c(file.path(data_dir, "Fucci/all/counts/RPFv4-Fucci-G011"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-Fucci-G012"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-Fucci-Mit01"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-Fucci-Mit02"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-Fucci-Mit03"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-Fucci-Mit04"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-Fucci-Mit06"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-Fucci-Mit07"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-Fucci-Mit08"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-CDKN1AKO01"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-CDKN1AKO02"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-CDKN1AKO03"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-CDKN1AKO04"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-SR201"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-SR202"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-SR203"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-SR204"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-SR205"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-SR206"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-MGI01"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-MGI02"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-MGI03"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-MGI04"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-MGI05"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-MGI06"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-MGI07"),
                                         file.path(data_dir, "Fucci/all/counts/RPFv4-MGI08")),
                           biotypes = c(file.path(data_dir, "Fucci/biotypes/RPFv4-Fucci-G011"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-Fucci-G012"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-Fucci-Mit01"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-Fucci-Mit02"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-Fucci-Mit03"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-Fucci-Mit04"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-Fucci-Mit06"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-Fucci-Mit07"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-Fucci-Mit08"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-CDKN1AKO01"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-CDKN1AKO02"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-CDKN1AKO03"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-CDKN1AKO04"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-SR201"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-SR202"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-SR203"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-SR204"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-SR205"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-SR206"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-MGI01"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-MGI02"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-MGI03"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-MGI04"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-MGI05"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-MGI06"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-MGI07"),
                                        file.path(data_dir, "Fucci/biotypes/RPFv4-MGI08")),
                           allbiotypes = c(file.path(data_dir, "Fucci/all/biotypes/RPFv4-Fucci-G011"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-Fucci-G012"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-Fucci-Mit01"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-Fucci-Mit02"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-Fucci-Mit03"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-Fucci-Mit04"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-Fucci-Mit06"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-Fucci-Mit07"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-Fucci-Mit08"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-CDKN1AKO01"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-CDKN1AKO02"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-CDKN1AKO03"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-CDKN1AKO04"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-SR201"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-SR202"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-SR203"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-SR204"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-SR205"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-SR206"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-MGI01"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-MGI02"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-MGI03"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-MGI04"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-MGI05"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-MGI06"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-MGI07"),
                                           file.path(data_dir, "Fucci/all/biotypes/RPFv4-MGI08")),
                           contamination = c(file.path(data_dir, "Fucci/contamination/RPFv4-Fucci-G011_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-Fucci-G012_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-Fucci-Mit01_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-Fucci-Mit02_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-Fucci-Mit03_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-Fucci-Mit04_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-Fucci-Mit06_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-Fucci-Mit07_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-Fucci-Mit08_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-CDKN1AKO01_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-CDKN1AKO02_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-CDKN1AKO03_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-CDKN1AKO04_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-SR201_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-SR202_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-SR203_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-SR204_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-SR205_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-SR206_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-MGI01_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-MGI02_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-MGI03_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-MGI04_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-MGI05_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-MGI06_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-MGI07_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz"),
                                             file.path(data_dir, "Fucci/contamination/RPFv4-MGI08_contamination_Aligned.sortedByCoord.out_contamination_counts.csv.gz")),
                           facs = c(file.path(data_dir, "Fucci/facs/G011.csv"),
                                    file.path(data_dir, "Fucci/facs/G012.csv"),
                                    file.path(data_dir, "Fucci/facs/Mit01.csv"),
                                    file.path(data_dir, "Fucci/facs/Mit02.csv"),
                                    file.path(data_dir, "Fucci/facs/Mit03.csv"),
                                    file.path(data_dir, "Fucci/facs/Mit04.csv"),
                                    file.path(data_dir, "Fucci/facs/Mit06.csv"),
                                    file.path(data_dir, "Fucci/facs/Mit07.csv"),
                                    file.path(data_dir, "Fucci/facs/Mit08.csv"),
                                    file.path(data_dir, "Fucci/facs/KO01.csv"),
                                    file.path(data_dir, "Fucci/facs/KO02.csv"),
                                    file.path(data_dir, "Fucci/facs/KO03.csv"),
                                    file.path(data_dir, "Fucci/facs/KO04.csv"),
                                    file.path(data_dir, "Fucci/facs/SR201.csv"),
                                    file.path(data_dir, "Fucci/facs/SR202.csv"),
                                    file.path(data_dir, "Fucci/facs/SR203.csv"),
                                    file.path(data_dir, "Fucci/facs/SR204.csv"),
                                    file.path(data_dir, "Fucci/facs/SR205.csv"),
                                    file.path(data_dir, "Fucci/facs/SR206.csv"),
                                    file.path(data_dir, "Fucci/facs/MGI01.csv"),
                                    file.path(data_dir, "Fucci/facs/MGI02.csv"),
                                    file.path(data_dir, "Fucci/facs/MGI03.csv"),
                                    file.path(data_dir, "Fucci/facs/MGI04.csv"),
                                    file.path(data_dir, "Fucci/facs/MGI05.csv"),
                                    file.path(data_dir, "Fucci/facs/MGI06.csv"),
                                    file.path(data_dir, "Fucci/facs/MGI07.csv"),
                                    file.path(data_dir, "Fucci/facs/MGI08.csv")),
                           reads = c(file.path(data_dir, "Fucci/reads/RPFv4-Fucci-G011_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-Fucci-G012_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-Fucci-Mit01_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-Fucci-Mit02_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-Fucci-Mit03_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-Fucci-Mit04_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-Fucci-Mit06_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-Fucci-Mit07_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-Fucci-Mit08_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-CDKN1AKO01_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-CDKN1AKO02_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-CDKN1AKO03_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-CDKN1AKO04_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-SR201_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-SR202_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-SR203_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-SR204_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-SR205_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-SR206_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-MGI01_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-MGI02_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-MGI03_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-MGI04_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-MGI05_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-MGI06_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-MGI07_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz"),
                                     file.path(data_dir, "Fucci/reads/RPFv4-MGI08_Aligned.toTranscriptome.out.dedup.sorted.csv.predicted.csv.gz")),
                           stringsAsFactors = FALSE)


## cells meeting cutoff filter, but are much lower than cells in similar phase of cell cycle 
s_low <- c("G011_GTACTCCGTC", "Mit01_ATCCTAACAG", "Mit01_TTGTCTCTAC", "Mit02_GACCTACGAG", "Mit02_TAAGCCTCTA", "Mit02_TGGTGTTGAC", "Mit03_ACAGAGACCG",
           "Mit03_AGTTCTTCCG", "Mit03_CCGACAAGGA", "Mit03_TCAGAGGCAC", "Mit04_AACCGTAACA", "Mit04_CTAATCGCGG", "Mit04_GTCTTGGAAT", "Mit04_TCTCACTCTA",
           "Mit06_CCGTGAGAAC", "Mit06_CGATACTGAA", "Mit06_GAACCTCCTA", "Mit06_GCAGCGTTAG", "Mit06_TAAGCCTCTA", "Mit07_CAACTGTCAT", "Mit07_TTGCTCACGG",
           "Mit08_AATTCGGTGG", "Mit08_CAACTGTCAT", "KO01_CCGACAAGGA", "SR201_CGAGACTTCA", "SR201_GAAGTCGTGG", "SR201_TAGGTGGTGG", "SR201_TCAAGGTCGC",
           "SR202_ACGCGCAACA", "SR202_AGTCACAACA", "SR202_CGTCGGTAAG", "SR202_TTGGTGTGTC", "SR204_TATCCGATCT", "KO04_AATGCTCTAG", "SR204_CGTTGTTGCG",
           "MGI01_TTCATCCGTG", "MGI02_GTGTTCTGTT", "MGI03_AATCTGTTGG", "MGI03_ACTGGACGTT", "MGI03_GTTAGCACAA", "MGI04_AACTCCGAGA", 
           "MGI04_ACTTGCGTGT", "MGI04_ATGCTTCTAC", "MGI04_CGGCTGTTAG", "MGI05_AATCTGTTGG", "MGI05_ACGTATCTTG", "MGI05_AGAACAATCG", 
           "MGI05_CCGTAACCTT", "MGI05_GCTCTCCTTC", "MGI06_GATCACTGGT", "MGI07_ATCTGTGAGG", "MGI07_CCTCATCGTT", "MGI07_GTTCTAATGG") 

fucciFacs <- readFACS(fucciSamples, barcodes = barcodes) %>%
  mutate(green = X..488..530.40,
         red = X..561..585.29)

fucciReads <- canonical_pc %>%
  inner_join(fread(file.path(data_dir, "Fucci/reads/mergedFilteredReads.csv.gz")), by = "transcript_id")

fucciCounts <- readMergeCounts(fucciSamples)

fucciBiotypes <- readMergeBiotypes(fucciSamples)


fucciMeta <- data.frame(row.names = colnames(fucciCounts),
                        CB = colnames(fucciCounts),
                        well = barcodes[gsub(".*_(.*?)$","\\1",colnames(fucciCounts),perl=TRUE),"well"],
                        plate = gsub("^(.*?)[_].*","\\1",colnames(fucciCounts),perl=TRUE),
                        cds = fucciCounts["cds",],
                        utr3 = fucciCounts["utr3",],
                        utr5 = fucciCounts["utr5",],
                        cds_frac = fucciCounts["cds_frac",],
                        cell_total = fucciCounts["cell_total",],
                        sort_population = rep("Interphase",ncol(fucciCounts)),
                        type = rep("RPF", ncol(fucciCounts)),
                        stringsAsFactors = FALSE)

fucciMeta[ (grepl("^G0", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(1:8, LETTERS[1:16])), "sort_population"] = "G0"
fucciMeta[ (grepl("^Mit", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(1:8, LETTERS[1:16])), "sort_population"] = "Mit"
fucciMeta[ (grepl("^Mit", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(21:24, LETTERS[15:16])), "sort_population"] = "NTC"
fucciMeta[ (grepl("^G0", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(21:24, LETTERS[15:16])), "sort_population"] = "NTC"

fucciMeta[ (grepl("^KO", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(1:4, LETTERS[1:16]) ), "sort_population"] = "KO17"
fucciMeta[ (grepl("^KO", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(5:8, LETTERS[1:16]) ), "sort_population"] = "KO15"
fucciMeta[ (grepl("^KO", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(9:12, LETTERS[1:16]) ), "sort_population"] = "KO14"
fucciMeta[ (grepl("^KO", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(13:16, LETTERS[1:16]) ), "sort_population"] = "KO13"
fucciMeta[ (grepl("^KO", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(17:24, LETTERS[1:16]) ), "sort_population"] = "WT"
fucciMeta[ (grepl("^KO", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(22:24, LETTERS[15:16]) ), "sort_population"] = "NTC"

fucciMeta[ (grepl("^SR2", fucciMeta$plate) & fucciMeta$well %in% c(columnFeatures(1:4, LETTERS[1:16]), columnFeatures(5, LETTERS[1:8])) ), "sort_population"] = "2h_post"
fucciMeta[ (grepl("^SR2", fucciMeta$plate) & fucciMeta$well %in% c(columnFeatures(6:9, LETTERS[1:16]), columnFeatures(5, LETTERS[9:16])) ), "sort_population"] = "6h_post"
fucciMeta[ (grepl("^SR2", fucciMeta$plate) & fucciMeta$well %in% c(columnFeatures(10:13, LETTERS[1:16]), columnFeatures(14, LETTERS[1:8])) ), "sort_population"] = "10h_post"
fucciMeta[ (grepl("^SR2", fucciMeta$plate) & fucciMeta$well %in% c(columnFeatures(15:18, LETTERS[1:16]), columnFeatures(14, LETTERS[9:16])) ), "sort_population"] = "20h_post"
fucciMeta[ (grepl("^SR2", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(19:21, LETTERS[1:16]) ), "sort_population"] = "cycling_post"
fucciMeta[ (grepl("^SR2", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(22:24, LETTERS[1:14]) ), "sort_population"] = "cycling_pre"
fucciMeta[ (grepl("^SR2", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(22:24, LETTERS[15:16]) ), "sort_population"] = "NTC"

fucciMeta[ (grepl("^MGI", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(1:15, LETTERS[1:16])), "sort_population"] = "Interphase"
fucciMeta[ (grepl("^MGI", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(16:21, LETTERS[1:16])), "sort_population"] = "Mit"
fucciMeta[ (grepl("^MGI", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(22:24, LETTERS[1:14])), "sort_population"] = "G0"
fucciMeta[ (grepl("^MGI", fucciMeta$plate) & fucciMeta$well %in% columnFeatures(22:24, LETTERS[15:16])), "sort_population"] = "NTC"

fucciMeta[ (grepl("^KO", fucciMeta$plate)), "batch"] = "RPF_KO"
fucciMeta[ (grepl("Mit|G0", fucciMeta$plate)), "batch"] = "RPF_First"
fucciMeta[ (grepl("^SR2", fucciMeta$plate)), "batch"] = "RPF_SR2"
fucciMeta[ (grepl("^MGI", fucciMeta$plate)), "batch"] = "RPF_MGI"

## Are there any cells with two FACS entries? Discard
multipleCells <- pull(fucciFacs %>%
  group_by(CB) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  filter(n>1) %>%
  select(CB))

fucciFacs <- fucciFacs %>%
  filter(!(CB %in% multipleCells))
fucciMeta <- fucciMeta %>% 
  filter(!(CB %in% multipleCells))


## add FACS information to fucci metadata
fucciMeta <- fucciMeta %>%
  full_join(fucciFacs %>% filter(CB %in% fucciMeta$CB) %>% select(-Well, -sample), by = "CB")
rownames(fucciMeta) <- fucciMeta$CB
fucciMeta <- fucciMeta[colnames(fucciCounts),]

fucciCounts <- fucciCounts[!(rownames(fucciCounts) %in% c("cds","cds_frac","cell_total","utr3","utr5")),]

fucciCells <- pull(fucciMeta %>% 
                 filter(cds_frac >= 0.87, cell_total > 4500, cell_total < 50000) %>%
                 filter(!(sort_population %in% c("NTC"))) %>%
                 filter(!(batch %in% c("Sync"))) %>%
                 filter(!(sort_population %in% c("KO14", "KO13", "KO15", "KO17",
                                       "9h_post", "5h_post", "2h_post", "20h_post", "10h_post", "6h_post"))) %>%
                 filter(!(CB %in% s_low)) %>%
                 select(CB))

saveRDS(fucciMeta, file.path(data_dir, "Fucci/compiled/allmeta.rds"))

fucciReads <- fucciReads %>%
  filter(CB %in% fucciCells)
fucciCounts <- fucciCounts[,fucciCells]
fucciBiotypes <- fucciBiotypes[,fucciCells]
fucciMeta <- fucciMeta[fucciCells,]

saveRDS(fucciReads, file.path(data_dir, "Fucci/compiled/reads.rds"))
saveRDS(fucciCounts, file.path(data_dir, "Fucci/compiled/counts.rds"))
saveRDS(fucciMeta, file.path(data_dir, "Fucci/compiled/meta.rds"))
saveRDS(fucciBiotypes, file.path(data_dir, "Fucci/compiled/biotypes.rds"))
saveRDS(fucciSamples, file.path(data_dir, "Fucci/compiled/samples.rds"))
