# scRibo-seq
This repository contains the collection of scripts to analyze raw single-cell ribosomal profiling data and generate the figures presented in the manuscript.

The repository is organized in three subfolders:
#### 1. data_processing
Contains the NextFlow pipeline to process raw fastqs, creating QC plots, count tables, pre-processed alignments, and apply the random-forest correction.

#### 2. random_forest
Contains the R scripts to prepare, train, and tune a random-forest model to correct for the MNase sequence bias.

#### 3. figures
Contains the R scripts to main the subfigures for the manuscript.
