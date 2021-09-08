# scRibo-seq
This repository contains the collection of scripts to analyze raw single-cell ribosomal profiling data and generate the figures presented in the manuscript.

VanInsberghe, M., van den Berg, J., Andersson-Rolf, A., Clevers, H., van Oudenaarden, A., Single-cell Ribo-seq reveals cell cycle-dependent translational pausing. *Nature* (2021). https://doi.org/10.1038/s41586-021-03887-4



The repository is organized in three subfolders:
#### 1. data_processing
Contains the NextFlow pipeline to process raw fastqs, creating QC plots, count tables, pre-processed alignments, and apply the random-forest correction.

#### 2. random_forest
Contains the R scripts to prepare, train, and tune a random-forest model to correct for the MNase sequence bias. Additionally contains the tuned model used to correct the MNase sequence bias.

#### 3. figures
Contains the R scripts to main the subfigures for the manuscript.

#### 4. processing_config
Contains the nextflow configuration files and any modifications to the standard workflow for each data set.
