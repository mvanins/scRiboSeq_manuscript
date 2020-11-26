# RPF
Nextflow pipeline for processing (sc)RPF data

### Basic usage
Create a nextflow.config in your data directory, setting, e.g.,:
```
params {
  genomeDir = "/hpc/genome_directory"
  parseReads = "--id"
  ranger_model = "fucci.tuned.model_30-45.RData"
}
```

and running `nextfow run main.nf -resume`


## References
Some of the scripts for tRNA masking are from Hoffmann, A. et al. Accurate mapping of tRNA reads. Bioinformatics 34, 1116-1124, doi:10.1093/bioinformatics/btx756 (2018).
https://github.com/AnneHoffmann/tRNA-read-mapping
