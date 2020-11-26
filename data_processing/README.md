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
