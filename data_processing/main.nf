#!/usr/bin/env nextflow

/*
Ribosome footprint analysis pipeline
*/

def check_max(obj, type) {
  if (type == 'memory') {
    try {
      if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
        return params.max_memory as nextflow.util.MemoryUnit
      else
        return obj
    } catch (all) {
      println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
        return obj
    }
  } else if (type == 'time') {
    try {
      if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
        return params.max_time as nextflow.util.Duration
      else
        return obj
    } catch (all) {
      println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
        return obj
    }
  } else if (type == 'cpus') {
    try {
      return Math.min( obj, params.max_cpus as int )
    } catch (all) {
      println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
        return obj
    }
  }
}

if(!params.genomeDir) {
  // build geneome directory if not supplied
  
  referenceOut = "${params.outdir}/reference"
  annotations = file(params.annotations, checkIfExists: true)
  genomeFasta = file(params.genomeFasta, checkIfExists: true)
  contaminationFasta = file(params.contaminationFasta, checkIfExists: true)

  process index_genome {

    publishDir "${referenceOut}/original/", mode:'copy', overwrite: true

    tag "indexing $genome"

    input:
    file(genomeFasta)
    file(annotations)

    output:
    set file(genomeFasta), file("*.fai") into genome_euk, genome_mit, genome_prep, genome_scribo, genome_star
    file("${annotations}")
    file("*.fai") into genomeFai
    file("${genomeFasta}") into genomeFasta

    script:
    """
    samtools faidx ${genomeFasta}
    """
  }

  process eukaryotic_tRNAscan {
    // follows the recommended paramaters for the updated tRNAscan-SE 
    // from https://dx.doi.org/10.1007%2F978-1-4939-9173-0_1
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6768409/
    
    cpus = { check_max( 36 * task.attempt, 'cpus' ) }
    memory = { check_max( 10.GB * task.attempt, 'memory' ) }
    time = { check_max( 16.h * task.attempt, 'time' ) }

    publishDir "${referenceOut}/tRNAScan/eukaryotic/", mode:'copy', overwrite: true
    
    tag "Eukaryotic tRNAScan on $fasta"
  
    input:
    set file(genome), file(fai) from genome_euk
  
    output:
    file("${genome.baseName}.Eukaryotic_tRNAs_bp.bed") into eukar_tRNAs_bed
    file("*")
  
    script:
    """
    tRNAscan-SE \
      -HQ \
      -o# \
      -f# \
      -s# \
      -m# \
      -b# \
      -a# \
      -l# \
      --brief \
      --thread ${task.cpus} \
      -p ${genome.baseName}.Eukaryotic_tRNAs \
      ${genome}
  
    cat ${genome.baseName}.Eukaryotic_tRNAs.out | tRNAScanToBed.pl > ${genome.baseName}.Eukaryotic_tRNAs_bp.bed
    sed '/^MT\\|^chrM/d' ${genome.baseName}.Eukaryotic_tRNAs_bp.bed > temp.bed
    mv temp.bed ${genome.baseName}.Eukaryotic_tRNAs_bp.bed
    """
  }
  
  process mitochondrial_tRNAscan {
    // follows the recommended paramaters for the updated tRNAscan-SE 
    // from https://dx.doi.org/10.1007%2F978-1-4939-9173-0_1

    cpus = { check_max( 36 * task.attempt, 'cpus' ) }
    memory = { check_max( 10.GB * task.attempt, 'memory' ) }
    time = { check_max( 16.h * task.attempt, 'time' ) }
    
    publishDir "${referenceOut}/tRNAScan/mitochondrial/", mode:'copy', overwrite: true  
  
    tag "Mitochondrial tRNAScan on $fasta"
  
    input:
    set file(genome), file(fai) from genome_mit
  
    output:
    file("${genome.baseName}.Mitochondrial_tRNAs_bp.bed") into mito_tRNAs_bed
    file("*")
  
    script:
    """
    tRNAscan-SE \
      -M vert \
      -Q \
      -o# \
      -f# \
      -m# \
      -b# \
      -a# \
      -l# \
      --brief \
      --thread ${task.cpus} \
      -p ${genome.baseName}.Mitochondrial_tRNAs \
      ${genome}

    cat ${genome.baseName}.Mitochondrial_tRNAs.out | tRNAScanToBed.pl > ${genome.baseName}.Mitochondrial_tRNAs_bp.bed
    """
  }
  
  process create_genomes{

    cpus = { check_max( 1 * task.attempt, 'cpus' ) }
    memory = { check_max( 5.GB * task.attempt, 'memory' ) }
    time = { check_max( 6.h * task.attempt, 'time' ) }
    clusterOptions '--gres=tmpspace:75G' // for sort in mergeGTF
  
    tag "create tRNA-masked genomes for $genome";
  
    publishDir "${referenceOut}/tRNA-masked/", mode:'copy', overwrite: true
  
    input:
    file(eukaryotic) from eukar_tRNAs_bed
    file(mitochondrial) from mito_tRNAs_bed
    set file(genome), file(genome_index) from genome_prep
    file(annotations) from annotations
  
    output:
    file("${genome.baseName}.tRNA_masked_regions.bed")
    // file("${genome.baseName}.tRNA_masked.fa")
    file("${genome.baseName}.pre-tRNAs.bed")
    file("${genome.baseName}.pre-tRNAs.fa")
    file("${genome.baseName}.mature-tRNAs.bed")
    file("${genome.baseName}.mature-tRNAs.fa")
    file("${genome.baseName}.mature-tRNAs.cluster.fa")
    file("${genome.baseName}.mature-tRNAs.cluster.fa.fai")
    file("${genome.baseName}.tRNACluster.gtf")
    file("${genome.baseName}.tRNA-masked-withmature.fa")
    file("${genome.baseName}.tRNA-masked-withmature.fa.fai")
    // miR file("${annotations.baseName}.miRNA.tRNA.gtf")
    file("${annotations.baseName}.tRNA.gtf")
    set file("${genome.baseName}.tRNA-masked-withmature.fa"), file("${annotations.baseName}.tRNA.gtf") into star_genome
    set file("${genome.baseName}.mature-tRNAs.cluster.fa"), file("${genome.baseName}.tRNACluster.gtf") into star_tRNA_contaminant_genome
  
    script:
    """
    #### tRNAs
    # mask tRNAs in genome
    cat ${eukaryotic} ${mitochondrial} | sort -k 1,1 -k2,2n  > ${genome.baseName}.tRNA_masked_regions.bed
    bedtools maskfasta -fi ${genome} -fo ${genome.baseName}.tRNA_masked.fa -bed ${genome.baseName}.tRNA_masked_regions.bed
  
    # create genome file
    cut -f 1-2 ${genome_index} > ${genome_index.baseName}.genome

    ## pre tRNAs
    # remove pseudogenes, add 50 nt 5' and 3' flanking regions
    grep -v "pseudo\\s" ${eukaryotic} | expandBed12.pl --slop 50 --index ${genome_index} > ${eukaryotic.baseName}.pre-tRNAs.bed
  
    # select only MT tRNAs from the MT model, add flaking regions
    sed '/^MT\\|^chrM/!d' ${mitochondrial} > ${mitochondrial.baseName}.only_MT.bed
    grep -v "pseudo\\s" ${mitochondrial.baseName}.only_MT.bed | expandBed12.pl --slop 50 --index ${genome_index} > ${mitochondrial.baseName}.pre-tRNAs.bed
  
    # combine to pre-tRNAs.bed
    cat ${eukaryotic.baseName}.pre-tRNAs.bed ${mitochondrial.baseName}.pre-tRNAs.bed | sort -k 1,1 -k2,2n > ${genome.baseName}.pre-tRNAs.bed
  
    # extract pre-tRNA sequences
    bedtools getfasta -name -split -s -fi ${genome} -bed ${genome.baseName}.pre-tRNAs.bed -fo ${genome.baseName}.pre-tRNAs.fa
    
    ## mature tRNAs
    cat ${eukaryotic} ${mitochondrial.baseName}.only_MT.bed | sort -k 1,1 -k2,2n > ${genome.baseName}.mature-tRNAs.bed
    bedtools getfasta -name -split -s -fi ${genome} -bed ${genome.baseName}.mature-tRNAs.bed | appendFasta.pl --append cca > ${genome.baseName}.mature-tRNAs.fa
  
    ## cluster mature tRNAs
    collapseSequences.pl ${genome.baseName}.mature-tRNAs.fa > ${genome.baseName}.mature-tRNAs.cluster.fa
    # produces cluster_info.gtf
    mv cluster_info.gtf ${genome.baseName}.tRNACluster.gtf
    samtools faidx ${genome.baseName}.mature-tRNAs.cluster.fa
  
    ## assemble final genome
    cat ${genome.baseName}.tRNA_masked.fa ${genome.baseName}.mature-tRNAs.cluster.fa > ${genome.baseName}.tRNA-masked-withmature.fa
    samtools faidx ${genome.baseName}.tRNA-masked-withmature.fa
  
  
    #### Annotations
    mergeGTF.pl ${annotations} ${genome.baseName}.tRNACluster.gtf > ${annotations.baseName}.tRNA.gtf
  
    # geneinfo file for later analysis
    #extractGeneInfo.pl ${annotations.baseName}.tRNA.gtf > ${annotations.baseName}.tRNA.geneinfo
    """
  }

  process scRibo_annotations{

    cpus = { check_max( 1 * task.attempt, 'cpus' ) }
    memory = { check_max( 30.GB * task.attempt, 'memory' ) }
    time = { check_max( 6.h * task.attempt, 'time' ) }

    tag "scRibo annotation preprocessing for $genome"
    
    publishDir "${referenceOut}/annotations/", mode:'copy', overwrite: true

    input:
    set file(genome), file(genome_index) from genome_scribo
    file(annotations) from annotations

    output:
    file(annotations)
    file("*.{csv.gz,feather}")
    file("*.pickle")
    file("${annotations.baseName}.annotations.{csv.gz,feather}") into annotationFeather
    file("*.annotations.pickle") into annotationPickle
    

    script:
    """
    extractAnnotations.py -g ${annotations}
    extractCodons.py -f ${genome} -t ${annotations.baseName}.annotations.pickle
    """
}
  
  process star_reference{

    cpus = { check_max( 24 * task.attempt, 'cpus' ) }
    memory = { check_max( 75.GB * task.attempt, 'memory' ) }
    time = { check_max( 8.h * task.attempt, 'time' ) }
  
    tag "STAR reference"
  
    publishDir "${referenceOut}/", mode:'copy', overwrite: true
  
    input:
    set file(masked_genome), file(annotations) from star_genome
  
    output:
    file("star_tRNAmasked_${params.starOverhang}")
    file("star_tRNAmasked_${params.starOverhang}") into starIndex
  
    script:
    """
    mkdir ./star_tRNAmasked_${params.starOverhang}
  
    STAR \
      --runMode genomeGenerate \
      --runThreadN ${task.cpus} \
      --genomeDir ./star_tRNAmasked_${params.starOverhang} \
      --genomeFastaFiles ${masked_genome} \
      --limitGenomeGenerateRAM 75161927680 \
      --sjdbGTFfile ${annotations} \
      --sjdbOverhang ${params.starOverhang}
    """
  }
  
  process star_contamination_reference{

    cpus = { check_max( 8 * task.attempt, 'cpus' ) }
    memory = { check_max( 30.GB * task.attempt, 'memory' ) }
    time = { check_max( 8.h * task.attempt, 'time' ) }
  
    tag "STAR contamination reference"
  
    publishDir "${referenceOut}/contamination/", mode:'copy', overwrite: true
  
    input:
    set file(tRNA_genome), file(tRNA_annotations) from star_tRNA_contaminant_genome
    file(contamination_fasta) from contaminationFasta
  
    output:
    file("star_contamination_${params.starOverhang}")
    file('contamination.fa')
    file('contamination.gtf')
    file("star_contamination_${params.starOverhang}") into starContaminationIndex
  
    script:
    """
    # prepare contamination gtf
    fastaToGTF.pl ${contamination_fasta} > ${contamination_fasta.baseName}.gtf

    cat ${tRNA_genome} ${contamination_fasta} > contamination.fa
    samtools faidx contamination.fa

    mergeGTF.pl ${tRNA_annotations} ${contamination_fasta.baseName}.gtf > contamination.gtf

    mkdir ./star_contamination_${params.starOverhang}
  
    STAR \
      --runMode genomeGenerate \
      --runThreadN ${task.cpus} \
      --genomeDir ./star_contamination_${params.starOverhang} \
      --genomeFastaFiles contamination.fa \
      --limitGenomeGenerateRAM 32212254720\
      --sjdbGTFfile contamination.gtf \
      --sjdbOverhang ${params.starOverhang}
    """
  }
  
}else{
  // genome exists, set variables
  annotationFeather = file("${params.genomeDir}/annotations/*.annotations.{csv.gz,feather}", checkIfExists: true).first()
  annotationPickle = file("${params.genomeDir}/annotations/*.annotations.pickle", checkIfExists: true).first()
  genomeFasta = file("${params.genomeDir}/original/*.{fa,fasta}", checkIfExists: true).first()
  genomeFai =  file("${params.genomeDir}/original/*.fai", checkIfExists: true).first()
  starIndex = file("${params.genomeDir}/star_tRNAmasked_${params.starOverhang}", checkIfExists: true)
  starContaminationIndex = file("${params.genomeDir}/contamination/star_contamination_${params.starOverhang}", checkIfExists: true)
  
  //Channel
  //  .fromPath("${params.genomeDir}/annotations/*.annotations.{csv.gz,feather}")
  //  .first()
  //  .ifEmpty{ exit 1, "Cannot find annotation csv/feather in: ${params.genomeDir}/annotations\nGenome directory not properly formed" }
  //  .set{annotationFeather}

  //Channel
  //  .fromPath("${params.genomeDir}/annotations/*.annotations.pickle")
  //  .first()
  //  .ifEmpty{ exit 1, "Cannot find annotation pickle in: ${params.genomeDir}/annotations\nGenome directory not properly formed" }
  //  .set{annotationPickle}

  //Channel
  //  .fromPath("${params.genomeDir}/original/*.{fa,fasta}")
  //  .first()
  //  .ifEmpty{ exit 1, "Cannot find genome fasta in ${params.genomeDir}/original\nGenome directory not properly formed" }
  //  .set{genomeFasta}

  //Channel
  //  .fromPath("${params.genomeDir}/original/*.fai")
  //  .first()
  //  .ifEmpty{ exit 1, "Cannot find genome index in ${params.genomeDir}/original\nGenome directory not properly formed" }
  //  .set{genomeFai}

  //Channel
  //  .fromPath("${params.genomeDir}/star_tRNAmasked_${params.starOverhang}")
  //  .ifEmpty{ exit 1, "Cannot find STAR index with overhang ${params.starOverhang} in: ${params.genomeDir}\nGenome directory not properly formed" }
  //  .set{starIndex}


  //  if(params.deplete){
  //    Channel
  //      .fromPath("${params.genomeDir}/rRNA_depletion")
  //      .ifEmpty{ exit 1, "Cannot find rRNA depletion index in ${params.genomeDir}\nGenome directory not properly formed" }
  //      .set(rRNAIndex)

  //    Channel
  //      .fromrath("${params.genomeDir}/rRNA_depletion")
  //      .ifEmpty{ exit 1, "Cannot find rRNA depletion index in ${params.genomeDir}\nGenome directory not properly formed" }
  //      .set(tRNAIndex)
  //  }
}


if(!params.bulk){
  whitelist = file(params.whitelist, checkIfExists: true)
}

Channel
  .fromFilePairs (params.reads, size: params.bulk ? 1:2)
  .ifEmpty{ exit 1, "Cannot find any reads matching: ${params.reads}\nPath needs to be enclosed in quotes!\nPath requires at least one * wildcard!" }
  .into { raw_reads_fastqc; raw_reads_umi }


process rawFastqc {

  cpus { check_max( 8 * task.attempt, 'cpus' ) }
  memory { check_max( 5.GB * task.attempt, 'memory' ) }
  time { check_max( 6.h * task.attempt, 'time' ) }
  
  tag "FastQC on raw $library"
  publishDir "${params.outdir}/QC/raw_fastqc", mode:'copy', overwrite: true

  input:
  set name, file(reads) from raw_reads_fastqc

  output:
  file("*_fastqc.{zip,html}") into multiqc_fastqc

  script:
  library = name.toString().tokenize('_').get(0)
  """
  fastqc -t ${task.cpus} -f fastq -q ${reads} --nogroup
  """
}

if(params.extract_method){
  process extract_umi{

    cpus = { check_max( 2 * task.attempt, 'cpus' ) }
    memory = { check_max( 5.GB * task.attempt, 'memory' ) }
    time = { check_max( 6.h * task.attempt, 'time' ) }

    tag "extract UMI $library"
  
    input:
    set name, file(reads) from raw_reads_umi
  
    output:
    //set library, file("${library}_R1.fastq.gz"), file("${library}_R2.fastq.gz") into extracted_fastq
    set name, file("*.extracted.fastq.gz") into extracted_reads
  
    script:
    library = name.toString().tokenize('_').get(0)
    template "${params.extract_method}.sh"
  }
}else{
  raw_reads_umi.set{ extracted_reads }
}

if(params.trim_method){
  process trim {

    cpus = { check_max( 8 * task.attempt, 'cpus' ) }
    memory = { check_max( 25.GB * task.attempt, 'memory' ) }
    time = { check_max( 10.h * task.attempt, 'time' ) }
  
    tag "trimming $library"
  
    input:
    set name, file(reads) from extracted_reads
    // set library, file(read1), file(read2) from extracted_fastq
  
    output:
    set library, file("*.trimmed.fastq.gz") into trimmed_reads_alignment, trimmed_reads_contamination, trimmed_reads_fastqc
    //set library, file("${library}_R1_trimmed.fastq.gz"), file("${library}_R2_trimmed.fastq.gz") into trimmed_fastq
    file("${library}.trim_report.txt") into multiqc_cutadapt
  
    script:
    library = name.toString().tokenize('_').get(0)
    template "${params.trim_method}.sh"
  
  }
}else{
  extracted_reads.set{ trimmed_reads }
}

process trimmedFastqc {
  
  cpus = { check_max( 8 * task.attempt, 'cpus' ) }
  memory = { check_max( 5.GB * task.attempt, 'memory' ) }
  time = { check_max( 6.h * task.attempt, 'time' ) }
  
  tag "FastQC on trimmed $library"
  publishDir "${params.outdir}/QC/trimmed_fastqc", mode:'copy', overwrite: true

  input:
  set library, file(reads) from trimmed_reads_fastqc

  output:
  file("*_fastqc.{zip,html}") into multiqc_trimmed_fastqc

  script:
  """
  fastqc -t ${task.cpus} -f fastq -q ${reads} --nogroup
  """
}


// Contamination quantification
if(params.bulk){
  process STAR_contamination{
    cpus = { check_max( 16 * task.attempt, 'cpus' ) }
    memory = { check_max( 50.GB * task.attempt, 'memory' ) }
    time = { check_max( 5.h * task.attempt, 'time' ) }
    cache 'lenient'
    clusterOptions "--gres=tmpspace:50G"

    tag "STAR contamination alignment ${library}"

    publishDir "${params.outdir}/contamination/quantification", pattern: "*_ReadsPerGene.out.tab", mode:'copy', overwrite: true
    publishDir "${params.outdir}/contamination/alignment", pattern: "*.ba*", mode:'copy', overwrite: true

    input:
    file starContaminationIndex
    set library, file(reads) from trimmed_reads_contamination
    
    output:
    set file("*_Log.final.out"), file("*_ReadsPerGene.out.tab") into multiqc_star_contamination
    set library, file("${library}_contamination_Aligned.sortedByCoord.out.bam"), file("${library}_contamination_Aligned.sortedByCoord.out.bam.bai") into aligned_contamination
    file("*.bam")
    file("*.bai")

    script:
    """
    STAR \
      ${params.STARmapParams} \
      --outFilterMultimapNmax 100000000 \
      --outSAMmultNmax 1 \
      --outSAMunmapped Within \
      --limitBAMsortRAM 32212254720 \
      --runThreadN ${task.cpus} \
      --genomeDir ${starContaminationIndex} \
      --readFilesIn ${reads} \
      --outFileNamePrefix ${library}_contamination_ \
      --outSAMattributes ${params.STARoutSamAttributes}

    samtools index ${library}_contamination_Aligned.sortedByCoord.out.bam

    sambamba sort \
      -m 2G \
      -t ${task.cpus} \
      ${library}_contamination_Aligned.toTranscriptome.out.bam
    
    rm ${library}_contamination_Aligned.toTranscriptome.out.bam
    mv ${library}_contamination_Aligned.toTranscriptome.out.sorted.bam ${library}_contamination_Aligned.toTranscriptome.out.bam 
    mv ${library}_contamination_Aligned.toTranscriptome.out.sorted.bam.bai ${library}_contamination_Aligned.toTranscriptome.out.bam.bai
    
    """
  }
}else{
  process STARSolo_contamination{
    cpus = { check_max( 16 * task.attempt, 'cpus' ) }
    memory = { check_max( 50.GB * task.attempt, 'memory' ) }
    time = { check_max( 5.h * task.attempt, 'time' ) }
    cache 'lenient'
    clusterOptions "--gres=tmpspace:50G"

    tag "STARSolo contamination alignment ${library}"

    publishDir "${params.outdir}/contamination/quantification", pattern: "*_Solo.out", mode:'copy', overwrite: true
    publishDir "${params.outdir}/contamination/quantification", pattern: "_ReadsPerGene.out.tab", mode:'copy', overwrite: true
    publishDir "${params.outdir}/contamination/alignment", pattern: "*.ba*", mode:'copy', overwrite: true

    input:
    file whitelist
    file starContaminationIndex
    set library, file(reads) from trimmed_reads_contamination

    output:
    file("*_Solo.out")
    set file("*_Log.final.out"), file("*_ReadsPerGene.out.tab") into multiqc_star_contamination
    set library, file("${library}_contamination_Aligned.sortedByCoord.out.bam"), file("${library}_contamination_Aligned.sortedByCoord.out.bam.bai") into aligned_contamination
    file("*.bam")
    file("*.bai")

    script:
    if(params.cDNARead == "R2"){
      readOrder = "${reads[1]} ${reads[0]}"
    }else if(params.cDNARead == "R1"){
      readOrder = "${reads[0]} ${reads[1]}"
    }
    """
    STAR \
      ${params.STARmapParams} \
      --outFilterMultimapNmax 100000000 \
      --outSAMmultNmax 1 \
      --outSAMunmapped Within \
      --limitBAMsortRAM 32212254720 \
      --runThreadN ${task.cpus} \
      --genomeDir ${starContaminationIndex} \
      --readFilesIn ${readOrder} \
      --outFileNamePrefix ${library}_contamination_ \
      --outSAMattributes ${params.STARoutSamAttributes} ${params.STARsoloSamAttributes} \
      --soloType CB_UMI_Simple \
      --soloCBwhitelist ${whitelist} \
      --soloCBstart ${params.CBstart} \
      --soloCBlen ${params.CBlen} \
      --soloUMIstart ${params.UMIstart} \
      --soloUMIlen ${params.UMIlen} \
      --soloBarcodeReadLength 0 \
      --soloStrand ${params.readStrand} \
      --soloFeatures Gene GeneFull Velocyto \
      --soloUMIdedup 1MM_Directional \
      --soloCellFilter None

    removeEmptyCB.py \
      -b ${library}_contamination_Aligned.sortedByCoord.out.bam \
      -t ${task.cpus}

    rm ${library}_contamination_Aligned.sortedByCoord.out.bam
    mv ${library}_contamination_Aligned.sortedByCoord.out_stripped.bam ${library}_contamination_Aligned.sortedByCoord.out.bam

    samtools index ${library}_contamination_Aligned.sortedByCoord.out.bam
  
    addCB.py \
      -b ${library}_contamination_Aligned.toTranscriptome.out.bam \
      -t ${task.cpus} \
      -c ${whitelist} \
    
    sambamba sort \
      -m 2G \
      -t ${task.cpus} \
      ${library}_contamination_Aligned.toTranscriptome.out_CB.bam
    
    rm ${library}_contamination_Aligned.toTranscriptome.out.bam
    rm ${library}_contamination_Aligned.toTranscriptome.out_CB.bam
  
    mv ${library}_contamination_Aligned.toTranscriptome.out_CB.sorted.bam ${library}_contamination_Aligned.toTranscriptome.out.bam
    mv ${library}_contamination_Aligned.toTranscriptome.out_CB.sorted.bam.bai ${library}_contamination_Aligned.toTranscriptome.out.bam.bai
  
    find . -type f \\( -name '*.tsv' -o -name '*.mtx' \\) -exec gzip "{}" \\;
    """
  }
}

process count_contamination{
  cpus = { check_max( 2, 'cpus' ) }
  memory = { check_max( 10.GB * task.attempt, 'memory' ) }
  time = { check_max( 6.h * task.attempt, 'time' ) }

  tag "Quantifying contamination ${library}"

  publishDir "${params.outdir}/contamination/quantification", mode:'copy', overwrite: true

  input:
  set library, file(alignedBam), file(alignedBai) from aligned_contamination

  output:
  file("*.{csv,csv.gz}")

  script:
  """
  countContamination.py \
    --bam ${alignedBam} \
    --threads ${task.cpus}
  """
}


// Alignment 
if(params.bulk){
  process STAR{

    cpus = { check_max( 24 * task.attempt, 'cpus' ) }
    memory = { check_max( 50.GB * task.attempt, 'memory' ) }
    time = { check_max( 6.h * task.attempt, 'time' ) }
    cache 'lenient'
    clusterOptions "--gres=tmpspace:50G"

    tag "STAR Alignment ${library}"

    publishDir "${params.outdir}/quantification", pattern: "*_ReadsPerGene.out.tab", mode:'copy', overwrite: true
    publishDir "${params.outdir}/alignment", pattern: "*.ba*", mode:'copy', overwrite: true
  
    input:
    file starIndex
    set library, file(reads) from trimmed_reads_alignment 
  
    output:
    set file("*_Log.final.out"), file("*_ReadsPerGene.out.tab") into multiqc_star
    set library, file("${library}_Aligned.sortedByCoord.out.bam"), file("${library}_Aligned.sortedByCoord.out.bam.bai") into aligned_genome
    set library, file("${library}_Aligned.toTranscriptome.out.bam"), file("${library}_Aligned.toTranscriptome.out.bam.bai") into aligned_transcriptome, scRibo_transcriptome_all
      
    script:
    """
    STAR \
      ${params.STARmapParams} \
      --outFilterMultimapNmax 1 \
      --runThreadN ${task.cpus} \
      --genomeDir ${starIndex} \
      --readFilesIn ${reads} \
      --outFileNamePrefix ${library}_ \
      --outSAMattributes ${params.STARoutSamAttributes}
    
    samtools index ${library}_Aligned.sortedByCoord.out.bam
  
    sambamba sort \
      -m 2G \
      -t ${task.cpus} \
      ${library}_Aligned.toTranscriptome.out.bam
    
    rm ${library}_Aligned.toTranscriptome.out.bam
  
    mv ${library}_Aligned.toTranscriptome.out.sorted.bam ${library}_Aligned.toTranscriptome.out.bam
    mv ${library}_Aligned.toTranscriptome.out.sorted.bam.bai ${library}_Aligned.toTranscriptome.out.bam.bai
  
    """
  }
}else{
  process STARsolo {

    cpus = { check_max( 24 * task.attempt, 'cpus' ) }
    memory = { check_max( 50.GB * task.attempt, 'memory' ) }
    time = { check_max( 6.h * task.attempt, 'time' ) }
    cache 'lenient'
    clusterOptions "--gres=tmpspace:50G"
  
    tag "STARSolo Alignment $library"
  
    publishDir "${params.outdir}/alignment", pattern: "*.ba*", mode:'copy', overwrite: true
    publishDir "${params.outdir}/quantification", pattern: "*_Solo.out", mode:'copy', overwrite: true
    publishDir "${params.outdir}/quantification", pattern: "*_ReadsPerGene.out.tab", mode:'copy' , overwrite: true
  
    input:
    file whitelist
    file starIndex
    set library, file(reads) from trimmed_reads_alignment
    //set library, file(read1), file(read2) from trimmed_reads_alignment
  
    output:
    // file("*.bam")
    // file("*.bai")
    file("*_Solo.out")
    set file("*_Log.final.out"), file("*_ReadsPerGene.out.tab") into multiqc_star
    set library, file("${library}_Aligned.sortedByCoord.out.bam"), file("${library}_Aligned.sortedByCoord.out.bam.bai") into aligned_genome
    set library, file("${library}_Aligned.toTranscriptome.out.bam"), file("${library}_Aligned.toTranscriptome.out.bam.bai") into aligned_transcriptome, scRibo_transcriptome_all
      
    script:
    if(params.cDNARead == "R2"){
      readOrder = "${reads[1]} ${reads[0]}"
    }else if(params.cDNARead == "R1"){
      readOrder = "${reads[0]} ${reads[1]}"
    }
    """
    STAR \
      ${params.STARmapParams} \
      --outFilterMultimapNmax 1 \
      --runThreadN ${task.cpus} \
      --genomeDir ${starIndex} \
      --readFilesIn ${readOrder} \
      --outFileNamePrefix ${library}_ \
      --outSAMattributes ${params.STARoutSamAttributes} ${params.STARsoloSamAttributes} \
      --soloType CB_UMI_Simple \
      --soloCBwhitelist ${whitelist} \
      --soloCBstart ${params.CBstart} \
      --soloCBlen ${params.CBlen} \
      --soloUMIstart ${params.UMIstart} \
      --soloUMIlen ${params.UMIlen} \
      --soloBarcodeReadLength 0 \
      --soloStrand ${params.readStrand} \
      --soloFeatures Gene GeneFull Velocyto \
      --soloUMIdedup 1MM_Directional \
      --soloCellFilter None

    ## for STAR >= 2.7.8
    # removeEmptyCB.py \
    #   -b ${library}_Aligned.sortedByCoord.out.bam \
    #   -t ${task.cpus}

    # rm ${library}_Aligned.sortedByCoord.out.bam
    # mv ${library}_Aligned.sortedByCoord.out_stripped.bam ${library}_Aligned.sortedByCoord.out.bam

    samtools index ${library}_Aligned.sortedByCoord.out.bam
  
    addCB.py \
      -b ${library}_Aligned.toTranscriptome.out.bam \
      -t ${task.cpus} \
      -c ${whitelist} \
    
    sambamba sort \
      -m 2G \
      -t ${task.cpus} \
      ${library}_Aligned.toTranscriptome.out_CB.bam
    
    rm ${library}_Aligned.toTranscriptome.out.bam
    rm ${library}_Aligned.toTranscriptome.out_CB.bam
  
    mv ${library}_Aligned.toTranscriptome.out_CB.sorted.bam ${library}_Aligned.toTranscriptome.out.bam
    mv ${library}_Aligned.toTranscriptome.out_CB.sorted.bam.bai ${library}_Aligned.toTranscriptome.out.bam.bai
  
    find . -type f \\( -name '*.tsv' -o -name '*.mtx' \\) -exec gzip "{}" \\;
    """
  }
}

if(params.deduplicate_genome_method){
  process dedup_genome {
    memory { check_max( 50.GB + (100.GB*(task.attempt-1)), 'memory' ) }
    time { check_max( 40.h * task.attempt, 'time') }
    
    tag "umi_tools dedup on $library"
    
    input:
    set library, file(bam), file(bai) from aligned_genome
    
    output:
    set library, file("${bam.baseName}.dedup.bam") into sort_deduplicated_genome
    file("*_dedup.log") into multiqc_deduplicate_genome
    
    script:
    template "${params.deduplicate_genome_method}.sh"
  }
  
  process sort_dedup_genome {

    cpus = { check_max( 12 * task.attempt, 'cpus' ) }
    memory = { check_max( 10.GB * task.attempt, 'memory' ) }
    time = { check_max( 6.h * task.attempt, 'time' ) }
    clusterOptions "--gres=tmpspace:50G"
    
    tag "sorting deduped $library"
  
    publishDir "${params.outdir}/deduplicate", mode:'copy', overwrite: true
  
    input:
    set library, file(bam) from sort_deduplicated_genome
  
    output:
    file("${bam.baseName}.sorted.bam")
    file("${bam.baseName}.sorted.bam.bai")
    
    script:
    """
    sambamba sort \
      -m 4G \
      -t ${task.cpus} \
      ${bam}
  
    """
  }
}else{
  // aligned_genome .set { scRibo_transcriptome }
  multiqc_deduplicate_genome = Channel.empty()
}


if(params.deduplicate_transcriptome_method){
  process dedup_transcriptome {
    memory { check_max( 50.GB + (100.GB*(task.attempt-1)), 'memory' ) }
    time { check_max( 40.h * task.attempt, 'time') }
    
    tag "umi_tools dedup on $library"
    
    input:
    set library, file(bam), file(bai) from aligned_transcriptome
    
    output:
    set library, file("${bam.baseName}.dedup.bam") into sort_deduplicated_transcriptome
    file("*_dedup.log") into multiqc_deduplicate_transcriptome
    
    script:
    template "${params.deduplicate_transcriptome_method}.sh"
  }
  
  process sort_dedup_transcriptome {
    cpus = { check_max( 12 * task.attempt, 'cpus' ) }
    memory = { check_max( 10.GB * task.attempt, 'memory' ) }
    time = { check_max( 6.h * task.attempt, 'time' ) }
    clusterOptions "--gres=tmpspace:50G"
    
    tag "sorting deduped $library"
  
    publishDir "${params.outdir}/deduplicate", mode:'copy', overwrite: true
  
    input:
    set library, file(bam) from sort_deduplicated_transcriptome
  
    output:
    // file("${bam.baseName}.sorted.bam")
    // file("${bam.baseName}.sorted.bam.bai")
    set library, file("${bam.baseName}.sorted.bam"), file("${bam.baseName}.sorted.bam.bai") into scRibo_transcriptome 
    
    script:
    """
    sambamba sort \
      -m 4G \
      -t ${task.cpus} \
      ${bam}
  
    """
  }
}else{
  aligned_transcriptome .set { scRibo_transcriptome }
  multiqc_deduplicate_transcriptome = Channel.empty()
  scRibo_transcriptome_all = Channel.empty()
  scRibo_count_all = Channel.empty()
  scRibo_biotypes_all = Channel.empty()
}


scRibo_parse_all = (params.parseAllReads
                 ? scRibo_transcriptome_all
                 : Channel.empty() )


process scRibo_parseAllReads {
  cpus = { check_max( 2 * task.attempt, 'cpus' ) }
  memory = { check_max( 75.GB + (200.GB*(task.attempt-1)), 'memory' ) }
  time = { check_max( 8.h * task.attempt, 'time' ) }

  tag "scRibo ${library}"

  publishDir "${params.outdir}/scRibo", mode:'copy', overwrite: true

  input:
  set library, file(bam), file(bai) from scRibo_parse_all
  file(annotationPickle)
  file(genomeFasta)
  file(genomeFai)

  output:
  file("${bam.baseName}.{csv.gz,feather}")
  set library, file("${bam.baseName}.{csv.gz,feather}") into count_all, biotypes_all
 
  script:
  """
  parseReads.py \
    -t ${annotationPickle} \
    -f ${genomeFasta} \
    -b ${bam} \
    ${params.parseReads}
  """
}

scRibo_count_all = (params.parseAllReads
                 ? count_all
                 : Channel.empty() )

scRibo_biotypes_all = (params.parseAllReads
                 ? biotypes_all
                 : Channel.empty() )

process scRibo_count_all {
  memory { check_max( 50.GB * task.attempt, 'memory') }
  cpus = { check_max( 1 * task.attempt, 'cpus' ) }
  time = { check_max( 4.h * task.attempt, 'time' ) }
  clusterOptions '--gres=tmpspace:50G'

  tag "Counting raw alignments on ${inputAlignments}"
  publishDir "${params.outdir}/scRibo/all/quantification", mode:'copy', overwrite: true

  input:
  file(annotationFeather)
  set library, file(inputAlignments) from scRibo_count_all

  output:
  file("${library}")

  script:
  count_site = 'cut5'
  count_expansion = params.count_expansion
  template "${params.count_method}.sh"
}

process scRibo_biotypes_all {
  memory { check_max( 50.GB * task.attempt, 'memory') }
  cpus = { check_max( 1 * task.attempt, 'cpus' ) }
  time = { check_max( 4.h * task.attempt, 'time' ) }
  clusterOptions '--gres=tmpspace:50G'

  tag "Counting raw alignments on ${inputAlignments}"
  publishDir "${params.outdir}/scRibo/all/quantification/biotypes", mode:'copy', overwrite: true

  input:
  file(annotationFeather)
  set library, file(inputAlignments) from scRibo_biotypes_all

  output:
  file("${library}")

  script:
  template "${params.biotypes_method}.sh"
}


process scRibo_parsereads {
  cpus = { check_max( 2 * task.attempt, 'cpus' ) }
  memory = { check_max( 75.GB + (200.GB*(task.attempt-1)), 'memory' ) }
  time = { check_max( 8.h * task.attempt, 'time' ) }

  tag "scRibo ${library}"

  publishDir "${params.outdir}/scRibo", mode:'copy', overwrite: true

  input:
  set library, file(bam), file(bai) from scRibo_transcriptome
  file(annotationPickle)
  file(genomeFasta)
  file(genomeFai)

  output:
  // file("${bam.baseName}.{csv.gz,feather}")
  set library, file("${bam.baseName}.{csv.gz,feather}") into parsed_alignments_qc, parsed_alignments_predict
 
  script:
  """
  parseReads.py \
    -t ${annotationPickle} \
    -f ${genomeFasta} \
    -b ${bam} \
    ${params.parseReads}
  """
}

if(params.bulk){
  process scRibo_QCplots_bulk {
    memory { check_max( 75.GB + (100.GB*(task.attempt-1)), 'memory') }
    cpus = { check_max( 2 * task.attempt, 'cpus' ) }
    time = { check_max( 8.h * task.attempt, 'time' ) }
    clusterOptions '--gres=tmpspace:50G' // for data.table unzip
  
    tag "scRibo ${library}"
  
    publishDir "${params.outdir}/scRibo", mode:'copy', overwrite: true
  
    input:
    set library, file(readFeather) from parsed_alignments_qc
    file(annotationFeather)
  
    output:
    file("*.pdf")
  
    script:
    template "${params.qc_plots_method}.sh"
   
  }
}else{
  process scRibo_QCplots {
    memory { check_max( 75.GB + (100.GB*(task.attempt-1)), 'memory') }
    cpus = { check_max( 2 * task.attempt, 'cpus' ) }
    time = { check_max( 8.h * task.attempt, 'time' ) }
    clusterOptions '--gres=tmpspace:50G' // for data.table unzip
  
    tag "scRibo ${library}"
  
    publishDir "${params.outdir}/scRibo", mode:'copy', overwrite: true
  
    input:
    set library, file(readFeather) from parsed_alignments_qc
    file(annotationFeather)
    file whitelist from whitelist
  
    output:
    file("*.pdf")
  
    script:
    template "${params.qc_plots_method}.sh"
   
  }
}


if(params.ranger_model){
  rangerModel = file(params.ranger_model, checkIfExists: true)
  
  count_expansion = 10
  count_site = "psite"
  // params.count_expansion = 10
  
  process applyRandomForest {
    memory { check_max( 150.GB + (50.GB*(task.attempt-1)), 'memory') }
    cpus = { check_max( 24, 'cpus' ) }
    time = { check_max( 10.h * task.attempt, 'time' ) }
    clusterOptions '--gres=tmpspace:15G'
    
    tag "Applying Predictions on ${inputAlignments}"
    publishDir "${params.outdir}/scRibo/predicted", mode:'copy', overwrite: true
  
    input:
    rangerModel
    set library, file(inputAlignments) from parsed_alignments_predict
  
    output:
    // file("*.predicted.{csv.gz, feather}")
    set library, file("*.predicted.{csv.gz, feather}") into alignments_count, alignments_biotypes

    script:
    """
    predictModel.R \
      -m ${rangerModel} \
      -r ${inputAlignments}
    """
  }
}else{
  parsed_alignments_predict .into { alignments_count; alignments_biotypes }
  count_expansion = params.count_expansion
  count_site = params.count_site
}

process scRibo_count {
  memory { check_max( 50.GB * task.attempt, 'memory') }
  cpus = { check_max( 1 * task.attempt, 'cpus' ) }
  time = { check_max( 4.h * task.attempt, 'time' ) }
  clusterOptions '--gres=tmpspace:50G'

  tag "Counting alignments on ${inputAlignments}"
  publishDir "${params.outdir}/scRibo/quantification", mode:'copy', overwrite: true

  input:
  file(annotationFeather)
  set library, file(inputAlignments) from alignments_count
  val count_expansion from count_expansion
  val count_site

  output:
  file("${library}")

  script:
  template "${params.count_method}.sh"
}

process scRibo_biotypes {
  memory { check_max( 50.GB * task.attempt, 'memory') }
  cpus = { check_max( 1 * task.attempt, 'cpus' ) }
  time = { check_max( 4.h * task.attempt, 'time' ) }
  clusterOptions '--gres=tmpspace:50G'

  tag "Counting alignments on ${inputAlignments}"
  publishDir "${params.outdir}/scRibo/quantification/biotypes", mode:'copy', overwrite: true

  input:
  file(annotationFeather)
  set library, file(inputAlignments) from alignments_biotypes
  val count_expansion from count_expansion
  val count_site

  output:
  file("${library}")

  script:
  template "${params.biotypes_method}.sh"
}




process multiQC {
  memory { check_max( 15.GB * task.attempt, 'memory') }
  cpus = { check_max( 1 , 'cpus' ) }
  time = { check_max( 4.h * task.attempt, 'time' ) }

  publishDir "${params.outdir}/QC", mode:'copy', overwrite: true

  input:
  file('*') from multiqc_fastqc.collect().ifEmpty([])
  file('*') from multiqc_trimmed_fastqc.collect().ifEmpty([])
  file('*') from multiqc_cutadapt.collect().ifEmpty([])
  file('*') from multiqc_star.collect().ifEmpty([])
  file('*') from multiqc_deduplicate_genome.collect().ifEmpty([])
  file('*') from multiqc_deduplicate_transcriptome.collect().ifEmpty([])
  file('*') from multiqc_star_contamination.collect().ifEmpty([])

  output:
  file "multiqc_report.html"
  file "multiqc_data"

  script:
  """
  multiqc .
  """
}

process store_condaenv {
  
  publishDir "${params.tracedir}", pattern:"environment.yml", mode:'copy', overwrite:true

  input:
  genomeFasta
  
  output:
  file "environment.yml"

  script:
  """
  conda env export > environment.yml
  """
}
