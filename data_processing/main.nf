#!/usr/bin/env nextflow

/*
Ribosome footprint analysis pipeline
*/


if(!params.genomeDir) {
	// build geneome directory if not supplied
	
	referenceOut = "${params.outdir}/reference"
	annotations = file(params.annotations, checkIfExists: true)

	genomeFasta = file(params.genomeFasta, checkIfExists: true)

	process index_genome {
		cpus 1
		time '2h'
		memory '10G'

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
	
		cpus 36
		time '24h'
		memory '50G'
		
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
	
		tRNAscan2bed12.pl ${genome.baseName}.Eukaryotic_tRNAs.out ${genome.baseName}.Eukaryotic_tRNAs_bp.bed
		sed '/^MT\\|^chrM/d' ${genome.baseName}.Eukaryotic_tRNAs_bp.bed > temp.bed
		mv temp.bed ${genome.baseName}.Eukaryotic_tRNAs_bp.bed
		"""
	}
	
	process mitochondrial_tRNAscan {
		// follows the recommended paramaters for the updated tRNAscan-SE 
		// from https://dx.doi.org/10.1007%2F978-1-4939-9173-0_1
	
		cpus 36
		time '24h'
		memory '50G'
		
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
	
		tRNAscan2bed12.pl ${genome.baseName}.Mitochondrial_tRNAs.out ${genome.baseName}.Mitochondrial_tRNAs_bp.bed
		"""
	}
	
	process create_genomes{
	
		time '2h'
		memory '50G'
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
		file("${genome.baseName}.tRNACluster.gtf")
		file("${genome.baseName}.tRNA-masked-withmature.fa")
		file("${genome.baseName}.tRNA-masked-withmature.fa.fai")
		// miR file("${annotations.baseName}.miRNA.tRNA.gtf")
		file("${annotations.baseName}.tRNA.gtf")
		set file("${genome.baseName}.tRNA-masked-withmature.fa"), file("${annotations.baseName}.tRNA.gtf") into star_genome
	
		script:
		"""
		#### tRNAs
		# mask tRNAs in genome
		cat ${eukaryotic} ${mitochondrial} | sort -k 1,1 -k2,2n  > ${genome.baseName}.tRNA_masked_regions.bed
		bedtools maskfasta -fi ${genome} -fo ${genome.baseName}.tRNA_masked.fa -bed ${genome.baseName}.tRNA_masked_regions.bed
	
		## pre tRNAs
		# remove pseudogenes, add 50 nt 5' and 3' flanking regions
		modBed12.pl ${eukaryotic} ${eukaryotic.baseName}.pre-tRNAs.bed
	
		# select only MT tRNAs from the MT model, add flaking regions
		sed '/^MT\\|^chrM/!d' ${mitochondrial} > ${mitochondrial.baseName}.only_MT.bed
		modBed12.pl ${mitochondrial.baseName}.only_MT.bed ${mitochondrial.baseName}.pre-tRNAs.bed
	
		# combine to pre-tRNAs.bed
		cat ${eukaryotic.baseName}.pre-tRNAs.bed ${mitochondrial.baseName}.pre-tRNAs.bed | sort -k 1,1 -k2,2n > ${genome.baseName}.pre-tRNAs.bed
	
		# extract pre-tRNA sequences
		bedtools getfasta -name -split -s -fi ${genome} -bed ${genome.baseName}.pre-tRNAs.bed -fo ${genome.baseName}.pre-tRNAs.fa
		
		## from "best practices workflow", we will only align to clustered mature tRNAs
		##  # create pre-tRNAs and masked genome
		##  cat ${genome.baseName}.tRNA_masked.fa ${genome.baseName}.pre-tRNAs.fa > ${genome.baseName}.artificial.fa
		##  samtools faidx ${genome.baseName}.artificial.fa
	
		## mature tRNAs
		cat ${eukaryotic} ${mitochondrial.baseName}.only_MT.bed | sort -k 1,1 -k2,2n > ${genome.baseName}.mature-tRNAs.bed
		bedtools getfasta -name -split -s -fi ${genome} -bed ${genome.baseName}.mature-tRNAs.bed -fo temp_mature-tRNAs.fa
		addCCA.pl temp_mature-tRNAs.fa ${genome.baseName}.mature-tRNAs.fa 
	
		## cluster mature tRNAs
		collapseSequences.pl ${genome.baseName}.mature-tRNAs.fa > ${genome.baseName}.mature-tRNAs.cluster.fa
		# produces cluster_info.gtf
		mv cluster_info.gtf ${genome.baseName}.tRNACluster.gtf
	
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
		time '5h'
		memory '75G'

		tag "scRibo annotation preprocessing for $genome"
		
		beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/venv38/bin/activate"
		afterScript "deactivate"

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
	
		cpus 48
		time '12h'
		memory '100G'
	
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
	
}else{
	// genome exists, set variables
	annotationFeather = file("${params.genomeDir}/annotations/*.annotations.{csv.gz,feather}", checkIfExists: true).first()
	annotationPickle = file("${params.genomeDir}/annotations/*.annotations.pickle", checkIfExists: true).first()
	genomeFasta = file("${params.genomeDir}/original/*.{fa,fasta}", checkIfExists: true).first()
	genomeFai =  file("${params.genomeDir}/original/*.fai", checkIfExists: true).first()
	starIndex = file("${params.genomeDir}/star_tRNAmasked_${params.starOverhang}", checkIfExists: true)
}


if(!params.bulk){
	whitelist = file(params.whitelist, checkIfExists: true)
}

if(params.deplete){
	starDepleteIndex = file(params.starDepleteIndex, checkIfExists: true)
}


Channel
	.fromFilePairs (params.reads, size: params.bulk ? 1:2)
	.ifEmpty{ exit 1, "Cannot find any reads matching: ${params.reads}\nPath needs to be enclosed in quotes!\nPath requires at least one * wildcard!" }
	.into { raw_reads_fastqc; raw_reads_umi }


process rawFastqc {
	cpus 8
	time '5h'
	memory '5G'
	
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
		time '8h'
		memory '10G'
		cpus 2

                beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/venv38/bin/activate"
                afterScript "deactivate"
	
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
	process cutadapt {
		time '5h'
		cpus 8
		memory '10G'
	
		tag "cutadapt $library"
	
		beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/venv36/bin/activate"
		afterScript "deactivate"
	
		input:
		set name, file(reads) from extracted_reads
		// set library, file(read1), file(read2) from extracted_fastq
	
		output:
		set library, file("*.trimmed.fastq.gz") into trimmed_reads, trimmed_reads_fastqc
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
	cpus 8
	time '5h'
	memory '5G'
	
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

if (!params.deplete){
	trimmed_reads .set { depleted_reads }
	multiqc_star_depletion = Channel.empty()
	
}else{
	process STAR_deplete {
		cpus 24
		time '4h'
		memory '50G'
		cache 'lenient'
		clusterOptions "--gres=tmpspace:50G"

		tag "depleting contaminant sequences on $library"

		input:
		file starDepleteIndex from starDepleteIndex
		file whitelist from whitelist
		// set library, file(read1), file(read2) from trimmed_fastq
		set library, file(reads) from trimmed_fastq

		output:
		set library, file("*.depleted.fastq.gz") into depleted_reads
		set file("*_Log.final.out"), file("*_ReadsPerGene.out.tab") into multiqc_star_depletion
		
		script:
		"""
		STAR \
			--runThreadN ${task.cpus} \
			--genomeDir ${starDepleteIndex} \
			--readFilesIn ${read1} ${read2} \
			--readFilesCommand zcat \
			--outFileNamePrefix ${library}_depleted_ \
			--outSAMtype BAM SortedByCoordinate \
			--outSAMattributes NH HI AS nM NM MD jM jI XS MC ch CR UR CY UY GX GN sS sQ CB UB sM \
			--quantMode TranscriptomeSAM GeneCounts \
			--quantTranscriptomeBan Singleend \
			--seedSearchStartLmax 10 \
			--alignIntronMax 1000000 \
			--outFilterScoreMin 0 \
			--outFilterMultimapNmax 10000000 \
			--chimScoreSeparation 10 \
			--chimScoreMin 20 \
			--chimSegmentMin 15 \
			--outFilterMismatchNmax 5 \
			--soloType CB_UMI_Simple \
			--soloCBwhitelist ${whitelist} \
			--soloCBstart ${params.CBstart} \
			--soloCBlen ${params.CBlen} \
			--soloUMIstart ${params.UMIstart} \
			--soloUMIlen ${params.UMIlen} \
			--soloBarcodeReadLength 1 \
 			--soloStrand Forward \
 			--soloFeatures Gene GeneFull SJ Velocyto \
 			--soloUMIdedup 1MM_Directional \
 			--soloCellFilter None \
			--outReadsUnmapped Fastx \
			--limitBAMsortRAM 1020159082
		
		mv ${library}_depleted_Unmapped.out.mate1 ${library}_R1.depleted.fastq
		mv ${library}_depleted_Unmapped.out.mate2 ${library}_R2.depleted.fastq

		pigz *.fastq
		"""
	}
}

// Alignment 
if(params.bulk){
	process STAR{
		cpus 24
		time '4h'
		memory '75G'
		cache 'lenient'
		clusterOptions "--gres=tmpspace:50G"

		tag "STAR Alignment ${library}"

		publishDir "${params.outdir}/quantification", pattern: "*_ReadsPerGene.out.tab", mode:'copy', overwrite: true
		publishDir "${params.outdir}/alignment", pattern: "*.ba*", mode:'copy', overwrite: true
	
		input:
		file starIndex
		set library, file(reads) from depleted_reads 
	
		output:
		set file("*_Log.final.out"), file("*_ReadsPerGene.out.tab") into multiqc_star
		set library, file("${library}_Aligned.sortedByCoord.out.bam"), file("${library}_Aligned.sortedByCoord.out.bam.bai") into aligned_genome
		set library, file("${library}_Aligned.toTranscriptome.out.bam"), file("${library}_Aligned.toTranscriptome.out.bam.bai") into aligned_transcriptome
			
		script:
		"""
		STAR \
			${params.STARmapParams} \
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
		cpus 24
		time '4h'
		memory '85G'
		cache 'lenient'
		clusterOptions "--gres=tmpspace:50G"
	
		tag "STARSolo Alignment $library"
	
		beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/venv38/bin/activate"
		afterScript "deactivate"
	
		publishDir "${params.outdir}/alignment", pattern: "*.ba*", mode:'copy', overwrite: true
		publishDir "${params.outdir}/quantification", pattern: "*_Solo.out", mode:'copy', overwrite: true
		publishDir "${params.outdir}/quantification", pattern: "*_ReadsPerGene.out.tab", mode:'copy' , overwrite: true
	
		input:
		file whitelist
		file starIndex
		set library, file(reads) from depleted_reads
		//set library, file(read1), file(read2) from trimmed_reads_alignment
	
		output:
		// file("*.bam")
		// file("*.bai")
		file("*_Solo.out")
		set file("*_Log.final.out"), file("*_ReadsPerGene.out.tab") into multiqc_star
		set library, file("${library}_Aligned.sortedByCoord.out.bam"), file("${library}_Aligned.sortedByCoord.out.bam.bai") into aligned_genome
		set library, file("${library}_Aligned.toTranscriptome.out.bam"), file("${library}_Aligned.toTranscriptome.out.bam.bai") into aligned_transcriptome
			
		script:
		if(params.readStrand == "Reverse"){
			readOrder = "${reads[1]} ${reads[0]}"
		}else{
			readOrder = "${reads[0]} ${reads[1]}"
		}
		"""
		STAR \
			${params.STARmapParams} \
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
	 		--soloFeatures Gene GeneFull SJ Velocyto \
	 		--soloUMIdedup 1MM_Directional \
	 		--soloCellFilter None
		
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
		memory '200G'
		time '50h'
		// clusterOptions "--gres=tmpspace:50G"
		
		tag "umi_tools dedup on $library"
		
		beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/venv38/bin/activate"
		afterScript "deactivate"
		
		input:
		set library, file(bam), file(bai) from aligned_genome
		
		output:
		set library, file("${bam.baseName}.dedup.bam") into sort_deduplicated_genome
		file("*_dedup.log") into multiqc_deduplicate_genome
		
		script:
		template "${params.deduplicate_genome_method}.sh"
	}
	
	process sort_dedup_genome {
		memory '50G'
		cpus 12
		time '2h'
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
		memory '200G'
		time '50h'
		// clusterOptions "--gres=tmpspace:50G"
		
		tag "umi_tools dedup on $library"
		
		beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/venv38/bin/activate"
		afterScript "deactivate"
		
		input:
		set library, file(bam), file(bai) from aligned_transcriptome
		
		output:
		set library, file("${bam.baseName}.dedup.bam") into sort_deduplicated_transcriptome
		file("*_dedup.log") into multiqc_deduplicate_transcriptome
		
		script:
		template "${params.deduplicate_transcriptome_method}.sh"
	}
	
	process sort_dedup_transcriptome {
		memory '50G'
		cpus 12
		time '2h'
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
}
	

process scRibo_parsereads {
	cpus 2
	time '20h'
	memory { 75.GB + (200.GB*(task.attempt-1)) }
	// memory {200.GB + 50.GB * task.attempt}
	errorStrategy 'retry'

	tag "scRibo ${library}"

	publishDir "${params.outdir}/scRibo", mode:'copy', overwrite: true

	beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/venv38/bin/activate"
	afterScript "deactivate"

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

process scRibo_QCplots {
	memory { 75.GB + (100.GB*(task.attempt-1)) }
	// memory {200.GB + 50.GB * task.attempt}
	errorStrategy 'retry'
	cpus 2
	time '2h'
	clusterOptions '--gres=tmpspace:50G' // for data.table unzip

	tag "scRibo ${library}"

	publishDir "${params.outdir}/scRibo", mode:'copy', overwrite: true

	beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/venv38/bin/activate"
	afterScript "deactivate"

	input:
	set library, file(readFeather) from parsed_alignments_qc
	file(annotationFeather)
	file whitelist from whitelist

	output:
	file("*.pdf")
 
	script:
	"""
	RiboQC.R \
		-r ${readFeather} \
		-a ${annotationFeather} \
		-w ${whitelist}

	rm -f Rplots.pdf
	"""
}

if(params.ranger_model){
	rangerModel = file(params.ranger_model, checkIfExists: true)
	
	count_expansion = 10
	count_site = "psite"
	// params.count_expansion = 10
	
	process applyRandomForest {
		cpus 24
		memory '350G'
		clusterOptions '--gres=tmpspace:15G'
		time '8h'
	
		
		tag "Applying Predictions on ${inputAlignments}"
		publishDir "${params.outdir}/scRibo/predicted", mode:'copy', overwrite: true
	
		input:
		rangerModel
		set library, file(inputAlignments) from parsed_alignments_predict
	
		output:
		// file("*.predicted.{csv.gz, feather}")
		set library, file("*.predicted.{csv.gz, feather}") into alignments_count

		script:
		"""
		predictModel.R \
			-m ${rangerModel} \
			-r ${inputAlignments}
		"""
	}
}else{
	parsed_alignments_predict .set { alignments_count }
	count_expansion = params.count_expansion
	count_site = params.count_site
}

process scRibo_count {
	cpus 1
	memory '50G'
	clusterOptions '--gres=tmpspace:15G'
	time '4h'

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


process multiQC {
	time '1h'
	memory '5G'

	publishDir "${params.outdir}/QC", mode:'copy', overwrite: true

	beforeScript "source /hpc/hub_oudenaarden/mvanins/local/virtualEnvironments/mqcdev/bin/activate"
	afterScript "deactivate"
	
	input:
	file('*') from multiqc_fastqc.collect().ifEmpty([])
	file('*') from multiqc_trimmed_fastqc.collect().ifEmpty([])
	file('*') from multiqc_cutadapt.collect().ifEmpty([])
	file('*') from multiqc_star.collect().ifEmpty([])
	file('*') from multiqc_deduplicate_genome.collect().ifEmpty([])
	file('*') from multiqc_deduplicate_transcriptome.collect().ifEmpty([])
	file('*') from multiqc_star_depletion.collect().ifEmpty([])

	output:
	file "multiqc_report.html"
	file "multiqc_data"

	script:
	"""
	multiqc .
	"""
}

