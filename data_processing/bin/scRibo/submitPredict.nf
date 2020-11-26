#!/usr/bin/env nextflow

rangerModel = file(params.rangerModel, checkIfExists: true)

Channel
	.fromPath("${params.feather}")
	.ifEmpty { error "Cannot find any files"}
	.set { input_alignments }

process predictReads {
	cpus 24
	memory '350G'
	time '8h'

	
	tag "Applying Predictions on ${inputAlignments}"
	publishDir "${params.outdir}", mode:'copy', overwrite: true

	input:
	rangerModel
	file(inputAlignments) from input_alignments

	output:
	file("*.predicted.{csv.gz,feather}")

	script:
	"""
	Rscript /hpc/hub_oudenaarden/mvanins/local/avopipelines/RPF/bin/scRibo/predictModel.R \
		-m ${rangerModel} \
		-r ${inputAlignments}
	"""
}
