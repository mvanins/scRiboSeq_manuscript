cutadapt \
	--cores=${task.cpus} \
	-m ${params.minlength}: \
	-a ${params.threeadapter} \
	-o ${library}_R1.trimmed.fastq.gz \
	-p ${library}_R2.trimmed.fastq.gz \
	${reads[0]} ${reads[1]} > ${library}.trim_report.txt
