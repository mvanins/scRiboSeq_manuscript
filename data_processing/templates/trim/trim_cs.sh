cutadapt \
	--cores=${task.cpus} \
	-U 6 \
	-U -1 \
	--times 5 \
	-m :${params.minlength} \
	-A ${params.threeadapter} \
	-A AAAAAAAAAAAAAAA \
	-A TTTTTTTTTTTTTTT \
	-A CCCCCCCCCCCCCCC \
	-A GGGGGGGGGGGGGGG \
	-o ${library}_R1.trimmed.fastq.gz \
	-p ${library}_R2.trimmed.fastq.gz \
	${reads[0]} ${reads[1]} > ${library}.trim_report.txt
