# miRv6 move UMI to be compatable with STARSolo
moveUMI.pl \
	--length ${params.UMIlen} \
	--read1 ${reads[0]} \
	--read2 ${reads[1]} \
	--out1 ${library}_R1.extracted.fastq.gz \
	--out2 ${library}_R2.extracted.fastq.gz
