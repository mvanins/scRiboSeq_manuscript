cutadapt \
  --cores=${task.cpus} \
  -m ${params.minlength} \
  -a ${params.threeadapter} \
  -a AAAAAAAAAAAAAAA \
  -o ${library}_R1.trimmed.fastq.gz \
  ${reads[0]} > ${library}.trim_report.txt
