umi_tools dedup \
	-I ${bam} \
	-S ${bam.baseName}.dedup.bam \
	-L ${bam.baseName}_dedup.log \
	--spliced-is-unique \
	--per-cell \
	--read-length \
	--no-sort-output
