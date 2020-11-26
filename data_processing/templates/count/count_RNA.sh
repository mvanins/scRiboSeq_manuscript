countTables.R \
		--reads ${inputAlignments} \
		-a ${annotationFeather} \
		--rna \
		--site ${count_site} \
		--minlength 0 \
		--maxlength 1000

