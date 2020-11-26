countTables.R \
		--reads ${inputAlignments} \
		--annotations ${annotationFeather} \
		--site ${count_site} \
		--expansion ${count_expansion} \
		--minlength 0 \
		--maxlength 200
