RiboQC.R \
  --reads ${readFeather} \
  --annotations ${annotationFeather} \
  --whitelist ${whitelist} \
  --canonical

rm -f Rplots.pdf
