BEGIN{OFS="\t"}/^#/ {next} {$3=="miRNA_primary_transcript";print "chr"$1, $4, $5, $9,".", $7}END{} 

