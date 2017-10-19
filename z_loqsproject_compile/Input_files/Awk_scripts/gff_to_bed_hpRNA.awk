BEGIN{OFS="\t"}/^#/ {next} {print $1, $4, $5, $9,".", $7}END{} 

