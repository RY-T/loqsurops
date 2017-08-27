#helppage
#http://slowkow.com/notes/snakemake-tutorial

#Index3.fa obtained from NCBI M21017.1
 
CHR_LIST= [line.rstrip('\n') for line in open('Chr_list.txt')]
RM_GENETYPES=['LTR','LINE','DNA','Satellite','Low_complexity','RC','Simple_repeat','Other','Unknown','ARTEFACT','rRNA']

#requires installation
RepeatMasker_Dir='/usr/local/RepeatMasker/RepeatMasker'

rule all:
    input:
        'Index0.fa',
        'RM_output',
        'RM_output.bed'

rule extract_impt_chr:
    input:
        genome = 'dmel-all-chromosome-r6.16.fasta'
    output:
        '{Chr_list}.temp'
    shell:"""
    samtools faidx {input.genome}
    samtools faidx {input.genome} {wildcards.Chr_list} > {output}
    """
        
rule make_Index_0:
	input:
		expand('{Chr_list}.temp', Chr_list=CHR_LIST)
	output:
		'Index0.fa'
	shell:"""
	cat {input} > {output}
	rm {input}
	"""
rule repeatMasker:
	input: 
		genome = 'Index0.fa',
		RM_Dir = RepeatMasker_Dir
	output:
		'RM_output'
	shell:"""
	mkdir -p {output} && {input.RM_Dir} -dir {output} -species drosophila {input.genome}
	"""

rule extract_bed_from_RMout: 
#makes modified bedfiles chr start stop 
	input:
		'RM_output/Index0.fa.out' 
	output:
		'RM_output.bed'
	shell: 'awk -f ./RMout2bed.awk {input} > {output}'


#rule modified_bed_to_fasta
#	mkdir -p {output.out_dir}
#	awk '$11~ /{wildcards.genetype}/ {print $5 "\t" $6-1 "\t" $7-1 "\t"\
#	"{wildcards.genetype}\:"$10 "\t" "." "\t" $9}'{input.outfile}\
#	> {wildcards.genetype}dm6.bed