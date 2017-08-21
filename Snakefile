#helppage
#http://slowkow.com/notes/snakemake-tutorial

CHR_LIST= [line.rstrip('\n') for line in open('Chr_list.txt')]

RepeatMasker_Dir='/usr/local/RepeatMasker/RepeatMasker'


rule all:
    input:
        'Index0.fa',
        'RM_output'

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


#rule collate_outputs:
#    input:
#    	'RM_output'
#        
#    output:
#        'test.txt'
#    run:
#        with open(output[0], 'w') as out:
#            for i in input:
#                sample = i.split('.')[0]
#                for line in open(i):
#                    out.write(sample + ' ' + line)