#helppage
#http://slowkow.com/notes/snakemake-tutorial

#Index3.fa obtained from NCBI M21017.1
 
CHR_LIST= [line.rstrip('\n') for line in open('Chr_list.txt')]
RM_GENETYPES=['LTR','LINE','DNA','Satellite','Low_complexity','RC','Simple_repeat','Other','Unknown','ARTEFACT','RNA','rRNA','tRNA']
Flybase_GENETYPES2=['rRNA', 'tRNA', 'snRNA', 'snoRNA']
#requires installation
RepeatMasker_Dir='/usr/local/RepeatMasker/RepeatMasker'

import pandas as pd

rule all:
    input:
        'Index0.fa',
        'RM_output',
        'modified_bed/RM_output.bed',
#        'modified_bed/RNA.bed'
        expand('modified_bed/{RM_Genetype}RM.bed', RM_Genetype=RM_GENETYPES),
        'modified_bed/RNAonlyRM.bed'

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
#awk command obtained from https://www.biostars.org/p/128068/
	input:
		'RM_output/Index0.fa.out' 
	output:
		'modified_bed/RM_output.bed'
	shell: """
	mkdir -p 'modified_bed'
	awk -f ./RMout2bed.awk {input} > {output}
	"""

#rule pandas_search_from_modified_bed_part1:
#	input:'modified_bed/RM_output.bed'
#part 1 extract RNA out first, part2 then from RNA, extract tRNA, etc...
#	output: 
#		expand('modified_bed/{RM_Genetype}.bed', RM_Genetype=RM_GENETYPES)
#	run:"""
#		everything = pd.read_csv('{input}',sep="\t+",header=None)
#		something = everything[everything.iloc[:,3].str.contains({wildcards.RM_Genetype})]
#		something.to_csv('{output}',sep='\t',index=False)
#		"""

rule pandas_search_from_modified_bed_part1:
	input:
		bedfile = 'modified_bed/RM_output.bed'
#part 1 extract RNA out first, part2 then from RNA, extract tRNA, etc...
	output: expand('modified_bed/{RM_Genetype}RM.bed', RM_Genetype=RM_GENETYPES)
	run:
		import pandas as pd
		import os

		everything = pd.read_csv(input.bedfile ,sep="\t+",header=None)
		RM_Genetype	=['LTR','LINE','DNA','Satellite','Low_complexity','RC','Simple_repeat','Other','Unknown','ARTEFACT','RNA','rRNA','tRNA']
		for genetype in RM_Genetype:
			out_dir = "modified_bed/"
			out_file = genetype +"RM.bed"
			something = everything[everything.iloc[:,3].str.contains(genetype)]
			something.to_csv(os.path.join(out_dir, out_file),sep='\t',index=False)
	
rule pandas_search_from_modified_bed_part2:
	input:
		bedfile = 'modified_bed/RNARM.bed'
#part 2 remove tRNA and rRNA
	output:'modified_bed/RNAonlyRM.bed'
	run:
		import pandas as pd
		import os

		everything = pd.read_csv(input.bedfile ,sep="\t+",header=None)
		out_dir = "modified_bed/"
		out_file = "RNAonlyRM.bed"
		rRNA = everything[everything.iloc[:,3].str.contains('rRNA')]
		tRNA = everything[everything.iloc[:,3].str.contains('tRNA')]
		something = pd.concat([everything,rRNA,tRNA]).drop_duplicates(keep=False)
		something.to_csv(os.path.join(out_dir, out_file),sep='\t',index=False, header=False)
	


