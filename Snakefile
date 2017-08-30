#helppage
#http://slowkow.com/notes/snakemake-tutorial

#Index3.fa obtained from NCBI M21017.1
import os
 
CHR_LIST= [line.rstrip('\n') for line in open('Chr_list.txt')]
RM_GENETYPES=['LTR','LINE','DNA','Satellite','Low_complexity','RC','Simple_repeat','Other','Unknown','ARTEFACT','RNA','rRNA','tRNA']
Flybase_GENETYPES2=['rRNA', 'tRNA', 'snRNA', 'snoRNA']
#requires installation
RepeatMasker_Dir='/usr/local/RepeatMasker/RepeatMasker'
Modified_bed_files=[]

for file in os.listdir(os.path.join(os.getcwd(),"modified_bed")):
	if file.endswith("RM.bed"):
		Modified_bed_files.append(file[:-6])

import pandas as pd

rule all:
    input:
        'Index0.fa',
        'RM_output',
        'modified_bed/RM_output.bed',
        expand('modified_bed/{RM_Genetype}RM.bed', RM_Genetype=RM_GENETYPES),
        'modified_bed/RNAonlyRM.bed',
        expand('modified_bed/fasta/{genetype}.fa', genetype=Modified_bed_files),
        'modified_bed/fasta/bowtieIndexes'
        
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

##liftover does not work
#rule dme3gff3_to_bed:
#	input:
#		dme3gff3 = 'Liftover_files/dme.gff3' 
#	output:
#		dme3bed	='Liftover_files/mirbasedme3.bed'
#	shell: """
#	awk -f ./Gff3toBed.awk {input.dme3gff3} > {output.dme3bed}
#	"""

#install UCSC liftover, get chainfile from flybase, unzip with gzip -d	

#rule liftOver_mirbase21_dm3gff_todm6:
#	input:
#		dme3bed	='Liftover_files/mirbasedme3.bed',
#		chainfile = 'Liftover_files/dm3ToDm6.over.chain'
#	output:
#		newFile = 'Liftover_files/mirbase21_dm6.bed',
#		unMapped ='Liftover_files/liftover_unmapped',
#		
#	shell:"""
#	liftOver {input.dme3bed} {input.chainfile} {output.newFile} {output.unMapped}
#	"""

rule modified_bed_to_fasta:
#only for RM
	input:
		genome = 'Index0.fa'

	output:
		'modified_bed/fasta/{genetype}.fa'

	shell:"""
	mkdir -p 'modified_bed/fasta'
	bedtools getfasta -fi {input.genome} -bed 'modified_bed/{wildcards.genetype}RM.bed' -fo {output}
	"""
##rule to rename Indexes into numbers

#Index_num to change
#rule fasta_to_bowtieindex:
#	input:
#		fasta_folder = 'modified_bed/fasta'
#	output:
#		expand("modified_bed/fasta/Indexes{Index_num}.{index}.bt2", index=range(1,5), Index_num=Modified_bed_files),
#       expand("modified_bed/fasta/Indexes{Index_num}.rev.{index}.bt2", index=range(1,3),Index_num=Modified_bed_files)
#
#	run:
#		import os
#		fasta_file = []
#		for file in os.listdir(os.path.join(os.getcwd(),input.fasta_folder)):
#			if file.endswith(".fa"):
#				fasta_file.append(file)
#		for i in fasta_file:
#			bowtie-build os.path.join(os.getcwd(),input.fasta_folder,i) i[:-3]
	
	
