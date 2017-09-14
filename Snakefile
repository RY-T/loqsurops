#helppage
#http://slowkow.com/notes/snakemake-tutorial
#add comments

#snakemake --dag | dot -Tsvg > dag.svg
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

fasta_folder='modified_bed/fasta'
fasta_file_path = []

for file in os.listdir(os.path.join(os.getcwd(),fasta_folder)):
	if file.endswith(".fa"):
		fasta_file_path.append(os.path.join(fasta_folder,file))	
#tab or space matters
rule all:
    input:
        'Index0.fa',
        'RM_output',
        'modified_bed/RM_output.bed',
        expand('modified_bed/{RM_Genetype}RM.bed', RM_Genetype=RM_GENETYPES),
        'modified_bed/RNAonlyRM.bed',
        'modified_bed/pri_miRNA_mirbase21_dm6.bed',
        expand('modified_bed/fasta/{genetype}.fa', genetype=Modified_bed_files),
        'modified_bed/fasta/bt_indexes',
        'modified_bed/fasta/pri_miRNA.fa'

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
	mkdir -p {output} 
	{input.RM_Dir} -dir {output} -species drosophila {input.genome}
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
			something = everything[everything.iloc[:,3].str.contains(genetype+'/')]
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

#liftover does not work
#Gff3toBed also selects for pri-miRNA
rule liftover_mirbase21_to_dm6bed:
	input:
		dme3gff3 = 'Liftover_files/dme.gff3',
		dme3bed	='Liftover_files/mirbasedme3.bed',
		chainfile = 'Liftover_files/dm3ToDm6.over.chain',
		genome = 'Index0.fa' 

	output:
		dme3bed	='Liftover_files/mirbasedme3.bed',
		newFile = 'Liftover_files/mirbase21_dm6.bed',
		unMapped ='Liftover_files/liftover_unmapped',
		final_bed= 'modified_bed/pri_miRNA_mirbase21_dm6.bed',
		fasta_out='modified_bed/fasta/pri_miRNA.fa'
	shell: """
	awk -f ./Gff3toBed.awk {input.dme3gff3} > {output.dme3bed}
	liftOver {input.dme3bed} {input.chainfile} {output.newFile} {output.unMapped}
	awk -f ./RemoveChr.awk {output.newFile} > {output.final_bed}
	bedtools getfasta -name -fi {input.genome} -bed {output.final_bed} -fo {output.fasta_out}
	"""

#install UCSC liftover, get chainfile from flybase, unzip with gzip -d	

	
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
#do we want the individual fasta entry to be named?

##rule to rename Indexes into numbers

#Index_num to change

rule make_bt_indexes:
	input:
		fasta_file=expand(fasta_file_path)
	output:
		'modified_bed/fasta/bt_indexes'
	run:
		import os
		if os.path.isfile('modified_bed/fasta/bt_indexes') != True:
			os.mkdir('modified_bed/fasta/bt_indexes')
		fasta_file=fasta_file_path
		Fasta_folder=fasta_folder
		for file in fasta_file:
			prefix=os.path.join('modified_bed/fasta/bt_indexes',file[(len(Fasta_folder)+1):-3])
			shell('bowtie-build {file} {prefix}')