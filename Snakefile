#helppage
#http://slowkow.com/notes/snakemake-tutorial
#add comments

#snakemake --dag | dot -Tsvg > dag.svg
#LOGS=Logs_out
#snakemake -s snakefile_part_2 --jobname '$FOLDER/s.{rulename}.{jobid}' --stats $LOGS/snakemake.txt >& $LOGS/snakemake.log

import os
 
CHR_LIST= [line.rstrip('\n') for line in open('Chr_list.txt')]
RM_GENETYPES=['LTR','LINE','DNA','Satellite','Low_complexity','RC','Simple_repeat','Other','Unknown','ARTEFACT','RNA','rRNA','tRNA']
Flybase_GENETYPES2=['rRNA', 'tRNA', 'snRNA', 'snoRNA']
#requires installation
RepeatMasker_Dir='/usr/local/RepeatMasker/RepeatMasker'
Modified_bed_files=[]

fasta_file_path = []
	
#tab or space matters
rule all:
    input:
        'Index0.fa',
        'RM_output',
        'modified_bed/RM_output.bed',
        expand('modified_bed/{RM_Genetype}RM.bed', RM_Genetype=RM_GENETYPES),
        'modified_bed/RNAonlyRM.bed',
        'modified_bed/pri_miRNA_mirbase21_dm6.bed',
        'modified_bed/fasta/test.txt',
        'modified_bed/fasta/test2.txt',
        'modified_bed/fasta/Index0.fa',
        'modified_bed/fasta/pri_miRNA.fa',
        'modified_bed/fasta/flybasehpRNA.fa',
        'modified_bed/fasta/bt_indexes/test.txt'


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
		out_dir ='modified_bed',
		bedfiles='modified_bed/RM_output.bed'
	shell: """
	mkdir -p {output.out_dir}
	awk -f ./RMout2bed.awk {input} > {output.bedfiles}
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
			something = everything[everything.iloc[:,3].str.contains(genetype)]
			something.to_csv(os.path.join(out_dir, out_file),sep='\t',index=False)
		
#might cause a bug

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

#rule copy_bed_files from jinwee:
#	input:
#		source='cis_nat_snakefile_package/merged_bedfiles'
#	output:
#		outfile='modified_bed'
#	shell:"""
#	cp {input.source}/*merge.bed {output.outfile}
#	"""

rule merge_bed_to_fasta:
#only for GFF from jinwee
	input:
		genome = 'Index0.fa'

	output:
		test_text='modified_bed/fasta/test2.txt'
	run:
		import os
		from pathlib import Path
		Merge_bed=[]
		out_text_file=open(output.test_text,"w")
		out_dir=Path('modified_bed/fasta')
		if out_dir.is_dir() != True:
			os.mkdir('modified_bed/fasta')
		for file in os.listdir(os.path.join(os.getcwd(),"modified_bed")):
			if file.endswith("merge.bed"):
				Merge_bed.append(file)
		for i in range(len(Merge_bed)):
			bed_file=os.path.join('modified_bed/',Merge_bed[i])
			out_file=os.path.join('modified_bed/fasta',Merge_bed[i][18:-10]+'.fa')
			shell('''awk -f ./RemoveChr.awk {bed_file} > {bed_file}.nochr 
				bedtools getfasta -name -fi {input.genome} -bed {bed_file}.nochr -fo {out_file}''')
			out_text_file.write(Merge_bed[i]+'\n')
		out_text_file.close()
	
rule modified_bed_to_fasta:
#only for RM
	input:
		genome = 'Index0.fa'

	output:
		test_text='modified_bed/fasta/test.txt'
	run:
		import os
		from pathlib import Path
		out_text_file=open(output.test_text,"w")
		out_dir=Path('modified_bed/fasta')
		if out_dir.is_dir() != True:
			os.mkdir('modified_bed/fasta')
		for file in os.listdir(os.path.join(os.getcwd(),"modified_bed")):
			if file.endswith("RM.bed"):
				Modified_bed_files.append(file)
		for i in range(len(Modified_bed_files)):
			bed_file=os.path.join('modified_bed/',Modified_bed_files[i])
			out_file=os.path.join('modified_bed/fasta',Modified_bed_files[i][:-6]+'.fa')
			shell('bedtools getfasta -name -fi {input.genome} -bed {bed_file} -fo {out_file}')
			out_text_file.write(Modified_bed_files[i]+'\n')
		out_text_file.close()
	
	
#do we want the individual fasta entry to be named?
#make modified_bed
rule remove_files:
	output: 'modified_bed/RNAoRM.bed'
	shell:'''
	rm modified_bed/RNARM.bed
	mv modified_bed/RNAonlyRM.bed {output}
	rm modified_bed/tRNARM.bed
	'''

rule move_other_fastas:
	input:
		genome = 'Index0.fa'
	output:
		'modified_bed/fasta/Index0.fa'
	shell:'''
	cp {input.genome} {output}
	'''

rule hairpin_fa:
	input:
		genome = 'Index0.fa',
		main='cis_nat_snakefile_package/dm6.17_flybase.gff'
	output:
		bed_file_temp='modified_bed/flybasehpRNA.bed.temp',
		bed_file_temp_2='modified_bed/flybasehpRNA.bed.2.temp',
		bed_file='modified_bed/flybasehpRNA.bed',
		fasta_file='modified_bed/fasta/flybasehpRNA.fa'
	params:
		grep='Name=hpRNA'

	shell:'''
	grep {params.grep} {input.main} > {output.bed_file_temp}
	awk -f ./gff_to_bed_hpRNA.awk {output.bed_file_temp} > {output.bed_file_temp_2}
	awk -f ./RemoveChr.awk {output.bed_file_temp_2} > {output.bed_file}
	bedtools getfasta -name -fi {input.genome} -bed {output.bed_file} -fo {output.fasta_file}
	'''

rule make_bt_indexes:
	input:
		fasta_file=expand(fasta_file_path)
	output:
		test_text='modified_bed/fasta/bt_indexes/test.txt'
	run:
		import os
		out_text_file=open(output.test_text,"w")
		fasta_folder='modified_bed/fasta'
		fasta_file_path=[]
		for file in os.listdir(os.path.join(os.getcwd(),fasta_folder)):
			if file.endswith(".fa"):
				fasta_file_path.append(os.path.join(fasta_folder,file))
		rename_dict={'Index0.fa':'Index0','rRNA.fa':'Index1','exon_PT_rRNA.fa':'Index2','RNAonly.fa':'Index4','exon_PT_tRNA.fa':'Index5','exon_PT_snRNA.fa':'Index6','exon_PT_snoRNA.fa':'Index7','pri_miRNA.fa':'Index8','flybasehpRNA.fa':'Index9','LTR.fa':'Index10','LINE.fa':'Index11','DNA.fa':'Index12','Satellite.fa':'Index13','Low_complexity.fa':'Index14','RC.fa':'Index15','Simple_repeat.fa':'Index16','Other.fa':'Index17','Unknown.fa':'Index18','transposable_element.fa':'Index19','ARTEFACT.fa':'Index20'}
		for file in fasta_file_path:
			if file[len(fasta_folder)+1:] not in rename_dict:
				continue
			else:
				prefix=os.path.join('modified_bed/fasta/bt_indexes',rename_dict[file[(len(fasta_folder)+1):]])
				shell('bowtie-build {file} {prefix}')
				out_text_file.write(file+'\n')
		out_text_file.close()