#Snakefile for building Indexes
#raw fasta in up in 1 directory

#---------demo-----------
#rule quantify_genes:
#    input:
#        genome = 'genome.fa',
#        r1 = 'fastq/{sample}.R1.fastq.gz',
#        r2 = 'fastq/{sample}.R2.fastq.gz'
#    output:
#        '{sample}.txt'
#    shell:
#        'echo {input.genome} {input.r1} {input.r2} > {output}'
#-----------demo------------


#remove unmapped scaffolds from genome
rule extra_chr_from_genome:
	input:
		genome = "../FastaLoqs/dmel-all-chromosome-r6.16.fasta",
		#list of chr to include
		#cat Chr_list.txt | xargs > Chrlist.txt
		Chr_list = "Chr_list.txt"
	output:
		"../FastaLoqs/2L.temp"
	shell:
		"""
		samtools faidx {input.genome}
#		cat {input.Chr_list} | while read line
#		do
#		"samtools faidx {input.genome} $line > {output}.$line.temp"
#		done
		samtools faidx {input.genome} 2L > {output}
		"""

#combine and remove temp to make index0
rule make_index0:
	input:
		"../FastaLoqs/{sample}.temp"
	output:
		"../FastaLoqs/Index0.fa"
		#List = "output.done"
	shell:
		"""
		cat {wildcards.sample} > {output}
		rm {wildcards.sample}
		"""

