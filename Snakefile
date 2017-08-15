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
rule remove_U_from_genome:
	input:
		genome = "../FastaLoqs/dmel-all-chromosome-r6.16.fasta",
		#list of chr to include
		#cat Chr_list.txt | xargs > Chrlist.txt
		#Chr_list = "Chrlist.txt"
	output:
		"../FastaLoqs/Index0.fa"
		#List = "output.done"
	shell:
		"""
		samtools faidx {input.genome}
		
		samtools faidx {input.genome} \
2L
 2R
 3L
 3R
 4
 X
 Y
 rDNA
 2Cen_mapped_Scaffold_10_D1684
 2Cen_mapped_Scaffold_43_D1668
 2R2_mapped_Scaffold_56_D1828
 3Cen_mapped_Scaffold_1_D1896_D1895
 3Cen_mapped_Scaffold_27_D1777
 3Cen_mapped_Scaffold_31_D1643_D1653_D1791
 3Cen_mapped_Scaffold_36_D1605
 3Cen_mapped_Scaffold_41_D1641
 3Cen_mapped_Scaffold_50_D1686
 X3X4_mapped_Scaffold_14_D1732
 X3X4_mapped_Scaffold_6_D1712
 XY_mapped_Scaffold_42_D1648
 XY_mapped_Scaffold_7_D1574
 Y_mapped_Scaffold_12_D1771
 Y_mapped_Scaffold_15_D1727
 Y_mapped_Scaffold_18_D1698
 Y_mapped_Scaffold_20_D1762_D1719
 Y_mapped_Scaffold_21_D1683_D1693
 Y_mapped_Scaffold_23_D1638
 Y_mapped_Scaffold_26_D1717
 Y_mapped_Scaffold_30_D1720
 Y_mapped_Scaffold_34_D1584
 Y_mapped_Scaffold_53_D1765
 Y_mapped_Scaffold_5_D1748_D1610
 Y_mapped_Scaffold_9_D1573 > {output}
		"""


