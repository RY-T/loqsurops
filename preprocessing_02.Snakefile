import pandas as pd
#Env list: loqsproject  panda_env  snake_tutorial

run_list = pd.read_csv('temp_runlist.txt',sep='\t', header=None,index_col=None)
lib_list = run_list[0]
adapter = 'TGGAATTCTCGGGTGCCAAGG'
qualscore = 'Q33'
min_length = 18
max_length = 30

rule all:
    input:
        expand('02preprocessouput/{LIB}_final_col_filter.fasta', LIB=lib_list)

rule clip_convert_to_fasta_and_identifiers:
    input:
        '{LIB}.fastq'
    output:
        '{LIB}clip.fasta'
    shell:"""
        fastx_clipper -a {adapter} -c -v -i {input} |
        fastq_to_fasta -v -r | sed "s/^>/>{wildcards.LIB}_/" > {output}
        """

rule collapsing_repeats:
    input:
        '{LIB}clip.fasta'
    output:
        '{LIB}clipcolid_unlabeled.fasta'
    shell:"""
        fastx_collapser -v -i {input} -o {output}
        rm {input}
        """
#NOTE, here usage of wildcards is key

rule labelling_converting_tab:
    input:
        '{LIB}clipcolid_unlabeled.fasta'
    output:
        '{LIB}.tab'
    shell:"""
        fasta_formatter -t -i {input} |
        sed "s/-/_/g ; s/_/_{wildcards.LIB}_/" > {output}
        """

rule filter_length:
    input:
        intab = '{LIB}.tab'
    output:
        outtab = '{LIB}len.tab'
    run:
        import pandas as pd

        data = pd.read_csv(input.intab, sep='\t', names = ["Identifier", "Read"])
        data['length'] = data['Read'].str.len()
        filtered_data = data.query('18 <= length <= 30')
        filtered_data.to_csv(output.outtab, sep='\t', header=None, index=False)

#Note:
#1) superimpt to define output and input as variable
#2) pandas only runs on disk, so we have to write all our stuff to a file
#3) subtle difference in index / index col between to_csv and read_csv

rule convert_to_fasta:
    input:
        intab = '{LIB}len.tab'
    output:
        outfasta = '{LIB}len.fasta'
    run:
        import pandas as pd
        import re

        text_file = open(output.outfasta,"w")

        data = pd.read_csv(input.intab, sep='\t', names = ["Identifier", "Read", "length"])
        fasta_data = data[['Read', 'Identifier']]
        fasta_data['Read'] = '-' + fasta_data['Read'].astype(str)
        fasta_data['Identifier'] = '!' + fasta_data['Identifier'].astype(str)
        fasta_data['fasta'] = fasta_data['Identifier'] + fasta_data['Read']
        fasta = fasta_data['fasta']
        string = fasta.to_string(header=False)

        for line in string.split('\n'):
            fields = re.split('!|-',line)
            text_file.write((">{}\n{}\n".format(fields[1],fields[2])))

# oh my goodness

rule artifact_filter:
    input:
        '{LIB}len.fasta'
    output:
        '{LIB}_final_col_filter.fasta'
    shell:"""
        fastx_artifacts_filter -v -i {input} > {output}
        echo "Final read identiier count:"
        grep -c {wildcards.LIB} {output}
        """

rule moving_files:
    input:
        '{LIB}_final_col_filter.fasta'
    output:
        '02preprocessouput/{LIB}_final_col_filter.fasta'
    shell:"""
        mkdir -p '02preprocessouput'
        cp {input} {output}
        """
    















#fastx_artifacts_filter -v -i - -o {output}
