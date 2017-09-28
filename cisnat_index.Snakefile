
FT_list= [line.rstrip('\n') for line in open('feature_type_list.txt')]
PT_list= [line.rstrip('\n') for line in open('parent_type_list.txt')]
Exon_PT_list= [line.rstrip('\n') for line in open('exon_parent_type_list.txt')]
Strand =['plus', 'minus']
Intersect_type=['exon','gene']
Intersect_samples =['dm6.17_flybase_FT_exon_PTID.merge', 'dm6.17_flybase_FT_gene_GeneID' ]

rule all:
    input:
        'dm6.17_flybase.gff','dm6.17_flybase_FTcounts.txt',
        expand('dm6.17_flybase_FT_{FT}.gff', FT=FT_list),
        expand('transcript_PT_identifiers/dm6.17_flybase_FT_{FT}.gff.PTID', FT=PT_list),
        expand('transcript_PT_identifiers/dm6.17_flybase_FT_{FT}.gff.geneID', FT=Exon_PT_list),
        expand('transcript_PT_identifiers/dm6.17_flybase_FT_{PT}.gff.PTID_FT_exon', PT=Exon_PT_list),
        'dm6.17_flybase_FT_exon_PTID.gff','dm6.17_flybase_FT_exon_PTID.bed',
        expand('transcript_PT_identifiers/dm6.17_flybase_FT_gene_GeneType_{GT}.gff', GT=Exon_PT_list),
        'dm6.17_flybase_FT_gene_GeneID.gff','dm6.17_flybase_FT_gene_GeneID.bed',
        'dm6.17_flybase_FT_exon_PTID.merge.bed',
        expand('dm6.17_flybase_FT_exon_PTID.merge.{strand}.bed', strand=Strand),
        expand('dm6.17_flybase_FT_gene_GeneID.{strand}.bed', strand=Strand),
        expand('{sample}.intersect.bed' , sample=Intersect_samples),
        expand('{sample}.intersect.ID.bed',sample=Intersect_samples),

rule extract_flybase:
    input:
        in_gff = 'dmel-all-r6.17.gff'
    output:
        out_gff = 'dm6.17_flybase.gff', out_FTcounts = 'dm6.17_flybase_FTcounts.txt'
    run:
        import pandas as pd

        #Note, end limit here has to be hardcoded? # line 1874 for DM6.17
        #Start: End of fasta, 23278924:25080161
        #First fasta is >211000022280479
        #last fasta is >211000022280192
        shell("sed '1,1874d' {input} | sed '23278924,25080161d' > dm6.17_temp.gff")
        gff = pd.read_csv('dm6.17_temp.gff', sep='\t', header=None, index_col=None)
        flybase = gff[gff[1].str.contains("^FlyBase", na=False)]
        #flybase[0] = flybase[0].str.replace('dmel_mitochondrion_genome', 'M')
        flybase[0] = 'chr' + flybase[0].astype(str)
        FTcounts = flybase.groupby([2]).size()
        #Note, fix this line such that it also gives names

        flybase.to_csv(output.out_gff, sep='\t', header=None, index=False)
        FTcounts.to_csv(output.out_FTcounts, sep='\t', header=None, index=False)

rule split_by_FT:
    input:
        in_gff = 'dm6.17_flybase.gff'
    output:
        out_gff = expand('dm6.17_flybase_FT_{FT}.gff', FT=FT_list)
    run:
        import pandas as pd

        gff = pd.read_csv(input.in_gff, sep='\t', header=None, index_col=None)

        for feature in FT_list:
            FT_gff = gff.loc[gff[2] == feature]
            FT_gff.to_csv("dm6.17_flybase_FT_" + str(feature) + ".gff", sep='\t',header=None, index=False)

rule extract_transcript_gene_ID:
    input:
        in_gff = expand('dm6.17_flybase_FT_{FT}.gff', FT=PT_list),
        in_gene_gff = expand('dm6.17_flybase_FT_{FT}.gff', FT=Exon_PT_list)
    output:
        out_ID = expand('transcript_PT_identifiers/dm6.17_flybase_FT_{FT}.gff.PTID', FT=PT_list),
        out_gene_ID = expand('transcript_PT_identifiers/dm6.17_flybase_FT_{FT}.gff.geneID', FT=Exon_PT_list)
    run:
        import pandas as pd
        shell("mkdir -p 'transcript_PT_identifiers'")

        for gff in input.in_gff:
            gffdata = pd.read_csv(gff,sep='\t', header=None, index_col=None)
            split_PTID_field = (gffdata[8].str.split(';',expand=True))[0]
            split_geneID_field = (gffdata[8].str.split(';',expand=True))[2]
            PTID = (split_PTID_field.str.split('=',expand=True))[1]
            geneID = (split_geneID_field.str.split('=',expand=True))[1]
            geneID = geneID.drop_duplicates()
            geneID = 'ID=' + geneID
            PTID.to_csv(os.path.join('transcript_PT_identifiers/', gff + '.PTID'),sep='\t', header=None, index=False)
            geneID.to_csv(os.path.join('transcript_PT_identifiers/', gff + '.geneID'),sep='\t', header=None, index=False)

#Note, for .geneID, transposableelement.gff and gene.gff are junk and therefore not in Exon_PT_list
#Note, problem for snoRNA.geneID? We are short of 1? Ok, nvm i think we ignore this problem first,
#fix it later.
rule define_exon_PT:
    input:
        in_ID = expand('transcript_PT_identifiers/dm6.17_flybase_FT_{PT}.gff.PTID', PT=Exon_PT_list),
        in_db = 'dm6flybase.db'
    output:
        out_gff = expand('transcript_PT_identifiers/dm6.17_flybase_FT_{PT}.gff.PTID_FT_exon', PT=Exon_PT_list)
    run:
        import pandas as pd
        import gffutils
        db = gffutils.FeatureDB(input.in_db, keep_order=True)

        for PTID in input.in_ID:
            IDs = [line.rstrip('\n') for line in open(PTID)]
            sys.stdout = open( PTID + "_FT_exon", "w")
            for ID in IDs:
                children = db.children(ID, featuretype='exon')
                for feature in children:
                    print (feature)


rule combine_exon_files_with_PTID_and_convert_to_bed:
    input:
        in_gff= expand('transcript_PT_identifiers/dm6.17_flybase_FT_{PT}.gff.PTID_FT_exon', PT=Exon_PT_list)
    output:
        out_gff= 'dm6.17_flybase_FT_exon_PTID.gff',
        out_bed= 'dm6.17_flybase_FT_exon_PTID.bed'
    run:
        import pandas as pd
        counter = 0
        df = pd.DataFrame()
        for gff in input.in_gff:
            currentPT = Exon_PT_list[counter]
            gffdata = pd.read_csv(gff,sep='\s+',header=None, index_col=None)
            dropdata = gffdata.drop_duplicates()
            dropdata[9] = "Parenttype=" + currentPT
            df = df.append(dropdata)
            counter += 1
        tempdf = df
        tempdf[10] = tempdf[8].astype(str) + ";" + tempdf[9].astype(str)
        bed = tempdf.loc[:, [0,3,4,10,5,6]]
        bed[7] = bed[bed.columns[2]].astype(int) - bed[bed.columns[1]].astype(int)
        bed = bed.drop(bed[bed[7] <= 0].index)
        bed[10] = bed[10].astype(str) + ";Exonlength=" +  bed[7].astype(str)
        finalbed = bed[bed.columns[:-1]]
        df.to_csv('dm6.17_flybase_FT_exon_PTID.gff',sep='\t', header=None, index=False)
        finalbed.to_csv('dm6.17_flybase_FT_exon_PTID.bed',sep='\t', header=None, index=False)

#bedformat = chr start stop ID . Strand
#Note, removing exon whose exonlength==0
#oh, my, goodness, the .loc shit.
#Error on line 8546,30528 in dm6.17_flybase_FT_exon_PTID.bed.
#2L      FlyBase exon 11142647  11142647 . - .  Parent=FBtr0343625,FBtr0346610
#Note, this is present in the original gff file.
#Note, it is at this point that i cheat, i am so fucking done.



rule label_gene_with_Genetype:
    input:
        in_ID = expand('transcript_PT_identifiers/dm6.17_flybase_FT_{FT}.gff.geneID', FT=Exon_PT_list),
        in_gff = 'dm6.17_flybase_FT_gene.gff'
    output:
        out_gff = expand('transcript_PT_identifiers/dm6.17_flybase_FT_gene_GeneType_{GT}.gff', GT=Exon_PT_list)
    run:
        import pandas as pd
        counter = 0
        gff = [line.rstrip('\n') for line in open(input.in_gff)]

        for in_file in input.in_ID:
            current_GeneType = Exon_PT_list[counter]
            ID_list = [line.rstrip('\n') for line in open(in_file)]
            out_file = open(os.path.join('transcript_PT_identifiers/','dm6.17_flybase_FT_gene_GeneType_' + current_GeneType + '.gff'),"w")
            counter += 1
            for ID in ID_list:
                for entry in gff:
                    if ID in entry:
                        entry = entry + "\t" + "GeneType=" + current_GeneType
                        out_file.write(entry + '\n')

#sanity check at this point seems good
#Note my naming conventions are horrible

rule combine_Genetypes_convert_to_bed:
    input:
        in_gff = expand('transcript_PT_identifiers/dm6.17_flybase_FT_gene_GeneType_{GT}.gff', GT=Exon_PT_list)
    output:
        out_gff = 'dm6.17_flybase_FT_gene_GeneID.gff', out_bed = 'dm6.17_flybase_FT_gene_GeneID.bed'
    run:
        import pandas as pd
        c_gff = pd.DataFrame()

        for gff in input.in_gff:
            gffdata = pd.read_csv(gff,sep='\t',header=None, index_col=None)
            c_gff = c_gff.append(gffdata)

        field = (c_gff[8].str.split(';',expand=True))
        ID = field[field.columns[0]]
        Strand= ";GnStrand=" + c_gff[c_gff.columns[6]].astype(str)
        Start = ";GnStart=" + c_gff[c_gff.columns[3]].astype(str)
        Length = ";Genelength=" + (c_gff[c_gff.columns[4]].astype(int) - c_gff[c_gff.columns[3]].astype(int)).astype(str)
        GeneType = ";" + c_gff[c_gff.columns[-1]].astype(str)
        ID_field = ID.astype(str) + Strand.astype(str) + Start.astype(str) + Length.astype(str) + GeneType.astype(str)
        bed = c_gff[[0,3,4,5,6]]
        bed.insert(3, column = 'ID', value = ID_field)

        c_gff.to_csv(output.out_gff,sep='\t', header=None, index=False)
        bed.to_csv(output.out_bed,sep='\t', header=None, index=False)

#As expected im 1 gene short in the snoRNA

rule merge_exon_convert_bed:
    input:
        in_exon = 'dm6.17_flybase_FT_exon_PTID.bed'
    output:
        out_merge = 'dm6.17_flybase_FT_exon_PTID.merge.bed',
    run:
        import pandas as pd
        shell("sort-bed {input} | bedtools merge -s -c 4 -o collapse -i - | \
        sort-bed - > temp.merge")

        merge_exon = pd.read_csv('temp.merge',sep='\t',header=None, index_col=None)
        merge_exon_bed = merge_exon[[0,1,2,4,3]]
        merge_exon_bed.insert(4, column = 'temp', value = '.')
        merge_exon_bed.to_csv(output.out_merge,sep='\t', header=None, index=False)

        shell("rm temp.merge")

rule extract_strand_for_exon_gene:
    input:
        in_exon = 'dm6.17_flybase_FT_exon_PTID.merge.bed',
        in_gene = 'dm6.17_flybase_FT_gene_GeneID.bed'
    output:
        out_exon = expand('dm6.17_flybase_FT_exon_PTID.merge.{strand}.bed', strand=Strand),
        out_gene = expand('dm6.17_flybase_FT_gene_GeneID.{strand}.bed', strand=Strand)
    run:
        import pandas as pd

        exon_bed = pd.read_csv(input.in_exon,sep='\t',header=None, index_col=None)
        gene_bed = pd.read_csv(input.in_gene,sep='\t',header=None, index_col=None)
        counter = 0
        for bed in (exon_bed, gene_bed):
            plus_features, minus_features = bed.loc[bed[5] == '+'], bed.loc[bed[5] == '-']
            counter += 1
            print(counter)
            if counter == 1 :
                plus_features.to_csv('dm6.17_flybase_FT_exon_PTID.merge.plus.bed',sep='\t', header=None, index=False)
                minus_features.to_csv('dm6.17_flybase_FT_exon_PTID.merge.minus.bed',sep='\t', header=None, index=False)
            else:
                plus_features.to_csv('dm6.17_flybase_FT_gene_GeneID.plus.bed',sep='\t', header=None, index=False)
                minus_features.to_csv('dm6.17_flybase_FT_gene_GeneID.minus.bed',sep='\t', header=None, index=False)

rule intersect_within_gene_exon_remove_gene_duplicates:
    input:
        plus = "{sample}.plus.bed",
        minus = "{sample}.minus.bed"

    output:
        '{sample}.intersect.bed'
    run:
        import pandas as pd
        import os.path

        shell("bedtools intersect -S -a {input.plus} -b {input.minus} -bed | \
        sort-bed - > {output}")
        #Note, here is how you coordinate your input files, and you dont need wildcards.
        #however, you do need to expand a list of your output files in rule all though.
        gene_bed = 'dm6.17_flybase_FT_gene_GeneID.intersect.bed'
        if os.path.isfile(gene_bed):
            intersect_data = pd.read_csv(gene_bed,sep='\t',header=None, index_col=None)
            intersect_data = intersect_data.drop_duplicates()
            intersect_data.to_csv(gene_bed,sep='\t', header=None, index=False)
        else:
            pass

#nte, the resultant gene.intersect file has duplicate features why?
#But the resultant exon file does not.
#Oh my god, editing one file in the early chain -> changes everything dowstream
#Note .columns is for single columns, while [[]] is for multiple columnds

rule map_to_original_bed_for_intersect_ID:
    input:
        intersect = "{sample}.intersect.bed",
        bed = "{sample}.bed"
    output:
        '{sample}.intersect.ID.bed'

    shell:"""
        bedtools map -c 4 -o collapse -a {input.intersect} -b {input.bed} | \
        sort-bed - > {output}
        """

#Initial to do list for tmr:
#First thing, reformat the original exon file to include additional details.
#Next, re-rerun the exon part, to see if everything changes after that part.
#Then address the unsorted geneID issue.
#Sort order was unspecified, and file dm6.17_flybase_FT_gene_GeneID.bed is not sorted lexicographically.












#Note how bedtools merge will always mess up the bed format.


#wait a second, how come when i merge the gene.bed, there are merges?
#that means, there are actually overlapping genes?
#rule
















#for each exon PT, eg, mRNA, go to the FT gff, split ID field to obtain parent_gene_ID,
#collapse duplicates, then at this point, we have for each possible exon_PT, the list
#of gene IDS for each parenttype -> k done.

#next, input = dm6.17_flybase_FT_gene.gff. For each gene ID list, we iterate through
#input, for example, if geneID in {PT} ==  input.split(';')[0], add "Genetype={PT}".




#Note, dm6.17_flybase_FT_gene.gff.PTID, contains all the gene IDS? Yup i think so
