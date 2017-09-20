
FT_list= [line.rstrip('\n') for line in open('feature_type_list.txt')]
PT_list= [line.rstrip('\n') for line in open('parent_type_list.txt')]
Exon_PT_list= [line.rstrip('\n') for line in open('exon_parent_type_list.txt')]


rule all:
    input:
        'dm6.17_flybase.gff',
        'dm6.17_flybase_FTcounts.txt',
        expand('dm6.17_flybase_FT_{FT}.gff', FT=FT_list),
        expand('transcript_PT_identifiers/dm6.17_flybase_FT_{FT}.gff.PTID', FT=PT_list),
        expand('transcript_PT_identifiers/dm6.17_flybase_FT_{FT}.gff.geneID', FT=Exon_PT_list),
        expand('transcript_PT_identifiers/dm6.17_flybase_FT_{PT}.gff.PTID_FT_exon', PT=Exon_PT_list),
        'dm6.17_flybase_FT_exon_PTID.gff',
        expand('transcript_PT_identifiers/dm6.17_flybase_FT_gene_GeneType_{GT}.gff', GT=Exon_PT_list)

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


rule label_exon_files_with_PT_and_convert_to_bed:
    input:
        in_gff= expand('transcript_PT_identifiers/dm6.17_flybase_FT_{PT}.gff.PTID_FT_exon', PT=Exon_PT_list)
    output:
        out_gff= 'dm6.17_flybase_FT_exon_PTID.gff',
        #out_bed= 'dm6.17_flybase_FT_exon_PTID.bed'
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
        df.to_csv('dm6.17_flybase_FT_exon_PTID.gff',sep='\t', header=None, index=False)

#Note, we should also simultanously convert to .bed here


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
                        entry = entry + ";GeneType=" + current_GeneType
                        out_file.write(entry + '\n')

#  sanity check at this point seems good

#rule collapse_Genetypes_convert_to_bed:





















#for each exon PT, eg, mRNA, go to the FT gff, split ID field to obtain parent_gene_ID,
#collapse duplicates, then at this point, we have for each possible exon_PT, the list
#of gene IDS for each parenttype -> k done.

#next, input = dm6.17_flybase_FT_gene.gff. For each gene ID list, we iterate through
#input, for example, if geneID in {PT} ==  input.split(';')[0], add "Genetype={PT}".




#Note, dm6.17_flybase_FT_gene.gff.PTID, contains all the gene IDS? Yup i think so







#FBgn0259817
#FBti0059812
#db = gffutils.FeatureDB('dm6flybase.db', keep_order=True)
