
import os
gff_file_names = ["../Databases/Oryza_sativa.IRGSP-1.0.44.chr.gff3","../Databases/MSU7.gff3"]

for gff_file_name in gff_file_names:
    f = open(gff_file_name,"r")
    outfolder_location = "../Databases/cached_gff_chunks/"
    if not os.path.exists(outfolder_location):
        os.makedirs(outfolder_location)

    #1	RAP2018-11-26	CDS	4357	4455	.	+	0	ID=CDS:Os01t0100100-01;Parent=transcript:Os01t0100100-01;protein_id=Os01t0100100-01
    #1	RAP2018-11-26	exon	5457	5560	.	+	.	Parent=transcript:Os01t0100100-01;Name=Os01t0100100-01-E4;constitutive=1;ensembl_end_phase=2;ensembl_phase=0;exon_id=Os01t0100100-01-E4;rank=4
    #1	RAP2018-11-26	CDS	5457	5560	.	+	0	ID=CDS:Os01t0100100-01;Parent=transcript:Os01t0100100-01;protein_id=Os01t0100100-01
    #1	RAP2018-11-26	exon	7136	7944	.	+	.	Parent=transcript:Os01t0100100-01;Name=Os01t0100100-01-E5;constitutive=1;ensembl_end_phase=1;ensembl_phase=2;exon_id=Os01t0100100-01-E5;rank=5
    #1	RAP2018-11-26	CDS	7136	7944	.	+	1	ID=CDS:Os01t0100100-01;Parent=transcript:Os01t0100100-01;protein_id=Os01t0100100-01
    #Chr1	MSU_osa1r7	CDS	4357	4455	.	+	.	ID=LOC_Os01g01010.1:cds_2;Parent=LOC_Os01g01010.1
    #Chr1	MSU_osa1r7	CDS	5457	5560	.	+	.	ID=LOC_Os01g01010.1:cds_3;Parent=LOC_Os01g01010.1
    #Chr1	MSU_osa1r7	CDS	7136	7944	.	+	.	ID=LOC_Os01g01010.1:cds_4;Parent=LOC_Os01g01010.1
    #Chr1	MSU_osa1r7	CDS	8028	8150	.	+	.	ID=LOC_Os01g01010.1:cds_5;Parent=LOC_Os01g01010.1
    #Chr1	MSU_osa1r7	CDS	8232	8320	.	+	.	ID=LOC_Os01g01010.1:cds_6;Parent=LOC_Os01g01010.1

    counter = 0
    prev_transcript = "NULL"
    f_out = open("temp_file","w")#will never write here

    for line in f:
        line = line[:-1]
        cells = line.split("\t")

        if len(cells) > 5:
            feature_type = cells[2]
            if feature_type == "CDS":

                transcript_id = cells[8].replace("CDS:","").replace("ID=","").replace(":",";").split(";")[0] #ID=CDS:Os01t0100100-01; ID=LOC_Os01g01010.1:
                if transcript_id != prev_transcript:    #new locus to process
                    f_out.close()
                    f_out = open(outfolder_location+transcript_id+".gff","w")
                    prev_transcript = transcript_id
                f_out.write(line+"\n")
                counter += 1

        #for testing
        #if counter > 100:
        #    exit()
    f_out.close()
    f.close()





