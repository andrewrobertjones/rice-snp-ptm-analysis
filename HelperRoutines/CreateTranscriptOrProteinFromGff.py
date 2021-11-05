#Simple script to extract coordinates for a gene model from GFF
#to make transcript sequence and protein sequence

from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord

PATH_TO_DATABASE_FOLDER = "../../RiceDatabases/"


### Argument 1: input_gff file
### 2: BioSeq object containing single chromosome
### 3: Locus ID to extract from GFF
### 4: Returns a list of exons, ready for BLAST search
def extract_list_of_transcript_seq(input_gff,bio_chr_seq,locus_id):
    transcript_seqs = []

    f= open(input_gff,"r")

    chr_seq = bio_chr_seq.seq
    start_to_end = {}

    strand = "+"
    for line in f:
        line = line[:-1]
        cells = line.split("\t")
        strand = cells[6]

        if cells[2].upper() == "CDS" and locus_id in line:
            start = int(cells[3])
            end = int(cells[4])
            #transcript_seq += str(chr_seq[start:end])
            start_to_end[start] = end

    #Need to make sure exons are sorted small to large before next step
    sorted_starts = sorted(start_to_end.keys())

    cds_counter = 1
    for start in sorted_starts:
        end = start_to_end[start]
        new_cds_record = SeqRecord.SeqRecord(chr_seq[start-1:end],id=locus_id + "_cds" + str(cds_counter))
        cds_counter += 1
        transcript_seqs.append(new_cds_record)

    #bio_transcript_seq = Seq.Seq(transcript_seq,generic_dna)
    #if strand == "-":   #reverse
    #    bio_transcript_seq = bio_transcript_seq.reverse_complement()

    return transcript_seqs

### Argument 1: input_gff file
### 2: Locus ID to extract from GFF
### Function just returns the start and end regions of the coding sequence
def get_start_end_positions_from_gff(input_gff,locus_id):
    f = open(input_gff, "r")

    coding_start_and_end = []

    strand = "+"
    coding_start = 1E50
    coding_end = -1
    for line in f:
        line = line[:-1]
        cells = line.split("\t")
        strand = cells[6]

        if cells[2].upper() == "CDS" and locus_id in line:
            start = int(cells[3])
            end = int(cells[4])
            # transcript_seq += str(chr_seq[start:end])
            if start < coding_start:
                coding_start = start
            if end > coding_end:
                coding_end = end

    coding_start_and_end.append(coding_start)
    coding_start_and_end.append(coding_end)
    return coding_start_and_end



### Argument 1: input_gff file
### 2: BioSeq object containing single chromosome
### 3: Locus ID to extract from GFF
def make_transcript_seq(input_gff,bio_chr_seq,locus_id):
    #print("Running method with no offset")
    transcript_seq = ""

    f= open(input_gff,"r")

    #print("opening:",input_gff)

    chr_seq = bio_chr_seq.seq
    start_to_end = {}

    strand = "+"
    for line in f:
        line = line[:-1]
        cells = line.split("\t")
        strand = cells[6]

        if cells[2].upper() == "CDS" and locus_id in line:
            start = int(cells[3])
            end = int(cells[4])
            #transcript_seq += str(chr_seq[start:end])
            start_to_end[start] = end

    #Need to make sure exons are sorted small to large before next step
    sorted_starts = sorted(start_to_end.keys())

    exon_counter=1
    for start in sorted_starts:
        end = start_to_end[start]
        new_start = start-1
        new_end = end
        #print("exon:",exon_counter," ",chr_seq[new_start:new_end]," start:",new_start," new end:",new_end )
        transcript_seq += str(chr_seq[start-1:end]) #I don't understand why this works to shift the start back one position, perhaps 1 base chromosome pos, zero based in python?
        exon_counter += 1

    bio_transcript_seq = Seq.Seq(transcript_seq)
    if strand == "-":   #reverse
        bio_transcript_seq = bio_transcript_seq.reverse_complement()

    #print("transcript1:",bio_transcript_seq)
    return bio_transcript_seq


### Argument 1: input_gff file
### 2: BioSeq object containing single chromosome
### 3: Locus ID to extract from GFF
### This version of the function is supplied not with full chromosome but just with small chunk, and gff coords need to be shifted by offset
def make_transcript_seq_with_offset(input_gff,bio_chr_seq,locus_id,offset):
    #print("Running method with offset: ", offset)
    transcript_seq = ""

    f= open(input_gff,"r")
    #print("opening:",input_gff)

    chr_seq = bio_chr_seq.seq
    start_to_end = {}

    strand = "+"
    for line in f:
        line = line[:-1]
        cells = line.split("\t")
        strand = cells[6]

        if cells[2].upper() == "CDS" and locus_id in line:
            start = int(cells[3])
            end = int(cells[4])
            #transcript_seq += str(chr_seq[start:end])
            start_to_end[start] = end

    #Need to make sure exons are sorted small to large before next step
    sorted_starts = sorted(start_to_end.keys())

    exon_counter = 1
    for start in sorted_starts:
        end = start_to_end[start]
        new_start = start-offset-1
        new_end = end-offset
        #print("exon:",exon_counter," ",chr_seq[new_start:new_end]," start:",new_start," new end:",new_end )
        transcript_seq += str(chr_seq[new_start:new_end])
        exon_counter += 1

    bio_transcript_seq = Seq.Seq(transcript_seq)
    if strand == "-":   #reverse
        bio_transcript_seq = bio_transcript_seq.reverse_complement()

    #print("transcript2:",bio_transcript_seq)
    return bio_transcript_seq




### Calls make_transcript_seq and translates, returning protein sequence
###4th argument says whether to use the offset version of the function
def make_protein_seq(input_gff,chr_seq,locus_id,offset):

    if offset:
        bio_prot_seq = make_transcript_seq_with_offset(input_gff, chr_seq, locus_id,offset).translate()
    else:
        bio_prot_seq = make_transcript_seq(input_gff,chr_seq,locus_id).translate()
    return bio_prot_seq

### This creates a new Bio sequence object, with SNPs substituted
### Argument 1: Bio Seq object containing one chromosome
### Argument 2: List of SNPs in this format 27504966:C;27504982:G;27505010:C
### Argument 3: Offset, in case giving the function a shortened chromosome for speed

def substitute_snps_on_chromosome(chr,snp_list,offset):

    #todo - this method is rather inefficient for processing many data points, since it slice chromosomes many times
    coding_dna = chr
    log = ""
    temp_snp_list = snp_list.split(";")
    for pos_and_base in temp_snp_list:
        temp = pos_and_base.split(":")
        try:
            pos = int(temp[0])
            base = temp[1]

            if base != "A" and base != "G" and base != "T" and base != "C":
                log += "Inconclusive SNP:" + pos_and_base + "\n"
            else:

                offset_pos = pos - offset
                if offset_pos >= 0:
                    #print("SNP pos:",pos," offset:",offset, "first string pos",str(pos-offset - 1), " till ", str(pos-offset) )
                    coding_dna = coding_dna[:offset_pos - 1] + base + coding_dna[offset_pos:]
        except:
            print("error extracting SNP pos from: ",pos_and_base )

    return_data = []
    return_data.append(coding_dna)
    #print("log:",log)
    return_data.append(log)
    return return_data


###Add a call to here just to run example tests
def run_tests():
    transcript_ID = "LOC_Os01g48060.1"
    chr_num = "1"

    ### Change these for testing
    #transscript_ID = "LOC_Os06g48950.1"
    #chr_num = "6"

    gene_ID = transcript_ID.split(".")[0]
    gff3 = gene_ID + ".gff3"
    chr_file_stub  = PATH_TO_DATABASE_FOLDER + "ChromosomeSplit/Oryza_sativa_ensembl.IRGSP-1.0.dna.toplevel"
    chr_file_name = chr_file_stub + chr_num + ".fasta"
    records = list(SeqIO.parse(chr_file_name, "fasta"))
    chr = records[0]  # Code only works for fasta file with single chromosome

    transcript = make_transcript_seq(PATH_TO_DATABASE_FOLDER + "cached_gff_chunks/"+gff3,chr,transcript_ID)
    print (transcript)

    transcripts = extract_list_of_transcript_seq(PATH_TO_DATABASE_FOLDER + "cached_gff_chunks/"+gff3,chr,transcript_ID)
    print("transcripts:",str(transcripts))

    protein = make_protein_seq(PATH_TO_DATABASE_FOLDER + "cached_gff_chunks/"+gff3,chr,transcript_ID)
    print(protein)


#run_tests()

#TGATGAAGCAGGCGCAGCAGCAGCCGCCGCCGCCACCGGCGAGCTCTGCGGCGACGACGACCACCGCGATGGCAGCCGCTGCGGCGGCGGCGGTGGTGGGGAGCGGGTGCGAAGGGGAGAAGACGAAGGCGCCGGCGATCAACTCGGAGCTGTGGCACGCCTGCGCGGGGCCGCTGGTGTCGCTGCCGCCGGCGGGCAGCCTCGTCGTCTACTTCCCCCAGGGCCACAGCGAGCAGGCGGACCCAGAAACAGATGAAGTGTATGCACAAATGACTCTTCAGCCAGTTACTTCAGATGGGAAGGAGGCCCTGCAGTTATCAGAGCTTGCACTCAAACAAGCGAGACCACAGACAGAATTCTTTTGCAAGACACTGACTGCAAGTGATACAAGTACTCATGGAGGCTTCTCTGTGCCTCGTCGAGCTGCAGAAAAGATATTTCCTCCACTGGACTTCTCAATGCAACCACCTGCACAAGAACTACAGGCCAGGGATTTGCATGATAATGTGTGGACATTCCGTCACATATATCGGGGTCAGCCAAAAAGGCATCTGCTTACCACTGGCTGGAGTCTATTTGTAAGCGGCAAGAGGTTATTTGCTGGAGATTCTGTCATTTTTGTCAGGGATGAAAAGCAGCAACTTCTATTAGGAATCAGGCGTGCTAACCGACAGCCAACTAACATATCATCATCTGTCCTTTCAAGTGACAGCATGCACATAGGGATTCTTGCTGCTGCAGCCCATGCTGCTGCCAACAATAGCCCATTTACCATCTTTTATAACCCTAGGGCCAGTCCTACTGAATTTGTTATCCCATTTGCTAAGTATCAGAAGGCAGTCTATGGTAATCAAATATCTTTAGGGATGCGCTTTCGCATGATGTTTGAGACTGAGGAATTAGGAACACGAAGGTACATGGGAACAATAACTGGCATAAGTGATCTAGATCCAGTAAGATGGAAAAACTCGCAGTGGCGCAACTTACAGGTTGGTTGGGATGAATCCGCAGCCGGTGAAAGGCGAAATAGGGTTTCTATCTGGGAGATTGAACCGGTCGCTGCTCCATTTTTCATATGTCCTCCACCATTTTTTGGTGCGAAGCGGCCCAGGCAATTAGGTGACGAGTCCTCGGAAATGGAGAATCTCTTAAAGAGGGCTATGCCTTGGCTTGGTGAGGAAATATGCATAAAGGATCCTCAGACTCAGAACACCATAATGCCTGGGCTGAGCTTGGTTCAGTGGATGAACATGAACATGCAACAGAGCTCCTCATTTGCGAATACAGCCATGCAGTCTGAGTACCTTCGATCATTGAGCAACCCCAACATGCAAAATCTTGGTGCCGCCGATCTCTCTAGGCAATTATGCCTGCAGAACCAGCTTCTTCAACAGAACAATATACAGTTTAATACTCCCAAACTTTCTCAGCAAATGCAGCCAGTCAATGAGTTAGCAAAGGCAGGCATTCCGTTGAATCAGCTTGGTGTGAGCACCAAACCTCAGGAACAGATTCATGATGCTAGCAACCTTCAGAGGCAACAACCTTCCATGAACCATATGCTTCCTTTGAGCCAAGCTCAAACCAATCTTGGCCAAGCTCAGGTCCTTGTCCAAAATCAAATGCAACAGCAACATGCATCTTCAACTCAAGGTCAACAACCAGCTACCAGCCAGCCCTTGCTTCTGCCCCAGCAGCAGCAACAGCAGCAGCAGCAGCAGCAACAACAACAACAACAGCAACAACAACAAAAATTGCTACAACAGCAGCAGCAACAGCTTTTGCTCCAGCAACAGCAGCAATTGAGTAAGATGCCTGCACAGTTGTCAAGTCTGGCGAATCAGCAGTTTCAGCTAACTGATCAACAGCTTCAGCTGCAACTGTTACAAAAACTACAGCAACAACAGCAGTCATTGCTTTCACAACCTGCAGTCACCCTTGCACAATTACCTCTGATCCAAGAACAGCAGAAGTTACTTCTGGATATGCAACAGCAGCTGTCAAACTCCCAAACACTTTCCCAACAACAAATGATGCCTCAACAAAGTACCAAGGTTCCATCACAGAACACACCATTGCCACTGCCTGTGCAACAAGAGCCACAACAGAAGCTTCTACAGAAGCAAGCGATGCTAGCAGACACTTCAGAAGCTGCCGTTCCGCCGACCACATCAGTCAATGTCATTTCAACAACTGGAAGCCCTTTGATGACAACTGGTGCTACTCATTCTGTACTTACAGAAGAAATCCCTTCTTGTTCAACATCACCATCCACAGCTAATGGCAATCACCTTCTACAACCAATACTTGGTAGGAACAAACATTGTAGCATGATCAACACAGAAAAGGTTCCTCAGTCTGCTGCTCCTATGTCAGTTCCAAGCTCCCTTGAAGCTGTCACAGCAACCCCGAGAATGATGAAGGATTCACCAAAGTTGAACCATAATGTTAAACAAAGTGTAGTGGCTTCAAAATTAGCAAATGCTGGGACTGGTTCTCAAAATTATGTGAACAATCCACCTCCAACGGACTATCTGGAAACTGCTTCTTCCGCAACTTCAGTGTGGCTTTCCCAGAATGATGGACTTCTACATCAAAATTTCCCTATGTCCAACTTCAACCAGCCACAGATGTTCAAAGATGCTCCTCCTGATGCTGAAATTCATGCTGCTAATACAAGTAACAATGCATTGTTTGGAATCAATGGTGATGGTCCGCTGGGCTTCCCTATAGGACTAGGAACAGATGATTTCCTGTCGAATGGAATTGATGCTGCCAAGTACGAGAACCATATCTCAACAGAAATTGATAATAGCTACAGAATTCCGAAGGATGCCCAGCAAGAAATATCATCCTCAATGGTTTCACAGTCATTTGGTGCATCAGATATGGCATTTAATTCAATTGATTCCACGATCAACGATGGTGGCTTTTTGAACCGGAGTTCTTGGCCTCCTGCCGCTCCCTTAAAGAGGATGAGGACATTCACCAAGGTATATAAGCGAGGAGCTGTAGGCCGGTCCATTGACATGAGTCAGTTCTCTGGATATGATGAATTAAAGCATGCTCTGGCACGGATGTTCAGTATAGAGGGGCAACTTGAGGAACGGCAGAGAATTGGTTGGAAGCTCGTTTACAAGGATCATGAAGATGACATCCTACTTCTTGGCGACGACCCATGGGAGGAATTTGTCGGTTGCGTGAAATGCATTAGGATCCTTTCACCTCAAGAAGTTCAGCAGATGAGCTTGGAGGGTTGTGATCTCGGGAACAACATTCCCCCGAATCAGGCCTGCAGCAGCTCAGACGGAGGGAATGCATGGAGGGCTCGCTGCGATCAGAACTCCGGGGCCATTCTTAAGATCTCCATGATGAAATCAAAAGTTGAAGATGTCAGGTATTGGAATACTGCGTAAT
#MVGIDLNTVEEEEDEEEGGATGTVTAPAEARAGGAVCLELWHACAGPVAPLPRKGSAVVYLPQGHLEHLGAAPGSGPGAAVPPHVFCRVVDVSLHADAATDEVYAQVSLVADNEEVERRMREGEDGAACDGEGEDAVKRPARIPHMFCKTLTASDTSTHGGFSVPRRAAEDCFPPLDYSLQRPFQELVAKDLHGTEWRFRHIYRGQPRRHLLTTGWSGFINKKKLVSGDAVLFLRGEDGELRLGVRRAAQLKNASPFPALHNQISNTSSLSEVAHAVAVKSIFHIYYNPSCTHRLSQSEFIIPYWKFMRSFSQPFSVGMRFKLRYESEDASERRRTGIIIGSREADPMWHGSKWKCLVVKWDDDVECRRPNGVSPWEIELSGSVSGSHLSTPHSKRLKSCFPQVNPDIVLPNGSVSSDFAESARFHKVLQGQELLGLKTRDGTVNTASQATEARNFQYTDERSCSINMSNNILGVPRLGVKTPSGNPGFSYHCSGFGESQRFQEVLQGQEVFRPYRGGTLSDACIRGSGFRQPDGNHAPGAAFKWLAPQGCDHHGITTSVLPQASSPSSVLMFPQTSSKMPGLEYIYGCLDRNENSRHFKIGPTQDMTRTDQTLRLWPHLISGKVLDECTRNEKLHSPVSGAEHESNNKCLNTNGCKIFGISLTEKAQAGDEVDCGNASYHSRLQSLKPQMPKSLGSSCATVHEQRPVVGRVVDISAVNTMI*
