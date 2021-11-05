from Bio import SeqIO

from SNPSeek.Get_variety_snps_from_csv import *
from HelperRoutines.CreateTranscriptOrProteinFromGff import *
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from SNPSeek.RetrieveJsonForLocus import retrieve_json
from SNPSeek.ConvertOneJsonToCsv import convert_json_text
import os
from Extract_stats_from_protein_alignment import create_MSA_pos_stats
import os, shutil, gzip
#import _thread #needed for multi-threading
from threading import Thread

PATH_TO_DATABASE_FOLDER = "../../RiceDatabases/"

def process_one_variety(snp_string,short_chr,trans_id,offset,gff_file):
    #snp_string_values = "27501530:G;27501902:T;27501903:C;27501931:T;27501941:G;27501942:A;27501966:G;27501975:T;27501977:T;27501994:C;27502082:A;27502103:C;27502128:A;27502182:A;27502289:G;27502324:C;27502379:G;27502388:T;27502420:C;27502598:C;27502627:A;27502647:G;27502886:T;27502940:T;27503019:T;27503025:A;27503028:G;27503049:T;27503125:A;27503210:A;27503309:C;27503381:A;27503424:T;27503450:G;27503457:A;27503519:C;27503544:C;27503565:T;27503580:T;27503585:A;27503619:G;27503624:C;27503668:A;27503703:T;27503712:A;27503755:A;27503762:G;27503779:T;27503811:G;27503829:A;27503833:G;27503871:A;27503914:G;27503962:T;27503965:T;27503977:A;27503984:T;27503987:A;27504011:C;27504021:A;27504027:G;27504030:A;27504038:G;27504045:C;27504096:T;27504210:A;27504216:G;27504251:C;27504290:A;27504307:C;27504421:A;27504431:C;27504440:T;27504443:A;27504464:G;27504567:G;27504639:T;27504649:G;27504663:G;27504683:T;27504715:A;27504716:G;27504769:G;27504772:T;27504775:C;27504779:T;27504784:T;27504830:C;27504883:A;27504896:G;27504899:G;27504903:T;27504904:G;27504907:A;27504940:G;27504942:T;27504966:A;27504982:C;27505010:T;27505029:C;27505050:C;27505064:A;27505089:C;27505112:T;27505130:A;27505150:A;27505268:T;27505304:A;27505394:G;27505399:T;27505442:A;27505445:C;27505614:T;27505618:G;27505682:C;27505724:C;27505776:T;27505787:A;27505892:A;27506062:T;27506286:G;27506486:C;27506575:A;27506578:G;27506779:G;27506792:A;27506796:C;27506798:G;27506802:G"

    # transcript = make_transcript_seq("../Databases/cached_gff_chunks/"+gff3,chr,transscript_ID)
    # print (transcript)
    # protein = make_protein_seq("../../Databases/cached_gff_chunks/"+gff3,chr,transscript_ID)

    sequence_and_log = substitute_snps_on_chromosome(short_chr, snp_string,offset)
    changed_sequence = sequence_and_log[0]
    log = sequence_and_log[1]

    #sub_seq = sequence.seq
    #print("seq", sub_seq[27502003:27502068])
    #substituted_protein = make_protein_seq("../../Databases/cached_gff_chunks/" + gff_file, changed_sequence, trans_id,offset)
    substituted_protein = make_protein_seq(PATH_TO_DATABASE_FOLDER + "cached_gff_chunks/" + gff_file, changed_sequence, trans_id,offset)
    #print(protein2)
    return substituted_protein


def get_data_from_snpseek(transcript_ID):

    success = False
    if transcript_ID[4] == "t":
        #Os01t0614500-0  to Os01g0614500
        #print("Need to convert to gene ID for retrieval", transcript_ID)
        cells = transcript_ID.replace("t","g").split("-")
        gene_ID = cells[0] + "." + str(int(cells[1]))
        #print("Need to convert to gene ID for retrieval", transcript_ID, "geneID", gene_ID )
        transcript_ID = gene_ID


    returned_data = retrieve_json(transcript_ID)    #uses locus for SNPs but transcript ID for other things, this method converts and returns
    locus = returned_data[0]
    json_retrieved=returned_data[1]
    if json_retrieved:
        found_snps = convert_json_text(locus)

        if found_snps:

            #gene_ID = transcript_ID.split(".")[0]
            gff3 = ""
            if locus[0:3] == "LOC":
                #gff3 = transcript_ID + ".gff"
                success = True
            elif locus[0:2] == "Os" and locus[4] == "g":#i.e. gene locus provided, convert (back) to transcript
                transcript_ID = transcript_ID[0:4] + "t" + transcript_ID[5:].split(".")[0] + "-0" + transcript_ID[5:].split(".")[1] #Convert gene ID Os01g0625300.1 to t_Id: Os01t0625300-01
                #gff3 = transcript_ID + ".gff"
                success = True
            else:
                print("transcript ID not recognised yet, need to add code for BGI and n22",transcript_ID)
        else:
            print("No SNPs found for",locus," nothing to do")
    else:
        print("Cannot proceed with extraction of data, no response from SNP-Seek")

    if success:
        return(locus,transcript_ID)
    else:
        return(False)

def create_alignment(locus,transcript_id,chr_num,output_location):
    print("\tAbout to create alignment for",transcript_id)
    chr_file_stub = PATH_TO_DATABASE_FOLDER + "ChromosomeSplit/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa"
    chr_file_name = chr_file_stub + chr_num + ".fasta"  # Assumes already run / cached the chromosome splitter
    records = list(SeqIO.parse(chr_file_name, "fasta"))
    chr = records[0]  # Code only works for fasta file with single chromosome

    chr_seq = chr.seq

    gff3 = transcript_id + ".gff"
    coding_start_and_end = get_start_end_positions_from_gff(PATH_TO_DATABASE_FOLDER + "cached_gff_chunks/" + gff3, transcript_id)

    coding_start = coding_start_and_end[0] - 1  # Need to grab base before, based on gff vs python numbering system
    coding_end = coding_start_and_end[1]

    short_chromosome = SeqRecord(chr_seq[coding_start:coding_end])
    ref_protein = make_protein_seq(PATH_TO_DATABASE_FOLDER + "cached_gff_chunks/" + gff3, short_chromosome, transcript_id,
                                   coding_start)

    reference_record = SeqRecord(ref_protein, id=transcript_id + "_reference")
    input_text = PATH_TO_DATABASE_FOLDER + "temp_csv_files/" + locus + ".tsv"

    var_to_snp_dict = get_variety_to_snps_dictionary(input_text)

    if not os.path.exists(output_location):
        os.makedirs(output_location)
    output_fasta = output_location + transcript_id + "_proteins_in_varieties.fasta"
    output_fasta_diffs = output_location + transcript_id + "_different_proteins_in_varieties.fasta"

    out_records = []
    out_records.append(reference_record)

    out_diff_records = []
    out_diff_records.append(reference_record)

    counter = 1
    for variety in var_to_snp_dict:
        snps = var_to_snp_dict[variety]
        # print(variety," ->" ,snps )
        if counter < 5:
            variety_protein = process_one_variety(snps, short_chromosome, transcript_id, coding_start, gff3)
            # variety_protein = process_one_variety(snps, chr, transcript_ID, 0, gff3)
            variety_record = SeqRecord(variety_protein, id=transcript_id + "_" + variety)
            out_records.append(variety_record)
            if str(variety_protein) != str(reference_record.seq):
                out_diff_records.append(variety_record)
            # else:
            # print("match",variety,"and reference")
        # counter += 1   #Uncomment for testing

    out_gz = output_fasta + ".gz"
    #gzf = gzip.open(out_gz, "wb")

    with gzip.open(out_gz, "wt") as fout:
        SeqIO.write(sequences=out_records, handle=fout, format="fasta")

    out_gz2 = output_fasta_diffs + ".gz"    #records just showing which are different
    with gzip.open(out_gz2, "wt") as fout2:
        SeqIO.write(sequences=out_diff_records, handle=fout2, format="fasta")

    #out_records.close()
    #SeqIO.write(out_diff_records, output_fasta_diffs, "fasta")

    #out_gz = output_fasta + ".gz"
    #with open(output_fasta, "rb") as fin, gzip.open(out_gz, "wb") as fout:
    #    # Reads the file by chunks to avoid exhausting memory
    #    shutil.copyfileobj(fin, fout)
    #fin.close()
    #fout.close()

    #in_data = open(output_fasta , "rb").read()
    #out_gz = output_fasta + ".gz"
    #gzf = gzip.open(out_gz, "wb")
    #gzf.write(in_data)
    #gzf.close()
    #os.unlink(output_fasta) #remove unzipped version

    # If you want to delete the original file after the gzip is done:
    #os.unlink(in_file)


    return out_gz

#snp_string = read_csv_convert_to_string("temp_csv_files/LOC_Os01g48060.csv","COLOMBIA")

#Testing changing first ATG to CTG
#snp_string =  "27507024:G;27501530:G;27501902:T;27501903:C;27501931:T;27501941:G;27501942:A;27501966:G;27501975:T;27501977:T;27501994:C;27502082:A;27502103:C;27502128:A;27502182:A;27502289:G;27502324:C;27502379:G;27502388:T;27502420:C;27502598:C;27502627:A;27502647:G;27502886:T;27502940:T;27503019:T;27503025:A;27503028:G;27503049:T;27503125:A;27503210:A;27503309:C;27503381:A;27503424:T;27503450:G;27503457:A;27503519:C;27503544:C;27503565:T;27503580:T;27503585:A;27503619:G;27503624:C;27503668:A;27503703:T;27503712:A;27503755:A;27503762:G;27503779:T;27503811:G;27503829:A;27503833:G;27503871:A;27503914:G;27503962:T;27503965:T;27503977:A;27503984:T;27503987:A;27504011:C;27504021:A;27504027:G;27504030:A;27504038:G;27504045:C;27504096:T;27504210:A;27504216:G;27504251:C;27504290:A;27504307:C;27504421:A;27504431:C;27504440:T;27504443:A;27504464:G;27504567:G;27504639:T;27504649:G;27504663:G;27504683:T;27504715:A;27504716:G;27504769:G;27504772:T;27504775:C;27504779:T;27504784:T;27504830:C;27504883:A;27504896:G;27504899:G;27504903:T;27504904:G;27504907:A;27504940:G;27504942:T;27504966:A;27504982:C;27505010:T;27505029:C;27505050:C;27505064:A;27505089:C;27505112:T;27505130:A;27505150:A;27505268:T;27505304:A;27505394:G;27505399:T;27505442:A;27505445:C;27505614:T;27505618:G;27505682:C;27505724:C;27505776:T;27505787:A;27505892:A;27506062:T;27506286:G;27506486:C;27506575:A;27506578:G;27506779:G;27506792:A;27506796:C;27506798:G;27506802:G"

def process_multiple_loci(file_of_loci,output_location):

    loci_to_process_file =  open(file_of_loci,"r")

    counter = 1
    for line in loci_to_process_file:
        line = line[:-1]
        cells = line.split("\t")
        t_ID = cells[0]

        if "Os" in t_ID:
            #RAP example Os01t0588500-01
            #MSU example LOC_Os08g05540.1
            temp_transcript_ID = t_ID.replace("LOC_","")
            chromosome = str(int(temp_transcript_ID[2:4])) #get chromosomal position, convert to int then back to string to remove trailing zeroes
            returned_data = get_data_from_snpseek(t_ID)
            if returned_data: #False if failed
                locus = returned_data[0]
                transcript_id = returned_data[1]
                main_msa_fasta_gz = create_alignment(locus,transcript_id,chromosome,output_location)
                create_MSA_pos_stats(main_msa_fasta_gz,output_location)

        else:
            print("transcript not recognised, will not be processed (only RAP-DB and MSU IDs supporteD)",t_ID)

def process_multiple_loci_threaded(file_of_loci,output_location):

    loci_to_process_file =  open(file_of_loci,"r")

    counter = 1
    transcript_list = []
    for line in loci_to_process_file:
        line = line[:-1]
        cells = line.split("\t")
        t_ID = cells[0]
        transcript_list.append(t_ID)

    run_multi(transcript_list,output_location)

def run_multi(transcript_list,output_location):
    try:

        threads = []

        total_transcripts = len(transcript_list)

        for n in range(0, total_transcripts):
            transcript_ID = transcript_list[n]
            print("Starting thread",str(n+1), "for transcript",transcript_ID)
            t = Thread(target=process_one_transcript, args=(transcript_ID,output_location))
            threads.append(t)
            t.start()

        for t in threads:
            t.join()


        #thread_counter = 0
        #for transcript_ID in transcript_list:
        #    print("starting thread for ",transcript_ID)
        #    #_thread.start_new_thread(process_one_transcript, (transcript_ID,)) #the comma forces this to be a tuple
        #    new_thread = Thread(target=process_one_transcript, args=(transcript_ID,))
        #    new_thread.start()

    except:
        print("Error starting thread")



def process_one_transcript(transcript_ID,out_folder):

    if "Os" in transcript_ID:
        # RAP example Os01t0588500-01
        # MSU example LOC_Os08g05540.1
        temp_transcript_ID = transcript_ID.replace("LOC_", "")
        chromosome = str(int(temp_transcript_ID[
                             2:4]))  # get chromosomal position, convert to int then back to string to remove trailing zeroes

        returned_data = get_data_from_snpseek(transcript_ID)
        if returned_data:  # False if failed
            locus = returned_data[0]
            transcript_id = returned_data[1]
            main_msa_fasta_gz = create_alignment(locus, transcript_id, chromosome, out_folder)
            create_MSA_pos_stats(main_msa_fasta_gz, out_folder)
    else:
        print("transcript not recognised, will not be processed (only RAP-DB and MSU IDs supporteD)", transcript_ID)

def one_locus_test():
    #t_ID = "LOC_Os01g48060.1"
    t_ID = "Os01g0625300.1"
    t_ID = "LOC_Os05g45410.1"
    #t_ID = "LOC_Os02g35140.1"
    #t_ID = "LOC_Os04g43910.1"
    #chr = "4"

    #t_ID = "LOC_Os12g41380.1"   #SUMO protease
    chr = "5"

    get_data_from_snpseek(t_ID,chr)
    #transcript_ID = "LOC_Os06g09660.1"
    #chr_num = "6"


#one_locus_test()
#process_multiple_loci("ProcessHSF_Xinyang_project_output.txt")
#process_multiple_loci("ProcessARFs.txt")

#File format is one locus per line, can be MSU or RAP-DB gene or transcript IDs (although better tested with RAP-DB transcript IDs)
#data_loc = "D:/Dropbox/DocStore/ProteomicsSoftware/PTMExchange/Rice_build/Rice CSV files/Example_Rice_Modified_Proteins.txt"
#data_loc = "D:/Dropbox/DocStore/ProteomicsSoftware/PTMExchange/Rice_build/Rice CSV files/three_prots.txt"
#data_loc = "ProcessHSF_Xinyang_project_output.txt"
#output_loc = "D:/Dropbox/DocStore/ProteomicsSoftware/PTMExchange/Rice_build/Rice CSV files/out_saaps/"
#process_multiple_loci(data_loc,output_loc)


f = open("run_details.txt", "r")
data_loc = f.readline()[:-1]
output_loc = f.readline()[:-1]
process_multiple_loci_threaded(data_loc,output_loc)


