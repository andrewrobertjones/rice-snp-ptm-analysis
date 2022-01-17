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
PATH_TO_SNP_FOLDER = "../../RiceDatabases/Genes-3K-Base-CSV-20211130T093248Z-001/Genes-3K-Base-CSV/"


def process_one_variety(snp_string,short_chr,trans_id,offset,gff_file):
    #snp_string_values = "8_2968783:G;8_2968844:C;8_2969033:C',

    # transcript = make_transcript_seq("../Databases/cached_gff_chunks/"+gff3,chr,transscript_ID)
    # print (transcript)
    # protein = make_protein_seq("../../Databases/cached_gff_chunks/"+gff3,chr,transscript_ID)

    sequence_and_log = substitute_snps_on_chromosome(short_chr, snp_string,offset)
    changed_sequence = sequence_and_log[0]
    log = sequence_and_log[1]

    #sub_seq = sequence.seq
    #print("seq", sub_seq[27502003:27502068])
    #substituted_protein = make_protein_seq("../../Databases/cached_gff_chunks/" + gff_file, changed_sequence, trans_id,offset)
    #TODO - different csv format for local run compared to downloaded
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

def create_alignment(locus,transcript_id,chr_num,output_location,is_local_csv_files):
    print("\tAbout to create alignment for",transcript_id)
    out_gz = None
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


    if is_local_csv_files:
        #input_text = PATH_TO_DATABASE_FOLDER + "msu-Chr1-3kFilt/genes-MSU-filtered-3K/" + locus + "/" + locus+".csv"

        input_text = PATH_TO_SNP_FOLDER + "MSUgene-Base3K-chr"+chr_num+"/"+chr_num+"/"+ locus + "/" + transcript_id + ".csv"
        var_to_snp_dict = get_variety_to_snps_dictionary(input_text)
    else:
        input_text = PATH_TO_DATABASE_FOLDER + "temp_csv_files/" + locus + ".csv"
        var_to_snp_dict = get_variety_to_snps_dictionary(input_text)

    if len(var_to_snp_dict) > 1:    #Correctly processed, otherwise empty
        if not os.path.exists(output_location):
            os.makedirs(output_location)
        output_fasta = output_location + transcript_id + "_proteins_in_varieties.fasta"
        output_fasta_diffs = output_location + transcript_id + "_different_proteins_in_varieties.fasta"

        out_records = []
        out_records.append(reference_record)

        out_diff_records = []
        out_diff_records.append(reference_record)

        counter = 1

        cached_snps_prot = {}   #Only create new protein sequence if SNPs have changed

        for variety in var_to_snp_dict:
            snps = var_to_snp_dict[variety]
            # print(variety," ->" ,snps )
            if counter < 5:

                if snps not in cached_snps_prot:
                    variety_protein = process_one_variety(snps, short_chromosome, transcript_id, coding_start, gff3)
                    cached_snps_prot[snps] = variety_protein
                else:
                    variety_protein = cached_snps_prot[snps]

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
    else:
        print("\tFile not processed")

    return out_gz

#snp_string = read_csv_convert_to_string("temp_csv_files/LOC_Os01g48060.csv","COLOMBIA")

#Arguments are:
# 1) input file of transcript IDs to process,
# 2) location for output to be written;
# 3) Boolean of whether to do a local run when true, (files already in a local file store) or query SNP-Seq for Json download (false)
def process_multiple_loci_threaded(file_of_loci,output_location,do_local_run):

    loci_to_process_file =  open(file_of_loci,"r")

    counter = 1
    transcript_list = []
    for line in loci_to_process_file:
        line = line.rstrip('\n')
        cells = line.split("\t")
        t_ID = cells[0]
        transcript_list.append(t_ID)

    run_multi(transcript_list,output_location,do_local_run)

#Arguments are:
# 1) input file of transcript IDs to process,
# 2) location for output to be written;
# 3) Boolean of whether to do a local run when true, (files already in a local file store) or query SNP-Seq for Json download (false)
def run_multi(transcript_list,output_location,do_local_run):
    try:

        threads = []

        total_transcripts = len(transcript_list)

        for n in range(0, total_transcripts):
            transcript_ID = transcript_list[n]
            print("Starting thread",str(n+1), "for transcript",transcript_ID)
            t = Thread(target=process_one_transcript, args=(transcript_ID,output_location,do_local_run))
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


#Arguments are:
# 1) input file of transcript IDs to process,
# 2) location for output to be written;
# 3) Boolean of whether to do a local run when true, (files already in a local file store) or query SNP-Seq for Json download (false)
def process_one_transcript(transcript_ID,out_folder,do_local_run):

    if "Os" in transcript_ID:
        # RAP example Os01t0588500-01
        # MSU example LOC_Os08g05540.1
        temp_transcript_ID = transcript_ID.replace("LOC_", "")
        chromosome = str(int(temp_transcript_ID[
                             2:4]))  # get chromosomal position, convert to int then back to string to remove trailing zeroes

        locus_ID = transcript_ID
        if transcript_ID.find("."):  # Need to remove suffix for SNP-SEEK   #TODO does this work for RAP-DB transcript IDs?
            locus_ID = transcript_ID[:transcript_ID.find(".")]

        if do_local_run == False:
            returned_data = get_data_from_snpseek(transcript_ID)

        main_msa_fasta_gz = create_alignment(locus_ID, transcript_ID, chromosome, out_folder,do_local_run)
        if main_msa_fasta_gz != None:   #Error handling for case that csv is not present
            create_MSA_pos_stats(main_msa_fasta_gz, out_folder)
    else:
        print("transcript not recognised, will not be processed (only RAP-DB and MSU IDs supported)", transcript_ID)

def one_locus_test():
    #t_ID = "LOC_Os01g48060.1"
    #t_ID = "Os01g0625300.1"
    #t_ID = "LOC_Os05g45410.1"
    t_ID = "LOC_Os01g01100.1"
    #t_ID = "LOC_Os02g35140.1"
    #t_ID = "LOC_Os04g43910.1"
    #chr = "4"

    #t_ID = "LOC_Os12g41380.1"   #SUMO protease
    chr = "5"

    get_data_from_snpseek(t_ID)
    returned_data = get_data_from_snpseek(t_ID)
    temp_transcript_ID = t_ID.replace("LOC_", "")
    chromosome = str(int(temp_transcript_ID[
                         2:4]))  # get chromosomal position, convert to int then back to string to remove trailing zeroes

    output_location = "D:/Dropbox/DocStore/ProteomicsSoftware/PTMExchange/Rice_build/Rice CSV files/out_saaps/"
    if returned_data:  # False if failed
        locus = returned_data[0]
        transcript_id = returned_data[1]
        main_msa_fasta_gz = create_alignment(locus, transcript_id, chromosome, output_location,False)
        create_MSA_pos_stats(main_msa_fasta_gz, output_location)
    #transcript_ID = "LOC_Os06g09660.1"
    #chr_num = "6"


#process_multiple_loci("ProcessHSF_Xinyang_project_output.txt")
#process_multiple_loci("ProcessARFs.txt")

#File format is one locus per line, can be MSU or RAP-DB gene or transcript IDs (although better tested with RAP-DB transcript IDs)
#data_loc = "D:/Dropbox/DocStore/ProteomicsSoftware/PTMExchange/Rice_build/Rice CSV files/Example_Rice_Modified_Proteins.txt"
#data_loc = "D:/Dropbox/DocStore/ProteomicsSoftware/PTMExchange/Rice_build/Rice CSV files/three_prots.txt"
#data_loc = "ProcessHSF_Xinyang_project_output.txt"
#output_loc = "D:/Dropbox/DocStore/ProteomicsSoftware/PTMExchange/Rice_build/Rice CSV files/out_saaps/"
#process_multiple_loci(data_loc,output_loc)


f = open("run_details.txt", "r")
data_loc = f.readline().rstrip("\n")
output_loc = f.readline().rstrip("\n")
do_local_run = True
process_multiple_loci_threaded(data_loc,output_loc,do_local_run)


