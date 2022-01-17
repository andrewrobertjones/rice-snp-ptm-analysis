#Purpose of this script is to read protein MSA of 3K proteins for a given locus,
# and profile counts of each amino acid, per variety and grouped variety


from Bio import AlignIO
import os
import os, shutil, gzip

VAR_MAPPING_FILE = "3k_pop_details_UTF-8_v3.tsv"



var_to_original_var_name = {}

def gunzip_file(gzipped_file):

    input = gzip.GzipFile(gzipped_file, 'rb')
    s = input.read()
    input.close()

    gunzipped_output = os.path.splitext(gzipped_file)[0]
    print("Unzipping..." + gzipped_file)
    output = open(gunzipped_output, 'wb')
    output.write(s)
    output.close()
    return gunzipped_output


def read_alignment(alignment_gz_file):


    #align_handle = gzip.open(alignment_file, "rt")
    #for record in SeqIO.parse(handle, "fasta"):
    #    print(record.id, len(record))
    #handle.close()

    #alignment = AlignIO.read(alignment_file, "fasta")
    with gzip.open(alignment_gz_file, "rt") as fin:
        alignment = AlignIO.read(handle=fin, format="fasta")

    number_of_seqs = len(alignment)
    length_of_alignment = len(alignment[0].seq)

    #print("alen",length_of_alignment)
    #exit()

    amino_acid_subs = {}

    for j in range(0, length_of_alignment):
        chars_at_pos_to_varieties = {}
        for i in range(0,number_of_seqs):

            current_char=alignment[i].seq[j]
            varieties = []
            if current_char in chars_at_pos_to_varieties:
                varieties = chars_at_pos_to_varieties[current_char]

            #LOC_Os01g70270.1_LUBANG PUTI::IRGC 5429-1::[IRIS 313-10513] <unknown_description>
            #LOC_Os01g70270.1_QB_604::[CX210] <unknown description>
            #Os05t0392300-02_NIAO YAO::IRGC 5496-1::[IRIS 313-8743] <unknown description>
            description = alignment[i].description
            #if description[0:3] == "LOC":   #MSU
            #    variety_name = description.split("_",2)[2].split("<")[0][:-1].split("[")[0]    #grab LUBANG PUTI::IRGC 5429-1::
            #else:
            #    variety_name = description.split("_",1)[1].split("<")[0][:-1].split("[")[0]

            if description[0:3] == "LOC":  # MSU
                variety_name = description.split("_",2)[2].split(" ")[0]
            else:
                variety_name = description.split("_")[1].split(" ")[0]


            #Manual fix ("reference") also appears but can be safely ignored
            #if variety_name == "Reference":
            #    variety_name = "NIPPONBARE::"

            var_to_original_var_name[variety_name] = alignment[i].description
            varieties.append(variety_name)
            chars_at_pos_to_varieties[current_char] = varieties

        if len(chars_at_pos_to_varieties.keys()) > 1:
            #print("amino acid change found",repr(chars_at_pos_to_varieties))
            amino_acid_subs[j] = chars_at_pos_to_varieties
    return amino_acid_subs

def read_variety_group_mapping(mapping_file):
    variety_id_to_group = {}

    f = open(mapping_file,"r",encoding="Latin-1")

    for line in f:
        line = line[:-1]
        #print(line)
        cells = line.split("\t")
        #var_id = cells[0].split("[")[0] #PADI SIRANDAH KUNING::IRGC 73762-1::[IRIS 313-11904] take up to first [
        var_id = cells[1]   #Simpler to get the ID only
        var_group = cells[3]
        variety_id_to_group[var_id] = var_group

    return variety_id_to_group

#Count the different variety groups with each sub
def count_groups_per_sub(aa_subs,var_to_group,out_file):

    f_out = open(out_file,"w")
    f_out.write("pos,aa,vargroup,count\n")

    for pos in aa_subs:
        char_to_varieties = aa_subs[pos]
        for aa in char_to_varieties:
            vargroup_to_count = {}
            #print(repr(char_to_varieties[aa]))
            for variety in char_to_varieties[aa]:
                if variety in var_to_group:
                    group = var_to_group[variety]
                    count = 0

                    #print("group",group)
                    if group in vargroup_to_count:
                        count = vargroup_to_count[group]
                    count = count + 1
                    vargroup_to_count[group] = count  # increment counters
                else:
                    if variety != "reference" and variety != "Reference": #Safely ignore
                        print("variety not found",variety,"original name:",var_to_original_var_name[variety])


            for vargroup in vargroup_to_count:
                count = vargroup_to_count[vargroup]
                #print(pos,aa,vargroup,count)
                f_out.write(str(pos+1)+","+aa+","+vargroup+","+str(count)+"\n") #Positions changed to start counting at 1
            #print(repr(vargroup_to_count))
    print("Finished writing to ",out_file)

def create_MSA_pos_stats(three_k_prot_align_file,out_folder):
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    amino_acid_subs_vars = read_alignment(three_k_prot_align_file)
    var_id_to_group = read_variety_group_mapping(VAR_MAPPING_FILE)
    output_file = out_folder + os.path.basename(three_k_prot_align_file).split(".")[0] + "_stats.csv"
    count_groups_per_sub(amino_acid_subs_vars, var_id_to_group, output_file)


#ALIGNMENT_FILE = "variety_fasta/ARFs/LOC_Os01g70270.1_proteins_in_varieties.fasta"
def run_tests():
    ALIGNMENT_FILE = "variety_fasta/ARFs/LOC_Os11g32110.1_proteins_in_varieties.fasta"
    amino_acid_subs_vars = read_alignment(ALIGNMENT_FILE)
    var_id_to_group = read_variety_group_mapping(VAR_MAPPING_FILE)

    output_file = "temp_alignments_outs/" + os.path.basename(ALIGNMENT_FILE).split(".")[0] + "_stats.csv"
    count_groups_per_sub(amino_acid_subs_vars,var_id_to_group,output_file)