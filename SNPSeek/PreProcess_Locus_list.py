#Script to read in a list of MSU identifiers:
#Each row might have
# and check for those that are flagged as canoni

PATH_TO_MSU_LOCUS_FILE = "D:/Dropbox/DocStore/ProteomicsSoftware/PTMExchange/Rice_build/MSU_IDs_for_phospho.txt"
PATH_TO_MSU_OUTPUT_FILE = "D:/Dropbox/DocStore/ProteomicsSoftware/PTMExchange/Rice_build/MSU_IDs_for_phospho_cleaned.txt"
PATH_TO_NOT_FOUND_LOCUS_FILE = "D:/Dropbox/DocStore/ProteomicsSoftware/PTMExchange/Rice_build/MSU_IDs_found_but_not_canonical.txt"

f = open("all.locus_brief_info.7.0","r")

locus_to_is_canonical = {}
for line in f:
    line.rstrip("\n")
    cells = line.split("\t")
    locus_ID = cells[2]
    is_canonical = cells[8]
    if is_canonical == "Y":
        locus_to_is_canonical[locus_ID] = is_canonical


f.close()

#take the canonical one from each list, remove suffix
#LOC_Os01g01040.1_511;LOC_Os01g01040.3_511;LOC_Os01g01040.2_511
#LOC_Os01g01080.1_80
#LOC_Os01g01080.1_81
#LOC_Os01g01080.2_350;LOC_Os01g01080.1_501;LOC_Os01g01080.3_350
#LOC_Os01g01130.1_41
#LOC_Os01g01150.2_36;LOC_Os01g01150.1_36;LOC_Os01g01150.3_36
#LOC_Os01g01150.2_512;LOC_Os01g01150.1_512;LOC_Os01g01150.3_512
#LOC_Os01g01150.2_540;LOC_Os01g01150.1_540;LOC_Os01g01150.3_540
#LOC_Os01g01150.2_624;LOC_Os01g01150.1_624;LOC_Os01g01150.3_624
#LOC_Os01g01150.2_626;LOC_Os01g01150.1_626;LOC_Os01g01150.3_626
#LOC_Os01g01150.2_63;LOC_Os01g01150.1_63;LOC_Os01g01150.3_63
#LOC_Os01g01150.2_728;LOC_Os01g01150.1_728;LOC_Os01g01150.3_728


f = open(PATH_TO_MSU_LOCUS_FILE,"r")

found_canonical_MSU_IDs = {}
rows_without_match = {}

for line in f:
    msu_list = line.rstrip("\n")
    msu_ids = msu_list.split(";")

    canonical_found = False
    clean_snps = ""
    for msu_id in msu_ids:
        split_id = msu_id.split("_")
        clean_id = split_id[0] + "_" + split_id[1]  #reconstruct MSU_ID without the suffix (36) LOC_Os01g01150.3_36
        clean_snps += clean_id + ";"

        if clean_id in locus_to_is_canonical:
            canonical_found = True
            if clean_id not in found_canonical_MSU_IDs:
                found_canonical_MSU_IDs[clean_id] = 1   #Add once to output dictionary

    if canonical_found == False:
        clean_snps = clean_snps[:-1] #remove final ;
        rows_without_match[clean_snps] = 1
        print("row has no canonical ID:",clean_snps)

f.close()

f_out = open(PATH_TO_MSU_OUTPUT_FILE,"w")

for msu_id in found_canonical_MSU_IDs:
    f_out.write(msu_id + "\n")
f_out.close()

f_out = open(PATH_TO_NOT_FOUND_LOCUS_FILE,"w")
for msu_id in rows_without_match:
    f_out.write(msu_id + "\n")
f_out.close()

