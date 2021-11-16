import json
import os
import sys

import csv
import numpy as np

PATH_TO_DATABASE_FOLDER = "../../RiceDatabases/"
IN_FOLDER = PATH_TO_DATABASE_FOLDER + "cached_json_files/"     #location of JSON files cached from SNP-Seek
#snp_dictionary = {}


#Return the chromosome from the position dictionary entry
def get_chromosome(snp_dictionary):
    chr = None
    ref_positions = snp_dictionary["position"]
    if ref_positions != None and len(ref_positions)>0:
        first_position = ref_positions[0] #This returns another dictionary, inside each list cell
        chr = first_position["chromosome"]
    else:
        print("No SNPs found")
    return(chr)

#Extract positions from the dictionary entry
#key: position  value: [{'snpFeatureId': 3133832, 'typeId': 14, 'chromosome': 1, 'position': 40695776, 'refcall': 'C', 'alleleIndex': 540787, 'contig': 'chr01', 'refnuc': 'C', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133833, 'typeId': 14, 'chromosome': 1, 'position': 40695781, 'refcall': 'C', 'alleleIndex': 540788, 'contig': 'chr01', 'refnuc': 'C', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133834, 'typeId': 14, 'chromosome': 1, 'position': 40695844, 'refcall': 'G', 'alleleIndex': 540789, 'contig': 'chr01', 'refnuc': 'G', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133843, 'typeId': 14, 'chromosome': 1, 'position': 40696015, 'refcall': 'A', 'alleleIndex': 540790, 'contig': 'chr01', 'refnuc': 'A', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133847, 'typeId': 14, 'chromosome': 1, 'position': 40696206, 'refcall': 'G', 'alleleIndex': 540791, 'contig': 'chr01', 'refnuc': 'G', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133849, 'typeId': 14, 'chromosome': 1, 'position': 40696256, 'refcall': 'G', 'alleleIndex': 540792, 'contig': 'chr01', 'refnuc': 'G', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133863, 'typeId': 14, 'chromosome': 1, 'position': 40696365, 'refcall': 'C', 'alleleIndex': 540793, 'contig': 'chr01', 'refnuc': 'C', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133866, 'typeId': 14, 'chromosome': 1, 'position': 40696449, 'refcall': 'C', 'alleleIndex': 540794, 'contig': 'chr01', 'refnuc': 'C', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133870, 'typeId': 14, 'chromosome': 1, 'position': 40696600, 'refcall': 'C', 'alleleIndex': 540795, 'contig': 'chr01', 'refnuc': 'C', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133871, 'typeId': 14, 'chromosome': 1, 'position': 40696639, 'refcall': 'A', 'alleleIndex': 540796, 'contig': 'chr01', 'refnuc': 'A', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133907, 'typeId': 14, 'chromosome': 1, 'position': 40697439, 'refcall': 'C', 'alleleIndex': 540797, 'contig': 'chr01', 'refnuc': 'C', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133930, 'typeId': 14, 'chromosome': 1, 'position': 40697621, 'refcall': 'A', 'alleleIndex': 540798, 'contig': 'chr01', 'refnuc': 'A', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133936, 'typeId': 14, 'chromosome': 1, 'position': 40697717, 'refcall': 'C', 'alleleIndex': 540799, 'contig': 'chr01', 'refnuc': 'C', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133947, 'typeId': 14, 'chromosome': 1, 'position': 40697886, 'refcall': 'C', 'alleleIndex': 540800, 'contig': 'chr01', 'refnuc': 'C', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133959, 'typeId': 14, 'chromosome': 1, 'position': 40698335, 'refcall': 'A', 'alleleIndex': 540801, 'contig': 'chr01', 'refnuc': 'A', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133975, 'typeId': 14, 'chromosome': 1, 'position': 40698889, 'refcall': 'G', 'alleleIndex': 540802, 'contig': 'chr01', 'refnuc': 'G', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133976, 'typeId': 14, 'chromosome': 1, 'position': 40698908, 'refcall': 'C', 'alleleIndex': 540803, 'contig': 'chr01', 'refnuc': 'C', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133986, 'typeId': 14, 'chromosome': 1, 'position': 40699328, 'refcall': 'G', 'alleleIndex': 540804, 'contig': 'chr01', 'refnuc': 'G', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133989, 'typeId': 14, 'chromosome': 1, 'position': 40699488, 'refcall': 'G', 'alleleIndex': 540805, 'contig': 'chr01', 'refnuc': 'G', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133991, 'typeId': 14, 'chromosome': 1, 'position': 40699530, 'refcall': 'T', 'alleleIndex': 540806, 'contig': 'chr01', 'refnuc': 'T', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133993, 'typeId': 14, 'chromosome': 1, 'position': 40699668, 'refcall': 'A', 'alleleIndex': 540807, 'contig': 'chr01', 'refnuc': 'A', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3133994, 'typeId': 14, 'chromosome': 1, 'position': 40699683, 'refcall': 'C', 'alleleIndex': 540808, 'contig': 'chr01', 'refnuc': 'C', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3134014, 'typeId': 14, 'chromosome': 1, 'position': 40700055, 'refcall': 'C', 'alleleIndex': 540809, 'contig': 'chr01', 'refnuc': 'C', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3134015, 'typeId': 14, 'chromosome': 1, 'position': 40700140, 'refcall': 'G', 'alleleIndex': 540810, 'contig': 'chr01', 'refnuc': 'G', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3134019, 'typeId': 14, 'chromosome': 1, 'position': 40700351, 'refcall': 'T', 'alleleIndex': 540811, 'contig': 'chr01', 'refnuc': 'T', 'altnuc': '', 'chr': 1}, {'snpFeatureId': 3134021, 'typeId': 14, 'chromosome': 1, 'position': 40700499, 'refcall': 'T', 'alleleIndex': 540812, 'contig': 'chr01', 'refnuc': 'T', 'altnuc': '', 'chr': 1}]
def get_position_on_chr_as_list(snp_dictionary):
    snp_positions = []
    ref_positions = snp_dictionary["position"]
    for position_data in ref_positions: #position_data is dictionary like this {'snpFeatureId': 3133832, 'typeId': 14, 'chromosome': 1, 'position': 40695776, 'refcall': 'C', 'alleleIndex': 540787, 'contig': 'chr01', 'refnuc': 'C', 'altnuc': '', 'chr': 1}
        snp_pos=position_data["position"]
        snp_positions.append(snp_pos) #add just the position to the new dictionary
    return snp_positions

#Get reference bases, key: reference  value: ['C', 'C', 'G', 'A', 'G', 'G', 'C', 'C', 'C', 'A', 'C', 'A', 'C', 'C', 'A', 'G', 'C', 'G', 'G', 'T', 'A', 'C', 'C', 'G', 'T', 'T']
def get_reference_bases(snp_dictionary):
#    ref_bases = []
    ref_positions = snp_dictionary["position"]
    first_position = ref_positions[0] #This returns another dictionary, inside each list cell
    chr = first_position["chromosome"]
    return(chr)

#This function is going to loop through the 2D list - outer list is 3116 varieties, inner list is the count of positions for this locus
#The function will return a dictionary, where the keys are the variety names, the values of each dictionary entry is
def get_variety_alleles(snp_dictionary):
    var_alleles = {} #Create a dictionary where the key is the variety name, and the value is a list of alleles for each position (of length num_positions)

    varieties = snp_dictionary["varname"]
    var_allele_list = snp_dictionary["varalleles"]
    for i in range(0,len(var_allele_list)):  # there is a list of lists, we will loop through the outer list, just grabbing the list of alleles for of the 3116 varieties
        alleles_for_one_variety = var_allele_list[i]
        variety = varieties[i]

        if variety != None:
            if "－" in variety:
                variety = variety.replace("－","-")    #UTF encoding problem with this char in some names

            #COLOMBIA XXI::G1::[IRIS 313-15908]
            variety = variety.split("[")[1][:-1]
            var_alleles[variety] = alleles_for_one_variety ##++
        #print(variety,alleles_for_one_variety,end="\n")
    return var_alleles

## Converts json to separated text
##2021 fixing bug where one variety name has a comma in it, so need csv format
def convert_json_text(locus_to_process):
    success= True
    f_error_log = open("error_log.txt","w")
    outfolder = PATH_TO_DATABASE_FOLDER + 'temp_csv_files/'
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    filename = IN_FOLDER + locus_to_process + ".json"
    print("Converting:",filename ," to csv")
    with open(filename, "r") as json_file:

        outfile = outfolder + "/" + locus_to_process + ".csv"
        snp_dictionary = json.load(json_file)

        chromosome = get_chromosome(snp_dictionary)

        if chromosome != None:
            #print("chr: ", chromosome)
            snp_positions = get_position_on_chr_as_list(snp_dictionary)
            #print("positions,",snp_positions)
            ref_bases = snp_dictionary["reference"]
            #print("ref bases:" ,ref_bases)
            #varieties = snp_dictionary["varname"]  # loop-
            variety_alleles = get_variety_alleles(snp_dictionary)
            varieties = list(variety_alleles.keys())

            #print("varieties",varieties)
            snp_array = np.empty(shape=(len(varieties)+2,len(snp_positions)+1), dtype=('U50'))
            snp_array[0,0] = "Variety"

            ref_bases = snp_dictionary["reference"]
            #print("ref bases:" ,ref_bases)

            snp_array[1,0] = "Reference"
            for i in range(0,len(ref_bases)):
                row_pos = i+1
                snp_array[1,row_pos] = ref_bases[i]

            for row_variety in range(0,len(varieties)):

                variety = varieties[row_variety]
                variety_allele_list = variety_alleles[variety]
                snp_array[row_variety + 2, 0] = variety

                for col_snp in range(0,len(snp_positions)):
                    snp_pos = snp_positions[col_snp]

                    if row_variety == 0:
                        snp_array[0, col_snp+1] = snp_pos  #fill header row, first two positions taken
                    variety_allele = variety_allele_list[col_snp]
                    snp_array[row_variety+2,col_snp+1] = variety_allele

            json_file.close()
            np.savetxt(outfile,snp_array,delimiter=",",fmt='%s')
        else:
            print("No data for locus:",locus_to_process)
            success=False
    return success



#convert_json_text("LOC_Os01g48060")