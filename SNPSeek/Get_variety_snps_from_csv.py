###Script to read csv (from json) files

###To do - write method to convert for all varieties

#Converts csv of SNPs into a dictionary, keys = varieties; values = string of SNPs
def get_variety_to_snps_dictionary(snp_csv_file):

    varieties_to_snp_string = {}

    try:
        f = open(snp_csv_file, "r")
        snp_positions = []
        found_snps = []

        counter = 0
        found_var_name = False
        for line in f:
            snp_string = ""
            varname = ""
            if found_var_name == False:
                line = line[:-1]
                cells = line.split(",")
                if counter == 0:
                    snp_positions = cells[1:]  # grab all SNP positions
                else:
                    varname = cells[0]
                    found_snps = cells[1:]
                    if len(snp_positions) != len(found_snps):
                        print("Error for variety name ", varname, " bases and SNP position lengths do not match")
                        exit()

                    for i in range(0, len(snp_positions)):
                        pos = snp_positions[i]
                        base = found_snps[i]
                        snp_string += pos + ":" + base + ";"

                    if len(snp_string) > 1:
                        snp_string = snp_string[:-1]  # remove final ;
                    varieties_to_snp_string[varname] = snp_string
            counter += 1
    except IOError:
        print ("Error: File does not exist",snp_csv_file)

    return varieties_to_snp_string


#Converts tsv of SNPs into a dictionary, keys = varieties; values = string of SNPs, new method to handle new format from SNP-Seek
def get_variety_to_snps_dictionary_local(csv_from_json_file):

    varieties_to_snp_string = {}
    f = open(csv_from_json_file, "r")
    #print("opening file",tsv_from_json_file)
    snp_positions = []
    found_snps = []

    counter = 0
    found_var_name = False
    for line in f:
        snp_string = ""
        varname = ""
        if found_var_name == False:
            line = line[:-1]
            cells = line.split(",")
            if counter == 0:
                snp_positions = cells[1:]  # grab all SNP positions
            else:
                varname = cells[0]
                found_snps = cells[1:]
                if len(snp_positions) != len(found_snps):
                    print("Error for variety name ", varname, " bases and SNP position lengths do not match")
                    exit()

                for i in range(0, len(snp_positions)):
                    pos = snp_positions[i]
                    base = found_snps[i]
                    snp_string += pos + ":" + base + ";"

                if len(snp_string) > 1:
                    snp_string = snp_string[:-1]  # remove final ;
                varieties_to_snp_string[varname] = snp_string
        counter += 1


    return varieties_to_snp_string


#Method to extract the SNPs from a given CSV file for a given variety and convert to string format
def read_csv_convert_to_string(csv_from_json_file,variety_name):
    snp_string = ""

    f = open(csv_from_json_file,"r")
    snp_positions = []
    found_snps = []

    counter = 0
    found_var_name = False
    for line in f:

        if found_var_name==False:
            line = line[:-1]
            cells = line.split(",")
            if counter ==0:
                snp_positions = cells[1:] #grab all SNP positions
            else:
                varname = cells[0]
                if variety_name in varname:
                    found_var_name = True
                    found_snps = cells[1:]
        counter +=1

    if len(snp_positions) != len(found_snps):
        print("Error for variety name ",variety_name," bases and SNP position lengths do not match")

    for i in range(0,len(snp_positions)):
        pos = snp_positions[i]
        base = found_snps[i]
        snp_string += pos + ":" + base + ";"

    if len(snp_string) > 1:
        snp_string = snp_string[:-1]    #remove final ;
    return snp_string


#snp_string = read_csv_convert_to_string("temp_csv_files/LOC_Os01g48060.csv","COLOMBIA")
#print(snp_string)

