
import requests
import os
import json

PATH_TO_DATABASE_FOLDER = "../../RiceDatabases/"

def retrieve_json(locusname):

    success = True

    if locusname.find("."):  # Need to remove suffix for SNP-SEEK   #TODO does this work for RAP-DB transcript IDs?
        locusname = locusname[:locusname.find(".")]

    folder_name = PATH_TO_DATABASE_FOLDER + 'cached_json_files'
    outfile = folder_name + "/" + locusname + ".json"

    if not os.path.exists(outfile):

        print("requesting:",locusname, " from SNP-Seek")
        # r = requests.get('http://snp-seek.irri.org/ws/genotype/gettable?snp=true&locus=OS02G0141100')
        r = requests.get('http://snp-seek.irri.org/ws/genotype/gettable?snp=true&locus=' + locusname)


        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        if r:
            data = r.json()
            f_out = open(outfile, "w")
            json.dump(data, f_out)
        else:
            print("Error no response from SNP-Seek for",locusname)
            success=False
    else:
        print(outfile," already exists, no need to download again")

    return [locusname,success]

#retrieve_json("Os01g0625300.1")