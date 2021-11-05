# rice-snp-ptm-analysis
Analysis of SNPs and PTMs for rice

Setup.

Need to download zip file from https://drive.google.com/file/d/1ibemgNFD_2tT2rRUfxcRU51ZGnRN12ij/view?usp=sharing. 
This file should be unzipped one folder back from the main package, i.e. see parameter in code: PATH_TO_DATABASE_FOLDER = "../../RiceDatabases/"

The primary run for the code is done via CreateProteinAlignmentFromSNPSeek.py

This reads run_details.txt, which has two parameters, one per line:
1. Location of input file
2. Location of where to write output

The input file should have one RAP-DB or MSU transcript identifier per line, like this:

'
Os05t0392300-02
Os07t0694000-01
Os01t0549400-01
Os09t0103700-01
Os01t0614500-01
Os01t0588500-01
LOC_Os08g05540.1
Os08t0178300-01
Os01t0971900-01
LOC_Os07g30840.2
Os06t0232000-01
'

The code is lazily multi-threaded, so it will try to start as many threads as there are rows in the input file. This needs improving, so it should only start up a sensible number of threads.

Alternatively, to process lots of data, the analysis could be chunked in sets of 50 input transcripts, and multiple runs of the code be done (in series or parallel).

The code goes off to SNP-Seek and downloads SNP info in json format. It then extracts the coordinates and bases. It recreates the protein sequence using gff snippets (in Databases) for any given gene, then switches in the base changes to create a new protein molecule.

If the json has been downloaded before, the code uses the cached (local) version.


