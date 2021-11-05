Repo setup needs following

1. In Pycharm, make sure the folder Databases is excluded for index go to folder, right click - Mark Excluded - we are going to put loads of mini files in here
which will kill PyCharm
2. Need to run file_downloader.py first, to make sure all input database files are either downloaded fresh or unzipped from versions stored in GitHub
3. Need to run the chromosome splitter, so that the in the Databases/ChromosomeSplit folder, you get one file per chromosome, like this Oryza_sativa.IRGSP-1.0.dna.toplevel.fa1.fasta
4. then run gff_splitter. This splits up a big gff files in 10,000s of tiny files, so we can read them fast later on Databases/cached_gff_chunks
5. then run CreateProteinAlignmentFromSNPSeek - need to change the file referenced data_loc, which is the input data file.
This lists all the gene/transcript IDs (one per line)
The code then goes off to SNP_seek, gets SNPs from 3K rice genomes, code recreates protein sequence from GFF snippet and correct chromosome,
via CDS coordinates.

Nucleotide changes are then switched in to create amino acid polymorphisms.

This calls another routine which counts the groups of rice varieties within the 3K set that carry each amino acid change at each position




