'''
AVA-seq main pipeline. Set main parameters.
'''
##############################################
# INPUT FILES

# organization of files, replicates, libraries, etc 
# IMPORTANT: project name and lbrary name should not contain underscore "_"
organization_file = "organization1.txt"

# path to fastq.gz files location
fastaq_files_location = "../hpPRS/TestFastq/"

# path to edgeR script
edgeR_script = "scriptToRunR.R"

# path to diamond database
diamond_database = "nr.dmnd"

# positive controls file path
# every protein pair should be written in new line
positive_controls_location = "positive_controls_example.txt"

# negative controls file path
# every protein pair should be written in new line
negative_controls_location = "negative_controls_example.txt"

# path to frame shift names file location
# every FS fragment should be written in new line
fs_location = "frame_shift_example.txt"

##############################################
# OUTPUT FILES

# intermediate files in the initial pipeline (like .pep. .joinbp and .count files) will be here		
intermediate_files_location = "../hpPRS/WasteFiles/"	

# .counts files and original .diff files will be here
counts_files_location = "../hpPRS/CountsFilesTest/"

# output directory name
# .diff files after FS removal, reports and tables will be here
results_directory_location = "../hpPRS/Results/"
##############################################
# INTERACTION PARAMETERS

# cutoffs 
logFC_cutoff = 5
FDR_cutoff = 0.05
##############################################
# AVASEQ DESIGN DETAILS

ava_lambda  = "ACGTTTGGC" # original "CACAAGGG"
ava_rnap = "GAGGCGGCC"    # original with added T for fragments:  TCGTTTTGG
distance_to_fragment_start = 36
length_of_read_used = 75
##############################################
# STAGES

# Do you want fastq.gz --> counts files. Set True for 'Yes', False for 'No'
flag1 = True
# Do you want counts --> diff files. Set True for 'Yes', False for 'No'
flag2 = True
# Do you want diff --> controls and frame shift removal. Set True for 'Yes', False for 'No'
flag3 = True
# Do you want final PPI table. Set True for 'Yes', False for 'No'
flag4 = True

# Make sure not to set flag2 if flag1 is not set, 
# and there are not any previous .counts files to be processed.
# Make sure not to set flag3 and/or flag4 if flag1 and flag2 are not set,
# and there are not any previous .diff files to be processed.
# Combination of flag1, flag3 and/or flag4 - do not set this! 
# The computations depending on flag3 and flag4 will ignore new counts files.
##############################################
# DIAMOND 

# To create diamond nr.dmnd reference database for alignment
# Run in the main pipeline folder: "diamond makedb --in seqeunces_protein.faa -d nr"
# sequences_protein.faa is a protein database file in FASTA format
##############################################