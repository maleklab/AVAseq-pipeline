'''
#########################################################
Main processing file. 									#
There are four sub-parts.								#
														#
1. Initial pipeline translated to Python				#
Input: FastaQ files  Output: .counts files. 			#
														#
ORDER of functions         								#
1.1. fastqToFasta  	_p           						#
1.2. fastaToPep              							#			
1.3. pepToJoin  	_p             						#
1.4. pepJoinToDifflambdaFirst   _p						#
1.5. pepJoinDiffToFasta			_p      				#
1.6. forRevToBlastDIAM       							#
1.7. bpToJoin  		_p              					#
1.8. get_counts 	_p              					#
														#					
														#
2. edgeR_script                                     	#
Input: .counts files. Output: .diff files. 				#
														#
3. Additional processing (+/- controls, FS removal)		#
Input:. diff files  Output: Processed .diff files.		#
														#
4. Protein-Protein interactions table               	#
Input: .diff files  Output: one .txt file 				#
														#
														#
Define INPUT and OUTPUT parameters in parameters.py. 	#
Functions with _p next to the name have two versions: 	# 
Less pythonic without _p and more pythonic with _p.		#
#########################################################
'''

import os
import subprocess
import shutil
from utils import *        
from parameters import *   # parameters to be set by the user


def creation_of_output_directories(results_directory_location, intermediate_files_location, counts_files_location):

	# directory for intermediate files of the initial pipeline
	isExist = os.path.exists(intermediate_files_location)
	if not isExist:
		os.mkdir(intermediate_files_location)
		print("Directory '% s' created" % intermediate_files_location)
	else:
		print("Directory '% s' already exists!" % intermediate_files_location)


	# directory for .counts files 
	isExist = os.path.exists(counts_files_location)
	if not isExist:
		os.mkdir(counts_files_location)
		print("Directory '% s' created" % counts_files_location)
	else: 
		print("Directory '% s' already exists!" % counts_files_location)


	# directory where reports, tables and final diff files go
	isExist = os.path.exists(results_directory_location)
	if not isExist:
		os.mkdir(results_directory_location)
		print("Directory '% s' created" % results_directory_location)
	else:
		print("Directory '% s' already exists!" % results_directory_location)




def main(fastaq_files_location):

	organization_df = pd.read_csv(organization_file, sep="\t", header=0)

	# part 1
	if flag1:
		# Initial pipeline translated to Python
		print("FASTQ --> .counts --------------------------------------------------------------------")

		for i in organization_df.index:

			file_R1 = fastaq_files_location + str(organization_df.loc[i, 'R1'])
			file_R2 = fastaq_files_location + str(organization_df.loc[i, 'R2'])    

			print(f"Files in processing:\n{file_R1}\n{file_R2}")

			fasta_R1 = fastqToFasta_p(file_R1)
			fasta_R2 = fastqToFasta_p(file_R2)

			pep_R1 = fastaToPep(fasta_R1, ava_lambda, ava_rnap, distance_to_fragment_start, length_of_read_used)
			pep_R2 = fastaToPep(fasta_R2, ava_lambda, ava_rnap, distance_to_fragment_start, length_of_read_used)

			pep_join = pepToJoin_p(pep_R1, pep_R2)
			pep_diff = pepJoinToDifflambdaFirst_p(pep_join)
			join_diff_forr, join_diff_rev = pepJoinDiffToFasta_p(pep_diff)

			join_diff_forr_bp = forRevToBlastDIAM(join_diff_forr, diamond_database)
			join_diff_rev_bp = forRevToBlastDIAM(join_diff_rev, diamond_database)

			joinbp = bpToJoin_p(join_diff_forr_bp, join_diff_rev_bp)

			get_counts_p(joinbp)

		# moving all files except original .fastq.gz and .count to different folder
		for file in glob.iglob(fastaq_files_location + "*"):
			if file.endswith(".fastq.gz") or file.endswith(".count"):
				continue
			else:
				shutil.move(file, intermediate_files_location + os.path.basename(file))

		# move .count files to a specified folder
		for file in glob.iglob(fastaq_files_location + "*.count"):
			shutil.move(file, counts_files_location + os.path.basename(file))

		merge_counts(organization_df, counts_files_location)

		# after consolidation move .count files to intermediate folder 
		for file in glob.iglob(counts_files_location + "*.count"):
			shutil.move(file, intermediate_files_location + os.path.basename(file))

	print(".counts files created -----------------------------------------------------------------")



	# part 2 
	if flag2:
		print("edgeR statistical analysis ------------------------------------------------------------")
		subprocess.run(["Rscript", edgeR_script, counts_files_location])  
		print("edgeR completed -----------------------------------------------------------------------")
	diff_files_location = counts_files_location  # diff files are at the same location as .counts, can be changed





	# part 3
	if flag3:
		print("+/- controls --------------------------------------------------------------------------")
		controls_report(diff_files_location, positive_controls_location, negative_controls_location, logFC_cutoff, FDR_cutoff, results_directory_location)  # positive and negative controls
		print("Controls report created ---------------------------------------------------------------")
		removalFS(diff_files_location, fs_location, logFC_cutoff, FDR_cutoff, results_directory_location) # fragments which interact with frame shift removal
		print("Removal of interacting frame shift fragments completed---------------------------------")





	# part 4
	if flag4 and flag3:
		print("Creating Protein - Protein table of interactions----------------------------------------")
		create_protprot_table(results_directory_location, results_directory_location, logFC_cutoff, FDR_cutoff)
		print("Protein - Protein interactions table DONE\nPROCESSING FINISHED")
	if flag4 and not flag3:
		print("Creating Protein - Protein table of interactions----------------------------------------")
		create_protprot_table(diff_files_location, results_directory_location, logFC_cutoff, FDR_cutoff)
		print("Protein - Protein interactions table DONE\nPROCESSING FINISHED")




if __name__=="__main__":
	creation_of_output_directories(results_directory_location, intermediate_files_location, counts_files_location)
	main(fastaq_files_location)




