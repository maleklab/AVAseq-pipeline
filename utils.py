import sys
import os
import subprocess
import glob
import re
import pandas as pd
import numpy as np


# input .fastq.gz file --> output .fasta file
def fastqToFasta(file):
	if re.search(".gz$",file):   # matching .gz at the end of the file
		unziped_file = file[:-3]   # unpacked .gz 
		subprocess.run(f"gunzip -c {file} > {unziped_file}", shell=True, capture_output=True)
		out = unziped_file.replace('.fastq','.fasta')  # file to store .fasta, replace .fastq extension with .fasta
		# There are other options like sed
		s = subprocess.run(f"perl -ne 'y/@/>/;print($_.<>)&&<>&&<>' {unziped_file} > {out}", capture_output=True, text=True, shell=True)
		if s.returncode == 0:
			print(".fastq --> .fasta DONE")
	else:
		print("Files not in the right format!")
	return out 




def fastqToFasta_p(file):
	if re.search(".gz$",file):   # matching .gz at the end of the file
		unziped_file = file[:-3]   # unpacked .gz --> it has a name without .gz
		subprocess.run(f"gunzip -c {file} > {unziped_file}", shell=True, capture_output=True)
		output_file_name = unziped_file.replace('.fastq','.fasta')
		output_file = open(output_file_name, "w")

		for line in open(unziped_file, "r").readlines():
			if line[0]=="@":
				output_file.write(">" + line[1:])
			elif line[0] in ['A', 'C', 'T', 'G']:
				output_file.write(line)
		print(".fastq --> .fasta DONE")
	else:
		print("Files not in the right format!")
		output_file_name = "" # breaks the code
	return output_file_name




##############################################################################################################
# input .fasta file --> output .pep file
# Take fasta files of AVAseq and return proteins
def fastaToPep(file, ava_lambda, ava_rnap, distance_to_fragment_start, length_of_read_used):
	output_file_name = file.replace(".fasta", ".pep")
	output_file = open(output_file_name, "w")

	codon2aa = {'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TTC':'F','TTT':'F','TTA':'L',
				'TTG':'L','TAC':'Y','TAT':'Y','TAA':'*','TAG':'*','TGC':'C','TGT':'C',
				'TGA':'*','TGG':'W','CTA':'L','CTC':'L','CTG':'L','CTT':'L','CCA':'P',
				'CCC':'P','CCG':'P','CCT':'P','CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
				'CGA':'R','CGC':'R','CGG':'R','CGT':'R','ATA':'I','ATC':'I','ATT':'I',
				'ATG':'M','ACA':'T','ACC':'T','ACG':'T','ACT':'T','AAC':'N','AAT':'N',
				'AAA':'K','AAG':'K','AGC':'S','AGT':'S','AGA':'R','AGG':'R','GTA':'V',
				'GTC':'V','GTG':'V','GTT':'V','GCA':'A','GCC':'A','GCG':'A','GCT':'A',
				'GAC':'D','GAT':'D','GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G',
				'GGT':'G'}


	for line in open(file, "r").readlines():
		line = line.replace("\n", "")
		if line[0]==">":
			output_file.write(line.split(" ")[0] + "\t")
		else:
			locLam = line.rfind(ava_lambda) + distance_to_fragment_start   # rfind returns -1 if not found in the sequence
			locRnap = line.rfind(ava_rnap) + distance_to_fragment_start
			if locLam > locRnap:
				start = locLam
				note = "lambda"
			elif locRnap > locLam:
				start = locRnap
				note = "rnap"
			else:
				start = 0
				note = "neither"
			coding = line[start:(start+length_of_read_used)]
			DNA = coding
			protein = ''
			for i in range(0, len(DNA)-2, 3):
				codon = DNA[i:(i+3)]
				if codon.upper() in codon2aa.keys():
					protein += codon2aa[codon.upper()]
				else:
					protein += "X"
			output_file.write(note + "\t" + protein + "\n")
	print(".fasta --> .pep DONE")
	return output_file_name




##############################################################################################################
# take file of R1 and then join to R2 peps
def pepToJoin(file_r1_pep, file_r2_pep):
	out = file_r1_pep.replace('.pep','Join.pep').replace('_R1_','_R1R2_')
	s = subprocess.run(f"join -j1 {file_r1_pep} {file_r2_pep} > {out}", capture_output=True, shell=True, text=True)
	if s.returncode == 0:
				print("Join of .pep R1 and R2 DONE")
	return out



def pepToJoin_p(file_r1_pep, file_r2_pep):
	output_file_name = file_r1_pep.replace('.pep','Join.pep').replace('_R1_','_R1R2_')
	output_file = open(output_file_name, "w")
	df_r1 = pd.read_csv(file_r1_pep, sep="\t", names = ['Index', 'lr', 'AA'])
	df_r2 = pd.read_csv(file_r2_pep, sep="\t", names = ['Index', 'lr', 'AA'])
	df_out = pd.merge(df_r1, df_r2, on='Index')
	df_out.to_csv(output_file, sep = " ", index=False, header=False) 
	print("Join of .pep R1 and R2 DONE")
	return output_file_name



##############################################################################################################
# remove pairs that dont have lambdacI and RNAP both
def pepJoinToDifflambdaFirst(file):
	out = file.replace('.pep','Diff.pep')
	cmd = f'''cat {file} | grep -v 'neither' | awk '$2!=$4' | awk '{{if ($2=="lambda"){{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}} else if ($2 != "lambda"){{print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3}}}}' > {out}'''
	s = subprocess.run(cmd, capture_output=True, shell=True, text=True)
	if s.returncode==0:
		print(".pep --> Diff.pep DONE")
	return out



def pepJoinToDifflambdaFirst_p(file):
	output_file_name = file.replace('.pep','Diff.pep')
	output_file = open(output_file_name, "w")
	df = pd.read_csv(file, sep = " ", names = ['Index', 'lr1', 'AA1', 'lr2', 'AA2'])
	df = df[(df['lr1']!='neither') & (df['lr2']!='neither') & (df['lr1']!=df['lr2'])]
	df[['AA1','AA2']] = df[['AA2','AA1']].where(df['lr1'] == 'rnap', df[['AA1','AA2']].values)
	df[['lr1','lr2']] = df[['lr2','lr1']].where(df['lr1'] == 'rnap', df[['lr1','lr2']].values)
	df.to_csv(output_file, sep="\t", index=False, header=False)
	print(".pep --> Diff.pep DONE")
	return output_file_name



##############################################################################################################
# diff_pep file input --> .for and .rev output
def pepJoinDiffToFasta(file):
	forr = file.replace('.pep','.for')
	rev  = file.replace('.pep','.rev')
	cmd  = f"""perl -e 'open(FILE, "<{file}"); open(FORR, ">{forr}"); open(REV, ">{rev}"); while($li=<FILE>){{chomp; my @line=split(/\\s+/,$li); print FORR "$line[0]\\n$line[2]\\n"; print REV "$line[0]\\n$line[4]\\n";}}' """
	s = subprocess.run(cmd, capture_output=True, shell=True, text=True)
	if s.returncode==0:
		print("Diff.pep --> .for and .rev DONE")	
	return forr, rev



def pepJoinDiffToFasta_p(file):
	forr_name = file.replace('.pep','.for')
	forr = open(forr_name, 'w')
	rev_name = file.replace('.pep','.rev')
	rev = open(rev_name, 'w')
	df = pd.read_csv(file, sep="\t", names=['Index', 'lambda', 'AA1', 'rnap', 'AA2'])
	for i in df.index:
		forr.write(str(df.loc[i, 'Index']) + "\n" + str(df.loc[i, 'AA1']) + "\n")
		rev.write(str(df.loc[i, 'Index']) + "\n" + str(df.loc[i, 'AA2']) + "\n")
	print("Diff.pep --> .for and .rev DONE")
	return forr_name, rev_name



##############################################################################################################
# .for or .rev files inputs
# nr.dmnd is a database --> created with diamond 
def forRevToBlastDIAM(file, diamond_database):
	out = f"{file}" + ".bp"
	cmd = f"diamond blastp --db {diamond_database} --outfmt 6 -q {file} --matrix PAM30 --more-sensitive --masking 0 --max-target-seqs 1 | awk '$3>80&&$7==1' > {out}"
	s = subprocess.run(cmd, capture_output=True, shell=True, text=True)
	if s.returncode==0:
		print("Diamond alignment DONE")
	return out



##############################################################################################################

def bpToJoin(file_forbp, file_revbp):
	out = file_forbp.replace('.for.bp','.joinbp')
	forr_sorted = file_forbp.replace('.for.bp','sorted.for.bp')
	rev_sorted = file_revbp.replace('.rev.bp','sorted.rev.bp')
	s1 = subprocess.run(f"sort -k1,1 {file_forbp} > {forr_sorted}", capture_output=True, shell=True, text=True)
	s2 = subprocess.run(f"sort -k1,1 {file_revbp} > {rev_sorted}", capture_output=True, shell=True, text=True)
	cmd = f"join -1 1 -2 1 -t $'\\t' -o 1.1,1.2,1.9,2.2,2.9 {forr_sorted} {rev_sorted} > {out}" 
	s = subprocess.run(cmd, capture_output=True, shell=True, text=True)
	if s.returncode==0:
		print(".for and .rev --> .joinbp DONE")
	return out



def bpToJoin_p(file_forbp, file_revbp):
	output_file_name = file_forbp.replace('.for.bp','.joinbp')
	df_forr = pd.read_csv(file_forbp, sep="\t", header=None)
	df_forr.sort_values(by=df_forr.columns[0], inplace=True)
	df_rev = pd.read_csv(file_revbp, sep="\t", header=None)
	df_rev.sort_values(by=df_rev.columns[0], inplace=True)
	df = pd.merge(df_forr, df_rev, on=df_forr.columns[0]).iloc[:, [0,1,8,12,19]]
	df.to_csv(output_file_name, sep="\t", index=False, header=False)
	print(".for.bp and .rev.bp --> .joinbp DONE")
	return output_file_name



##############################################################################################################

def get_counts(file):
	dict_hash={}
	for line in open(file, "r").readlines():
		line = line.replace("\n", "")
		t = line.split("\t")[1:5]
		t = [t[0] + ':' + t[1], t[2] + ':' + t[3]] # [protA:x, protB:y]
		pair = t[0] + ':' + t[1]
		if pair in dict_hash.keys():
			dict_hash[pair] += 1      # dict{protA:x:protB:y}
		else:
			dict_hash[pair] = 1 
	output_file = open(f"{file}.count",'w')
	# sort the fields so that proteinA:fragmentx with proteinB:fragmenty is always reported as proteinA:fragmentx:proteinB:fragmenty
	dict_hash = dict(sorted(dict_hash.items()))
	for pair in dict_hash.keys():
		output_file.write(str(pair) + '\t' + str(dict_hash[pair]) + '\n')
	print(f"Counts DONE for file: {file}\n")



def get_counts_p(file):
	output_file = open(f"{file}.count",'w')
	df = pd.read_csv(file, sep="\t", names=['Index', 'fragment1', 'start1', 'fragment2', 'start2'])
	df['Pair'] = df['fragment1'] + ":" + df['start1'].astype(str) + ":" + df['fragment2'] + ":" + df['start2'].astype(str)
	df['Pair'].value_counts().sort_index().to_csv(output_file, sep="\t", header=False)
	print("Counts DONE for file: ", file)




##############################################################################################################
# Merge all separate replicate files for one library into one table
def merge_counts(counts_files_location):

	list_of_libraries = []
	for filepath in glob.iglob(counts_files_location + "*.count"):
		lib = filepath.split("/")[-1].split("_")[1]
		if lib not in list_of_libraries:
			list_of_libraries.append(filepath.split("/")[-1].split("_")[1])


	for lib in list_of_libraries:  # list of libraries
		for replicate in range(1,10):
			for filepath in glob.iglob(counts_files_location + f"*_{lib}_{replicate}*.count"):   # must be consistent with name formatting
				df = pd.read_csv(filepath, sep="\t", names=['FragPair', os.path.basename(filepath)])

				if replicate==1:
					df_main = df
				else:
					df_main = df_main.merge(df, how='outer', on='FragPair').fillna(0)
		df_main.iloc[:, 1:] = df_main.iloc[:, 1:].astype(int)
		df_main.to_csv(counts_files_location + "PRS_" + f"{lib}" + ".counts", sep="\t", index=False)




##############################################################################################################

# There is a list of expected positive and negative interactions 
# Extract .diff files information about interactions and create a report
def controls_report(diff_files_location, positive_controls_location, negative_controls_location, logFC_cutoff, FDR_cutoff, output_directory):

	# list of positive interactions
	positive_controls = [x.replace("\n","") for x in open(positive_controls_location, "r").readlines()]

	# list of negative interations
	negative_controls= [x.replace("\n","") for x in open(negative_controls_location, "r").readlines()]

	# reports file
	output_file = open(output_directory + "controls_report.txt", "w")
	list_pos_contr_data = []
	list_neg_contr_data = []


	for diff_file in glob.iglob(diff_files_location + "*diff.txt"):
		df = pd.read_csv(diff_file, sep="\t")
		# interacting dataframe
		df_int = df[(df['logFC']>logFC_cutoff) & (df['FDR']<FDR_cutoff)]
		broken_df = df_int.loc[:, 'genes'].str.split(pat=":", expand = True)
		df_int.loc[:, 'gene1'] = broken_df.iloc[:, 0]
		df_int.loc[:, 'frag1'] = broken_df.iloc[:, 0] + ':' + broken_df.iloc[:, 1]
		df_int.loc[:, 'gene2'] = broken_df.iloc[:, 2]
		df_int.loc[:, 'frag2'] = broken_df.iloc[:, 2] + ':' + broken_df.iloc[:, 3]

		df_int.loc[:, 'ProtPair'] = df_int.loc[:, 'gene1'] + ":" + df_int.loc[:, 'gene2']
		df_int.loc[:, 'ProtPairR'] = df_int.loc[:, 'gene2'] + ":" + df_int.loc[:, 'gene1']

		
		for ind in df_int.index:
			protpair = df_int.loc[ind, 'ProtPair']
			protpair_r = df_int.loc[ind, 'ProtPairR']

			if protpair in positive_controls:
				list_pos_contr_data.append([protpair, df_int.loc[ind, 'genes'], df_int.loc[ind, 'logFC'], df_int.loc[ind, 'FDR'], os.path.basename(diff_file)])
			elif protpair_r in positive_controls:
				list_pos_contr_data.append([protpair_r, df_int.loc[ind, 'genes'], df_int.loc[ind, 'logFC'], df_int.loc[ind, 'FDR'], os.path.basename(diff_file)])

			if protpair in negative_controls:
				list_neg_contr_data.append([protpair, df_int.loc[ind, 'genes'], df_int.loc[ind, 'logFC'], df_int.loc[ind, 'FDR'], os.path.basename(diff_file)])
			elif protpair_r in negative_controls:
				list_neg_contr_data.append([protpair_r, df_int.loc[ind, 'genes'], df_int.loc[ind, 'logFC'], df_int.loc[ind, 'FDR'], os.path.basename(diff_file)])


	df_positive_controls = pd.DataFrame(list_pos_contr_data, columns=["ProtPair", "FragPair", "logFC", "FDR", "File"])
	df_positive_controls.sort_values(by=['ProtPair'], inplace=True)

	df_negative_controls = pd.DataFrame(list_neg_contr_data, columns=["ProtPair", "FragPair", "logFC", "FDR", "File"])
	df_negative_controls.sort_values(by=['ProtPair'], inplace=True)


	# write to a file
	output_file.write("**********POSITIVE CONTROLS**********\n\n")
	output_file.write("EXPECTED INTERACTIONS obtained for the following protein pairs.\n\n")

	output_file.write("ProtPair\tFragPair\tlogFC\tFDR\tFile\n")
	for i in df_positive_controls.index:
		output_file.write(df_positive_controls.loc[i, 'ProtPair'] + "\t" + df_positive_controls.loc[i, 'FragPair'] + "\t" + df_positive_controls.loc[i, 'logFC'].astype("str") + "\t" + df_positive_controls.loc[i, 'FDR'].astype("str") + "\t" + df_positive_controls.loc[i, 'File'] + "\n")
		
	confirmed_positive = df_positive_controls['ProtPair'].unique().tolist()
	non_confirmed_positive = list(set(positive_controls) - set(confirmed_positive))
	if len(non_confirmed_positive)>0:
		output_file.write("\nMISSED EXPECTED INTERACTIONS:\n")
		for i in non_confirmed_positive:
			output_file.write(i + "\n")	


	output_file.write("\n\n\n**********NEGATIVE CONTROLS**********\n\n")
	output_file.write("UNEXPECTED INTERACTIONS obtained for the following protein pairs.\n\n")
	output_file.write("ProtPair\tFragPair\tlogFC\tFDR\tFile\n")

	for i in df_negative_controls.index:
		output_file.write(df_negative_controls.loc[i, 'ProtPair'] + "\t" + df_negative_controls.loc[i, 'FragPair'] + "\t" + df_negative_controls.loc[i, 'logFC'].astype("str") + "\t" + df_negative_controls.loc[i, 'FDR'].astype("str") + "\t" + df_negative_controls.loc[i, 'File'] + "\n")
	

	wrong_positive = df_negative_controls['ProtPair'].unique().tolist()
	confirmed_negative = list(set(positive_controls) - set(confirmed_positive))
	if len(confirmed_negative)>0:
		output_file.write("\nCORRECT NON-INTERACTIONS:\n")
		for i in confirmed_negative:
			output_file.write(i + "\n")	




##############################################################################################################
# There is a list of frame shift (FS) names. Remove all fragments from .counts files which 
# interact with FS. They are auto-activators.
def removalFS(diff_files_location, fs_location, logFC_cutoff, FDR_cutoff, output_directory):
	frame_shift = open(fs_location, "r")
	frame_shift_names = [x.replace("\n","") for x in frame_shift.readlines()]

	list_of_problematic_fragments = []

	for diff_file in glob.iglob(diff_files_location + "diff.txt"):
		df = pd.read_csv(diff_file, sep="\t")
		broken_df = df['genes'].str.split(pat=":", expand = True)
		df['gene1'] = broken_df[0]
		df['frag1'] = broken_df[0] + ':' + broken_df[1]
		df['gene2'] = broken_df[2]
		df['frag2'] = broken_df[2] + ':' + broken_df[3]
		# filter only interactions 
		df = df[(df['logFC']>logFC_cutoff) & (df['FDR']< FDR_cutoff)]
		# filter only interactions with either fragment in FS list
		df = df[df['frag1'].isin(frame_shift_names) | df['frag2'].isin(frame_shift_names)]

		for ind in df.index:
			fragment1 = df.loc[i, 'frag1']
			fragment2 = df.loc[i, 'frag2']
			gene1 = df.loc[i, 'gene1']
			gene2 = df.loc[i, 'frag2']

			if fragment1 not in frame_shift_names:
				list_of_problematic_fragments.append(fragment1)
			elif fragment2 not in frame_shift_names:
				list_of_problematic_fragments.append(fragment2)


	for diff_file in glob.iglob(diff_files_location + "*diff.txt"):
		df = pd.read_csv(diff_file, sep="\t")
		broken_df = df['genes'].str.split(pat=":", expand = True)
		df['gene1'] = broken_df[0]
		df['frag1'] = broken_df[0] + ':' + broken_df[1]
		df['gene2'] = broken_df[2]
		df['frag2'] = broken_df[2] + ':' + broken_df[3]
		df = df[~(df['frag1'].isin(list_of_problematic_fragments) | df['frag2'].isin(list_of_problematic_fragments))]
		df.drop(['gene1', 'gene2', 'frag1', 'frag2'], inplace=True, axis=1)
		df.to_csv(output_directory + os.path.basename(diff_file)[:-4] + "withFSremoval.txt", sep = "\t", index=False)




##############################################################################################################
# Create a table of interactions between proteins
def create_protprot_table(directory, logFC_cutoff, FDR_cutoff):  # directory is the location of input files, and for the output file

	# initialize dataframe
	# Orient1 Orient2 2mM 5mM FCmax FDRmin #Libs uniqFragPairs 
	df = pd.DataFrame(columns = ['gene1', 'gene2', 'Orient1', 'Orient2', '2mM', '5mM', 'FCmax', 
		'FDRmin', 'Libs', '#Libs', 'uniqFragPairs', '#uniqFragPairs'])

	for filepath in glob.iglob(output_directory + "*withFSremoval.txt"):  # DIFF AFTER FS REMOVAL
		print("Filepath: ", filepath)

		condition = os.path.basename(filepath).split(".")[-2][0]  

		libs  = os.path.basename(filepath).split("_")[1]  # DOES THIS ALWAYS HOLD?
		# print("Libs: ", libs)

		df_temp = pd.read_csv(filepath, sep = "\t", header = 0)
		df_temp.columns = 'genes logFC logCPM LR PValue FDR'.split()
		broken_df = df_temp['genes'].str.split(pat=":", expand = True)
		df_temp['gene1'] = broken_df[0]
		df_temp['frag1'] = broken_df[0] + ':' + broken_df[1]
		df_temp['gene2'] = broken_df[2]
		df_temp['frag2'] = broken_df[2] + ':' + broken_df[3]
		df_temp = df_temp[(df_temp['logFC']>logFC_cutoff) & (df_temp['FDR']< FDR_cutoff)]

		
		df_temp.reset_index(inplace=True)
		df_temp.drop(['index'], inplace=True, axis=1)

		for index in df_temp.index:
			k=0
			orient1 = 0
			orient2 = 0

			protein1 = df_temp.loc[index, 'gene1']
			protein2 = df_temp.loc[index, 'gene2']
			fragments = (df_temp.loc[index, 'frag1'], df_temp.loc[index, 'frag2'])
			fragments_orient2 = (df_temp.loc[index, 'frag2'], df_temp.loc[index, 'frag1'])

			df_subset_1 = df[(df['gene1']==protein1) & (df['gene2']==protein2)]
			df_subset_2 = df[(df['gene2']==protein1) & (df['gene1']==protein2)]

			len_df_subset_1 = len(df_subset_1.index)
			len_df_subset_2 = len(df_subset_2.index) 
			
			if len_df_subset_1 > 0:
				df_index = df[(df['gene1']==protein1) & (df['gene2']==protein2)].index.values[0]
				k = 1
				orient1 = 1
			elif len_df_subset_2 > 0:
				df_index = df[(df['gene2']==protein1) & (df['gene1']==protein2)].index.values[0]
				k = 1
				orient2 = 1
			
			if k==1:

				if orient1 == 1:
					df.loc[df_index, 'Orient1'] += 1
					if fragments not in df.loc[df_index, 'uniqFragPairs']:
						df.loc[df_index, 'uniqFragPairs'].append(fragments)
						df.loc[df_index, '#uniqFragPairs'] += 1 

				elif orient2 == 1:
					df.loc[df_index, 'Orient2'] += 1
					if fragments_orient2 not in df.loc[df_index, 'uniqFragPairs']:
						df.loc[df_index, 'uniqFragPairs'].append(fragments_orient2)
						df.loc[df_index, '#uniqFragPairs'] += 1 

				if condition=='2':
					df.loc[df_index, '2mM'] += 1
				elif condition=='5':
					df.loc[df_index, '5mM'] += 1

				if libs not in df.loc[df_index, 'Libs']:
					df.loc[df_index, 'Libs'].append(libs)
					df.loc[df_index, '#Libs'] += 1 

				df.loc[df_index, 'FCmax'] = max(df.loc[df_index, 'FCmax'], df_temp.loc[index, 'logFC'])
				df.loc[df_index, 'FDRmin'] = min(df.loc[df_index, 'FDRmin'], df_temp.loc[index, 'FDR'])		

			else:
				df.loc[len(df.index)] = [protein1, protein2, 1, 0, 0, 0, 0, 0.1, [], 0, [], 0] 
				marker = len(df.index)
				if condition == '2':
					df.loc[marker - 1, '2mM'] += 1
				elif condition == '5':
					df.loc[marker - 1, '5mM'] += 1

				df.loc[marker - 1, 'FCmax'] = df_temp.loc[index, 'logFC']
				df.loc[marker - 1, 'FDRmin'] = df_temp.loc[index, 'FDR']
				df.loc[marker - 1, 'uniqFragPairs'].append(fragments)
				df.loc[marker - 1, '#uniqFragPairs'] += 1
				df.loc[marker - 1, 'Libs'].append(libs)
				df.loc[marker - 1, '#Libs'] += 1



	# df = df[(df['Orient1']>0) | (df['Orient2']>0)]
	df.drop(['uniqFragPairs'], inplace=True, axis=1)
	df.drop(['Libs'], inplace=True, axis=1)
	df.columns = 'Protein1 Protein2 Orient1 Orient2 2mM 5mM logFCmax FDRmin #Libs #uniqFragPairs'.split()
	df.to_csv(directory + "protProt_table.txt", index = False, sep = "\t")









