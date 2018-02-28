#!/usr/bin/env python

### This script creates a new genome for models downloaded from Modbase

""" 
************************************************
***              READ MEEEEEEEE              ***
************************************************

Things to change if you want to re-use script...
-there is a variable called "genome_name"...change the value to the desired name you want
-the download all the models you want to use from ModBase to a directory...
the functions grab_pdb_info() and extract_modbase_fasta() capture the model info required from these files, so this should be the working directory for those
- all os.chdir...you need to change the directories in which the folders are created or where you navigate to for each function


FILES TO DOWNLOAD BEFORE PROCEEDING...
-all of your models from ModBase
-IFF YOU DECIDE TO USE UNIPROT MAPPINGS TABLE TO MAP, OTHERWISE DEFAULT IS TO JUST MAP FROM ENSEMBL's API...a reference table from Uniprot which maps your ModBase models to uniprots (rather, gene-id+isoform to UniProt)
	Steps:
		-go to http://www.uniprot.org/database/
		-search for "modbase" in the search field for "cross-referenced databases"
		-click on ModBase
		-on left, select under "Map To": UniProtKB
		-on left, filter by organism...so type name of organism in the search bar and select for the correct one...like "aedes aegypti", be sure to select "go"
		-then click "Download" at the top...paramters should be format = tab-separated and uncompressed
		-copy the downloaded file to folder of choice...make sure the script can find it...

to fix
-handle session timeouts when retreiving uniprots

"""

import os, glob, lzma, re, subprocess, pandas, requests, numpy, errno, shutil
import xml.etree.ElementTree as ET
from os.path import isfile 
from bs4 import BeautifulSoup

# change path below as appropriate in cluster
os.chdir("./")

genome_name = "kn_aedes_aegypti_modbase"


# create the directories in order to set up genome on cluster
def create_dirs():
	# print(os.getcwd())
	
	# if re-running on same genome, you need to delete the old genome
	# the following 2 if-statements deletes the scratch and data
	# directories of the genome

	if os.path.isdir("/ifs/data/c2b2/bh_lab/shares/databases/hfpd/genomes/{}".format(genome_name)):
		shutil.rmtree("/ifs/data/c2b2/bh_lab/shares/databases/hfpd/genomes/{}".format(genome_name), 
			ignore_errors=False, 
			onerror=None)
		shutil.rmtree("/ifs/scratch/c2b2/bh_lab/shares/hfpd/PrePPI/{}".format(genome_name), 
			ignore_errors=False, 
			onerror=None)
	else:
		os.chdir("../../genomes/")	
		os.makedirs("./{}".format(genome_name))	
		os.chdir("./{}".format(genome_name)) 
		dir_list = ['fasta', 'Seqs']
		for dir in dir_list:
			if not os.path.exists(dir):
				os.makedirs(dir, exist_ok = True)


# grab the seqid, mpqs score, structure/model, and geneid from the pdb files
# make a dataframe with these ... later will add fasta sequences via call to function on cluster
def grab_pdb_info():
	# print(os.getcwd())
	os.chdir("/ifs/home/c2b2/ss_lab/kn2404/proj/aegypti_modbase/aedes_aegypti_2016_214300/model")
	
	regex_mpqs = re.compile(r"^REMARK 220 MODPIPE QUALITY SCORE: +([+-]?\d*\.\d*)")
	regex_seqid = re.compile(r"^REMARK 220 MODPIPE SEQUENCE ID: +(\w*)")
	regex_atom = re.compile(r"^(ATOM.*|TER.*)")

	seqid_list = list()
	mpqs_list = list()
	struct_list = list()
	geneid_list = list()
	isoform_list = list()

	xzfile_names = [file for file in glob.glob('*.pdb.xz') if isfile(file)]
	for filename in xzfile_names: 
		with lzma.open(filename, 'r') as f:	
			mpqs, seqid = None, None
			struct = list()	
			file_content = f.readlines()
			for line in file_content:
				line = line.decode('utf-8').rstrip()
				m = regex_atom.match(line)
				if m:
					atom, = m.groups()
					struct.append(atom)
					continue
				m = regex_mpqs.match(line)
				if m:
					mpqs, = m.groups()
					continue
				m = regex_seqid.match(line)	
				if m:
					seqid, = m.groups()
					continue
		geneid = filename[:-7]
		isoform = filename[:-9]	
		seqid_list.append(seqid)
		mpqs_list.append(mpqs)
		struct_list.append(struct)
		geneid_list.append(geneid)
		isoform_list.append(isoform)

	df = pandas.DataFrame({	'geneid': geneid_list,
							'seqid': seqid_list,
							'mpqs': mpqs_list,
							'struct': struct_list,
							'isoform': isoform_list,
							})
	# print(df.head(5))
	# print(df.dtypes)
	# print(df.get_value(1, 'geneid'))
	return (df)

# grab the geneids and get fasta sequences for them	
def extract_modbase_fasta():	
	gene_id, fasta_seq = None, None
	regex_geneid = re.compile(r"^>(.*[^\s])")
	regex_fastaseq = re.compile(r"^([^>].*)")  
	geneid_fasta_dict = dict()

	xzfile_names = [file for file in glob.glob('*.pdb.xz') if isfile(file)]
	pdb_files = [file for file in glob.glob('*.pdb') if isfile(file)] 

	# extract all of the pdb files...make them readable
	for file in xzfile_names: 
		pdb_name = file[:-3] # this gets the geneid
		output_pdb = open(pdb_name, mode = 'w')
		with lzma.open(file, mode = 'r') as input_pdb_xz:
			file_content = input_pdb_xz.readlines()
			for line in file_content:
				line = line.decode('utf-8')
				output_pdb.write(line)
			output_pdb.close() 

	# call the function to get fasta sequences for each geneid
	for file in pdb_files:
		fasta = list()
		sshproc = subprocess.run(["/ifs/data/c2b2/bh_lab/shares/bin/seq {}".format(file)], shell=True, cwd="./", universal_newlines=True, stdout=subprocess.PIPE).stdout
		for line in sshproc.splitlines():
			m = regex_fastaseq.match(line)
			if m:
				fasta_line, = m.groups()
				fasta.append(fasta_line)
				continue
			m = regex_geneid.match(line)
			if m:
				gene_id, = m.groups()
				continue
		fasta = ''.join(fasta)
		geneid_fasta_dict[gene_id] = fasta
	return (geneid_fasta_dict)

# append the fasta corresponding to each geneid to the dataframe 
# filter the dataframe with all the info so that we get reliable models (filter based on mpqs score)
# create the hfpd id list and append to dataframe for each model
def post_processing(df, geneid_fasta_dict):
	df['fasta'] = df['geneid'].map(geneid_fasta_dict) # appending the fastas
	df['mpqs'] = df['mpqs'].astype(float)
	reliable_df = df.loc[df['mpqs'] > 1.1].copy() # filter for reliable models...defined as having mpqs > 1.1
	# print(reliable_df)
	return (reliable_df)

def add_uniprots(reliable_df):
	isoform_query_list = list()
	isoform_uniprot_dict = dict()
	for row in reliable_df.itertuples():
		isoform_query_list.append(row.geneid[:-2])
		isoform_query_list = list(set(isoform_query_list))
	# print (isoform_query_list)
	with requests.Session() as s:
		for isoform in isoform_query_list:
			url = 'https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=&protein={}'.format(isoform)
			response = s.get(url, stream=True)
			if response.status_code == 200:
				soup = BeautifulSoup(response.content, "lxml")
				uniprot_tag = soup.find_all('accession')
				if len(uniprot_tag) == 0:
					continue
				else:
					uniprot = uniprot_tag[0]
					isoform_uniprot_dict[isoform] = uniprot.text
		reliable_df['uniprot'] = reliable_df['isoform'].map(isoform_uniprot_dict) # appending the uniprots

""" 
	UNCOMMENT AND FIX BELOW IF ALSO WANT TO TRY TO GET UNIPROTS FROM PROTEOME FILE DOWNLOADED FROM UNIPROT KB...THIS CODE IS INCOMPLETE BUT IT'S...USEABLE...

	# uniprot_ref_table = pandas.read_csv("/ifs/home/c2b2/ss_lab/kn2404/proj/aegypti_modbase/aedes_aegypti_2016_214300/uniprot-database%3A%28type%3Amodbase%29.tab",
	# 	sep='\t',
	# 	header=(0))
	# uniprot_ref_table.drop(['Entry name', 'Status', 'Organism', 'Length'], 
	# 	axis = 1, 
	# 	inplace = True)
	# # print(uniprot_ref_table.get_value(3, 'Protein names'))
	# print(uniprot_ref_table)
"""

# add hfpds and delete rows with no uniprots
def post_post_processing(reliable_df):
	reliable_df = reliable_df.dropna()
	hfpd_list = list("%06d" % (i+1) for i in range(reliable_df.shape[0]))
	reliable_df = reliable_df.assign(hfpd = pandas.Series(hfpd_list).values)	
	return (reliable_df)

# write the files to desired folders to ready for PrePPI structural neighbor-search...all of the requisite info should be in the data frame at this point...the function below just writes the data frame rows/columns to the appropriate places now...
def write_files(reliable_df):
	fasta_dir = "../../../../genomes/{}/fasta/".format(genome_name)
	seq_dir = "../../../../genomes/{}/Seqs/".format(genome_name)

	# output the column of hfpds to make id_list file
	header_idlist = ['hfpd']
	reliable_df.to_csv('{}/id_list'.format(fasta_dir), columns = header_idlist, sep = " ", index = False, header = False)

	# output the columns of hfpds and corresponding uniprots to make map_list file
	header_maplist = ['hfpd', 'uniprot']
	reliable_df.to_csv('{}/map_list'.format(fasta_dir), columns = header_maplist, sep = "|", index = False, header = False)

	# edit the map_list file to have the weird delimiter they want "tab+>"
	with open('{}/map_list'.format(fasta_dir), 'r+') as maplist:
		content = maplist.read()
		content = re.sub(r'\|', r'\t>', content, flags = re.M)
		maplist.seek(0)
		maplist.truncate()
		maplist.write(content)
		maplist.close()
		

	# output the fasta sequences for every structure into the fasta directory
	for row in reliable_df.itertuples():
		fasta_filename = "{}".format(row.hfpd)
		with open(os.path.join(fasta_dir, fasta_filename), 'w') as fasta_outfile:
			fasta_outfile.writelines(">hfpd_{}|uniprot_{}|seqid_{}|mpqs_{}|geneid_{}|isoform_{}\n{}\n".format(row.hfpd, row.uniprot, row.seqid, row.mpqs, row.geneid, row.isoform, row.fasta))

	# make a separate directory for each structure/pdb whatever
	# in each directory created, output the structure with the name of the file being the corresponding hfpd
	# print(os.getcwd()) ## /ifs/home/c2b2/ss_lab/kn2404/proj/aegypti_modbase/aedes_aegypti_2016_214300/model
	
	for row in reliable_df.itertuples():
		# make a directory called by it's hfpd for each structure in Seqs
		if not os.path.exists("{}/{}".format(seq_dir, row.hfpd)):
			os.makedirs("{}/{}".format(seq_dir, row.hfpd))

		# make a Model directory under each hfpd directory
		if not os.path.exists("{}/{}/Models/".format(seq_dir, row.hfpd)):
			os.makedirs("{}/{}/Models/".format(seq_dir, row.hfpd))

		# output pdb structure to each hfpd/Model/ folder
		with open("{}/{}/Models/{}_m.pdb".format(seq_dir, row.hfpd, row.hfpd), 'w') as pdb_outfile:
			# pdb_outfile.writelines("{}\n".format(*row.struct, sep = '\n'))
			pdb_outfile.write('\n'.join(row.struct) + '\n')

def create_softlinks(reliable_df):

	# these are the paths for the directories to be made and linked
	source_dir = "/ifs/scratch/c2b2/bh_lab/shares/hfpd/PrePPI/{}/Pipeline".format(
		genome_name
		)
	target_dir = "/ifs/data/c2b2/bh_lab/shares/databases/hfpd/genomes/{}/Pipeline".format(
		genome_name
		)
	# if the source directory doesn't exist, make it
	if not os.path.exists(source_dir):
		os.makedirs(source_dir)
	
	# try to make the symbolic link btw the 2 
	# if the link already exists, remove the existing one and remake
	try:
		os.symlink(source_dir, target_dir)
	except OSError as e:
		if e.errno == errno.EEXIST:
			os.unlink(target_dir)
			os.symlink(source_dir, target_dir)
		else:
			raise e

	# repeat the symbolic link creation process as above 
	# this has to be done for all the hfpds, so below we're iterating 
	# through the dataframe for each hfpd

	for row in reliable_df.itertuples():
		pipeline_hfpd = "Pipeline_{}".format(row.hfpd)
		source_hfpd_dir = "/ifs/scratch/c2b2/bh_lab/shares/hfpd/PrePPI/{}/{}".format(
		genome_name, pipeline_hfpd
		)
		target_hfpd_dir = "/ifs/data/c2b2/bh_lab/shares/databases/hfpd/genomes/{}/Seqs/{}/Pipeline".format(
		genome_name, row.hfpd
		)

		nbr_hfpd_dir = "/ifs/data/c2b2/bh_lab/shares/databases/hfpd/genomes/{}/Seqs/{}/Nbr".format(
		genome_name, row.hfpd
		)
		
		if not os.path.exists(source_hfpd_dir):
			os.makedirs(source_hfpd_dir)

		if not os.path.exists(nbr_hfpd_dir):
			os.makedirs(nbr_hfpd_dir)

		try:
			os.symlink(source_hfpd_dir, target_hfpd_dir)
		except OSError as e:
			if e.errno == errno.EEXIST:
				os.unlink(target_hfpd_dir)
				os.symlink(source_hfpd_dir, target_hfpd_dir)
			else:
				raise e

def main():
	create_dirs()
	df = grab_pdb_info()
	geneid_fasta_dict = extract_modbase_fasta()
	reliable_df = post_processing(df, geneid_fasta_dict)
	add_uniprots(reliable_df)
	reliable_df = post_post_processing(reliable_df)
	write_files(reliable_df)
	create_softlinks(reliable_df)

if __name__ == "__main__":
    main()
	












