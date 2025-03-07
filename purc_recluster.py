#!/usr/bin/env python

logo = """
-------------------------------------------------------------
|                       purc_recluster                      |
|        Pipeline for Untangling Reticulate Complexes       |
|                        version 2.0                        |
|            https://bitbucket.org/crothfels/ppp            |
|															|
|      Fay-Wei Li & Carl J Rothfels & Peter W Schafran      |
|           see purc.py script for more information         |
-------------------------------------------------------------
"""

usage = """

Use this script to recluster the alleles/homeologs from a previous PURC run.

Example: ./purc_recluster.py -f purc_3_annotated.fa -o recluster -c 0.997 0.995 0.99 0.997 -s 1 4 --clean -m OTU

"""

citation = """
This script relies heavily on USEARCH, MUSCLE, and BLAST.
If this script assisted with a publication, please cite the following papers
(or updated citations, depending on the versions of USEARCH, etc., used).

PURC:
-Rothfels, C.J., K.M. Pryer, and F-W. Li. 2016. Next-generation polyploid
phylogenetics: rapid resolution of hybrid polyploid complexes using PacBio
single-molecule sequencing. New Phytologist

USEARCH/UCLUST:
-Edgar, R.C. 2010. Search and clustering orders of magnitude faster than BLAST.
Bioinformatics 26(19), 2460-2461.

UCHIME:
-Edgar, R.C., B.J. Haas, J.C. Clemente, C. Quince, R. Knight. 2011.
UCHIME improves sensitivity and speed of chimera detection, Bioinformatics 27(16), 2194-2200.

Cutadapt:
-Martin, M. 2011. Cutadapt removes adapter sequences from high-throughput sequencing reads.
EMBnet.journal 17:10-12.

MUSCLE:
-Edgar, R.C. 2004. MUSCLE: Multiple sequence alignment with high accuracy and high throughput.
Nucleic Acids Research 32:1792-1797.

BLAST:
-Camacho, C., G. Coulouris, V. Avagyan, N. Ma, J. Papadopoulos, et al. 2009.
BLAST+: Architecture and applications. BMC Bioinformatics 10: 421.
"""


import sys
import os
import re
import subprocess
import glob
import shutil
import time
import datetime
import argparse
import fileinput
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo

def parse_fasta(infile):
	"""Reads in a fasta, returns a list of biopython seq objects"""
	AllSeq = SeqIO.parse(infile, 'fasta')
	return [i for i in AllSeq]

def SplitBy(annotd_seqs_file, split_by = "locus-taxon", Multiplex_perBC_flag=False, Clust_particular=False, locus_to_clust='', specimen_to_clust=''):
	"""Uses the annotated sequences to split sequences into different files based on splits_list
	(could be by barcode or by locus, etc); returns a dictionary of seq counts for each subgroup; Clust_particular to control for specific specimen or locu to cluster"""

	#annotd_seqs = open(annotd_seqs_file, 'r')
	unsplit_seq = parse_fasta(annotd_seqs_file)

	try:
		os.chdir(masterFolder)
	except:
		pass

	splits_file_dic = {}
	splits_count_dic = {}
	LocusTaxonCountDict_unclustd = {}
	splits_list = []

	for each_seq in unsplit_seq:
		#finding the identifing annotation for the split of interest.
		# e.g., BC01, BC02, ... or GAP, PGI, ...

		if split_by == "taxon":
			split = str(each_seq.id).split('|')[0]
			if Clust_particular:
				if split != specimen_to_clust:
					continue
		elif split_by == "locus":
			split = str(each_seq.id).split('|')[1]
			if Clust_particular:
				if split != locus_to_clust:
					continue
		elif split_by == "taxon-locus":
			split = str(each_seq.id).split('|')[0] + "_" + str(each_seq.id).split('|')[1]
		elif split_by == "locus-taxon": # same as above, but folders/files labeled with locus name before taxon name
			split = str(each_seq.id).split('|')[1] + "_" + str(each_seq.id).split('|')[0]
		elif split_by == "barcode" and Multiplex_perBC_flag:
			split = str(each_seq.id).split('|')[3]
		elif split_by == "group" and Multiplex_perBC_flag: # where group is the let set to identify the taxonomic groups, e.g, A, B, C, ...
			split = str(each_seq.id).split('|')[2]
		elif split_by == "barcode" and not Multiplex_perBC_flag:
			split = str(each_seq.id).split('|')[2]
		else:
			sys.stderr.write('Error: Attempting to split sequences by an inappropriate (absent?) category.\n')

		if split not in splits_list: # add split identifier to the list of splits
			splits_list.append(split)
			os.makedirs(split)
			os.chdir(split)
			seq_file = split + '.fa'
			file_handle = open(seq_file, 'w') #FWL only have to do this once, for the first time that split occurs?
			splits_file_dic[split] = file_handle # e.g., {BC01:BC01.fa, BC02:BC02.fa, ...}
			# os.chdir('..')
		else:
			os.chdir(split)

		file_handle = splits_file_dic[split] #use that split as the key to find the corresponding output file
		file_handle.write('>' + str(each_seq.id) + '\n' + str(each_seq.seq) + '\n') #write the seq

		try:
			splits_count_dic[split] += 1
		except:
			splits_count_dic[split] = 1

		if split_by == "locus":
			try:
				LocusTaxonCountDict_unclustd[str(each_seq.id).split('|')[0], str(each_seq.id).split('|')[1]] += 1  # {('C_dia_5316', 'ApP'): 28} for example
				# The locus names are the same as locus_folder
			except:
				LocusTaxonCountDict_unclustd[str(each_seq.id).split('|')[0], str(each_seq.id).split('|')[1]] = 1

		os.chdir('..')
	if split_by == "locus":
		return splits_count_dic, LocusTaxonCountDict_unclustd #as {'BC01': 150, 'BC02': 156} for example
	else:
		return splits_count_dic

def clusterIt(file, clustID, round, previousClusterToCentroid_dict, verbose_level=1):
	"""The clustering step, using the clustID value"""
	outFile = re.sub(r"(.*)\.fa", r"\1C%s_%s.fa" %(round, clustID), file) # The rs indicate "raw" and thus python's escaping gets turned off
	outClustFile = re.sub(r"(.*).fa", r"\1clusts%s.uc" %(round), file)
	logFile = outFile + ".clusterIt.log"
	if round == 1:
		vsearch_cline = "%s --cluster_fast %s --id %f --gapopen 3I/1E --consout %s --uc %s --sizeout --threads %s --log %s" % (Vsearch, file, clustID, outFile, outClustFile, num_threads, logFile)
		#vsearch_cline = "%s -cluster_smallmem %s -id %f -gapopen 3I/1E -sortedby other -centroids %s -uc %s -sizeout" % (Vsearch, file, clustID, outFile, outClustFile)
		#vsearch_cline = "%s -cluster_smallmem %s -id %f -gapopen 3I/1E -usersort -consout %s -uc %s -sizeout" % (Vsearch, file, clustID, outFile, outClustFile) # Vsearch 7
        # Can add in "-cons_truncate" to the vsearch call, if the primer removal isn't effective, but note some problems with partial sequences results.
	elif round > 1:
		vsearch_cline = "%s --cluster_fast %s --id %f --gapopen 3I/1E --consout %s --uc %s --sizein --sizeout --threads %s --log %s" % (Vsearch, file, clustID, outFile, outClustFile, num_threads, logFile)
		#vsearch_cline = "%s -cluster_smallmem %s -id %f -gapopen 3I/1E -sortedby other -centroids %s -uc %s -sizein -sizeout" % (Vsearch, file, clustID, outFile, outClustFile)
		#vsearch_cline = "%s -cluster_smallmem %s -id %f -gapopen 3I/1E -usersort -consout %s -uc %s -sizein -sizeout" % (Vsearch, file, clustID, outFile, outClustFile) # Vsearch 7
	process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout
	if verbose_level == 2:
		#print '\n**Vsearch-clustering output on', file, '**\n'
		#print err
		log.write('\n**Vsearch-clustering output on' + str(file) + '**\n')
		log.write(str(err))

	# remove "centroid=" instances from output files
	with fileinput.input(files = outFile, inplace = True) as fa:
		for line in fa:
			line = line.strip("\n")
			if line.startswith(">"):
				line = line.replace("centroid=", "")
				print(line)
			else:
				print(line)

	uc = open(outClustFile, 'r')
	#if round == 1:
	ClusterToCentroid_dict = {}
	#	global ClusterToCentroid_dict
	for line in uc:
		if line.startswith('C'):
			cluster_name = 'Cluster' + str(line.split('\t')[1])
			centroid_seq_name = line.split('\t')[-2]

			if round > 1:
				#ClusterToCentroid_dict[cluster_name] = previousClusterToCentroid_dict[centroid_seq_name.split(';')[0]]
				ClusterToCentroid_dict[cluster_name] = centroid_seq_name
			elif round == 1:
				ClusterToCentroid_dict[cluster_name] = centroid_seq_name

	if round == 4:
		try:
			clustered_seqs = parse_fasta(outFile)
			renamed_clustered_seqs = open('temp', 'a')
			for seq in clustered_seqs: # seq as 'Cluster15;size=1;'
				#new_seq_id = ClusterToCentroid_dict[str(seq.id).split(';')[0]]
				new_seq_id = str(seq.id).replace(str(seq.id).split(';')[0], ClusterToCentroid_dict[str(seq.id).split(';')[0]])
				#new_seq_id = ClusterToCentroid_dict[str(seq.id).split(';')[0]] + ';' + str(seq.id).split(';')[0:]
				#print str(seq.id).split(';')
				renamed_clustered_seqs.write('>' + new_seq_id + '\n' + str(seq.seq) + '\n')
			os.remove(outFile)
			os.rename('temp', outFile)
		except:
			log.write(outFile + ' is empty; perhaps sizeThreshold too high\n')

	return outFile, ClusterToCentroid_dict, outClustFile

# This function looks for PCR chimeras -- those formed within a single amplicon pool
def deChimeIt(file, round, abskew=1.9, verbose_level=0):
	"""Chimera remover"""
	outFile = re.sub(r"(.*)\.fa", r"\1dCh%s.fa" %(round), file) # The rs indicate "raw" and thus python's escaping gets turned off
	outFile_uchime = re.sub(r"(.*)\.fa", r"\1dCh%s.uchime" %(round), file) # The rs indicate "raw" and thus python's escaping gets turned off
	logFile = outFile + ".deChimeIt.log"
	vsearch_cline = "%s --uchime_denovo %s --abskew %s --nonchimeras %s --uchimeout %s --log %s" % (Vsearch, file, abskew, outFile, outFile_uchime, logFile)
	process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout
	if verbose_level in [1, 2]:
		log.write('\n**Vchime output on' + str(file) + '**\n')
		log.write(str(err))

	# count number of chimera seq found by parsing 'outFile_uchime'
	chimera_count = 0
	uchime_out = open(outFile_uchime, 'r')
	for line in uchime_out:
		line = line.strip('\n')
		if line.split('\t')[-1] == 'Y':
			chimera_count = chimera_count + 1
	return outFile, chimera_count

def sortIt_length(file, verbose_level=0):
	"""Sorts clusters by seq length"""
	outFile = re.sub(r"(.*)\..*", r"\1_Sl.fa", file) # Make the outfile name by cutting off the extension of the infile name, and adding "_S1.fa"
	vsearch_cline = "%s --sortbylength %s --output %s --threads %s" %(Vsearch, file, outFile, num_threads)
	#vsearch_cline = "%s -sortbylength %s -output %s" %(Vsearch, file, outFile) # Vsearch 7
	process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout
	if verbose_level == 2:
		#print '\n**Vsearch-sorting output on', file, '**\n'
		#print err
		log.write('\n**Vsearch-sorting output on' + str(file) + '**\n')
		log.write(str(err))
	return outFile # having it spit out the outfile name, if necessary, so that it can be used to call downstream stuff and avoid complicated glob.globbing

def sortIt_size(file, thresh, round, verbose_level=0):
	"""Sorts clusters by size, and removes those that are smaller than a particular size
	(sent as thresh -- ie, sizeThreshold).
    "round" is used to annotate the outfile name with S1, S2, etc. depending on which sort this is"""
	outFile = re.sub(r"(.*)\.fa", r"\1Ss%s.fa" %(round), file)
	logFile = outFile + ".sortIt_size.log"
	vsearch_cline = "%s --sortbysize %s --output %s --minsize %d --log %s --threads %s" %(Vsearch, file, outFile, thresh, logFile, num_threads)
	#vsearch_cline = "%s -sortbysize %s -output %s -minsize %f" %(Vsearch, file, outFile, thresh) # Vsearch 7
	process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout
	if verbose_level == 2:
		#print '\n**Vsearch-sorting output on', file, '**\n'
		#print err
		log.write('\n**Vsearch-sorting output on' + str(file) + '**\n')
		log.write(str(err))
	return outFile

def align_and_consensus(inputfile, output_prefix):
	output_alignment = output_prefix.split(';')[0] + '_aligned.fa'
	output_consensus = output_prefix.split(';')[0] + '_consensus.fa'
	mafft_cline = '%s --auto %s > %s' % (Mafft, inputfile, output_alignment)
	process = subprocess.Popen(mafft_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout

	alignment = AlignIO.read(output_alignment, 'fasta')
	summary_align = AlignInfo.SummaryInfo(alignment)
	consensus = summary_align.gap_consensus(ambiguous='N',threshold=0.51)
	output = open(output_consensus, 'w')
	output.write('>' + output_prefix + '\n' + str(consensus).replace('-','') + '\n')
	return



def IterativeClusterDechimera(annotd_seqs_file, clustID, clustID2, clustID3, sizeThreshold, sizeThreshold2, verbose_level = 1): # M_p_barcode was set to FALSE, whcih led to splitting by taxon instead of barcode
	## Split sequences into separate files/folders for each locus ##
	sys.stderr.write('Splitting sequences into a folder/file for each locus...\n')
	if locus_to_cluster is not None: # specific target is given
		locusCounts, LocusTaxonCountDict_unclustd = SplitBy(annotd_seqs_file = annotd_seqs_file, split_by = "locus", Clust_particular = True, locus_to_clust = locus_to_cluster)
	else: # cluster all
		locusCounts, LocusTaxonCountDict_unclustd = SplitBy(annotd_seqs_file = annotd_seqs_file, split_by = "locus", Clust_particular = False, locus_to_clust = '')

	sys.stderr.write('Clustering/dechimera-izing seqs...\n')
	log.write('#Sequence clustering/dechimera-izing#\n')
	all_folders_loci = list(locusCounts.keys()) # SplitBy makes a dictionary where the keys are the subcategories (and thus also the
		# folders) and they correspond to the counts for each.
	LocusTaxonCountDict_clustd = {}
	LocusTaxonCountDict_chimera = {}

	## Go through each locus ##
	for locus_folder in all_folders_loci: # locus_folder = locus name
		os.chdir(locus_folder)
		sys.stderr.write('\nWorking on: ' + locus_folder + '...\n')
		if verbose_level in [1,2]:
			log.write('\nWorking on ' + str(locus_folder) + ' ...\n')

		if not os.stat(locus_folder + ".fa").st_size == 0: # ie, the file is not empty
			## Split sequences into separate taxon folders ##
			if specimen_to_cluster is not None: # specific target is given
				taxonCounts = SplitBy(annotd_seqs_file = locus_folder + ".fa", split_by = "taxon", Clust_particular = True, specimen_to_clust = specimen_to_cluster)
			else: # cluster all
				taxonCounts = SplitBy(annotd_seqs_file = locus_folder + ".fa", split_by = "taxon", Clust_particular = False, specimen_to_clust = '')
			all_folders_taxon = list(taxonCounts.keys())

			for taxon_folder in all_folders_taxon:
				if verbose_level in [1,2]:
					log.write("Working on " + taxon_folder + '\n')
				os.chdir(taxon_folder)

				if verbose_level in [1,2]:
					log.write("\tFirst clustering\n")
					log.write("\nAttempting to sort: " + taxon_folder + ".fa\n")

				sorted_length = sortIt_length(file = taxon_folder + ".fa", verbose_level = verbose_level)
				clustered1, previousClusterToCentroid_dict, outClustFile1 = clusterIt(file = sorted_length, previousClusterToCentroid_dict = '', clustID = clustID, round = 1, verbose_level = verbose_level)

				if verbose_level in [1,2]:
					log.write("\tFirst chimera slaying expedition\n")
				deChimered1, chimera_count1 = deChimeIt(file = clustered1, round = 1, abskew = abskew, verbose_level = verbose_level)

				if verbose_level in [1,2]:
					log.write("\tSecond clustering\n")
				sorted_size1 = sortIt_size(file = deChimered1, thresh = sizeThreshold, round = 1, verbose_level = verbose_level)
				clustered2, previousClusterToCentroid_dict, outClustFile2 = clusterIt(file = sorted_size1, previousClusterToCentroid_dict = previousClusterToCentroid_dict, clustID = clustID2, round = 2, verbose_level = verbose_level)

				if verbose_level in [1,2]:
					log.write("\tSecond chimera slaying expedition\n")
				deChimered2, chimera_count2 = deChimeIt(file = clustered2, round = 2, abskew = abskew, verbose_level = verbose_level)

				if verbose_level in [1,2]:
					log.write("\tThird clustering\n")
				sorted_size2 = sortIt_size(file = deChimered2, thresh = sizeThreshold, round = 2, verbose_level = verbose_level)
				clustered3, previousClusterToCentroid_dict, outClustFile3 = clusterIt(file = sorted_size2, previousClusterToCentroid_dict = previousClusterToCentroid_dict, clustID = clustID3, round = 3, verbose_level = verbose_level)

				if verbose_level in [1,2]:
					log.write("\tThird chimera slaying expedition\n")
				deChimered3, chimera_count3 = deChimeIt(file = clustered3, round = 3, abskew = abskew, verbose_level = verbose_level)

				if verbose_level in [1,2]:
					log.write("\Fourth clustering\n")
				sorted_size3 = sortIt_size(file = deChimered3, thresh = sizeThreshold, round = 2, verbose_level = verbose_level)
				clustered4, previousClusterToCentroid_dict, outClustFile4 = clusterIt(file = sorted_size3, previousClusterToCentroid_dict = previousClusterToCentroid_dict, clustID = clustID3, round = 4, verbose_level = verbose_level)

				if verbose_level in [1,2]:
					log.write("\tThird chimera slaying expedition\n")
				deChimered4, chimera_count4 = deChimeIt(file = clustered4, round = 4, abskew = abskew, verbose_level = verbose_level)

				sorted_size4 = sortIt_size(file = deChimered4, thresh = sizeThreshold2, round = 4, verbose_level = verbose_level)

				### Collect all sequences from each cluster and re-consensus ###
				ClusterToSeq_dict4 = {}
				with open(outClustFile4, 'r') as file_to_read:
					for line in file_to_read:
						if line.startswith("H") or line.startswith("C"):
							key = 'Cluster' + line.split('\t')[1]
							seq = line.split("\t")[8].split(";")[0]
							try:
								ClusterToSeq_dict4[key].append(seq)
							except:
								ClusterToSeq_dict4[key] = [seq]

				with open(outClustFile3, 'r') as file_to_read:
					for line in file_to_read:
						if line.startswith("H"):
							seq = line.split("\t")[8].split(";")[0]
							centroid = line.strip("\n").split("\t")[9].split(";")[0]
							for key in ClusterToSeq_dict4.keys():
								if centroid in ClusterToSeq_dict4[key]:
									ClusterToSeq_dict4[key].append(seq)

				with open(outClustFile2, 'r') as file_to_read:
					for line in file_to_read:
						if line.startswith("H"):
							seq = line.split("\t")[8].split(";")[0]
							centroid = line.strip("\n").split("\t")[9].split(";")[0]
							for key in ClusterToSeq_dict4.keys():
								if centroid in ClusterToSeq_dict4[key]:
									ClusterToSeq_dict4[key].append(seq)

				with open(outClustFile1, 'r') as file_to_read:
					for line in file_to_read:
						if line.startswith("H"):
							seq = line.split("\t")[8].split(";")[0]
							centroid = line.strip("\n").split("\t")[9].split(";")[0]
							for key in ClusterToSeq_dict4.keys():
								if centroid in ClusterToSeq_dict4[key]:
									ClusterToSeq_dict4[key].append(seq)

				# Go through the first clustering uc file
				#ClusterToSeq_dict1 = {}
				#file_to_read = open(outClustFile1, 'r')
				#for line in file_to_read:
				#	line = line.strip('\n')
				#	if line.startswith('H') or line.startswith('C'):
				#		key = 'Cluster' + line.split('\t')[1]
				#		seq = line.split('\t')[8]
				#		try:
				#			ClusterToSeq_dict1[key].append(seq)
				#		except:
				#			ClusterToSeq_dict1[key] = [seq]
				#file_to_read.close()

				# Go through the second clustering uc file
				#ClusterToSeq_dict2 = {}
				#file_to_read = open(outClustFile2, 'r')
				#for line in file_to_read:
				#	line = line.strip('\n')
				#	if line.startswith('H') or line.startswith('C'):
				#		key = 'Cluster' + line.split('\t')[1]
				#		seqs = ClusterToSeq_dict1[line.split('\t')[8].split(';')[0]] # use Cluster1 as key
				#		for seq in seqs:
				#			try:
				#				ClusterToSeq_dict2[key].append(seq)
				#			except:
				#				ClusterToSeq_dict2[key] = [seq]
				#file_to_read.close()

				# Go through the third clustering uc file
				#ClusterToSeq_dict3 = {}
				#file_to_read = open(outClustFile3, 'r')
				#for line in file_to_read:
				#	line = line.strip('\n')
				#	if line.startswith('H') or line.startswith('C'):
				#		key = 'Cluster' + line.split('\t')[1]
				#		seqs = ClusterToSeq_dict2[line.split('\t')[8].split(';')[0]] # use Cluster1 as key
				#		for seq in seqs:
				#			try:
				#				ClusterToSeq_dict3[key].append(seq)
				#			except:
				#				ClusterToSeq_dict3[key] = [seq]
				#file_to_read.close()

				# Go through the forth clustering uc file
				#ClusterToSeq_dict4 = {}
				#file_to_read = open(outClustFile4, 'r')
				#for line in file_to_read:
				#	line = line.strip('\n')
				#	if line.startswith('H') or line.startswith('C'):
				#		key = 'Cluster' + line.split('\t')[1]
				#		seqs = ClusterToSeq_dict3[line.split('\t')[8].split(';')[0]] # use Cluster1 as key
				#		for seq in seqs:
				#			try:
				#				ClusterToSeq_dict4[key].append(seq)
				#			except:
				#				ClusterToSeq_dict4[key] = [seq]
				#file_to_read.close()

				## Align and re-consensus of all constituent seq for each cluster ##
				SeqDict = SeqIO.index(taxon_folder + ".fa", 'fasta')
				for each_cluster in ClusterToSeq_dict4:
					if len(ClusterToSeq_dict4[each_cluster]) >= int(sizeThreshold2):
						cluster_seq_file = open(each_cluster, 'w')
						for seq in ClusterToSeq_dict4[each_cluster]:
							cluster_seq_file.write('>' + seq + '\n' + str(SeqDict[seq].seq) + '\n')
						cluster_seq_file.close()
						new_seq_name = taxon_folder + '_' + each_cluster + ';size=' + str(len(ClusterToSeq_dict4[each_cluster])) + ';'
						align_and_consensus(each_cluster, new_seq_name)

				## Put all consensus seq into one file *_Cluster_Finalconsensus.fa ##
				all_consensus_seq = open(taxon_folder + '_Cluster_Finalconsensus.fa', 'w')
				files_to_add_reconsensus = glob.glob("*_consensus.fa")
				for file in files_to_add_reconsensus:
					shutil.copyfileobj(open(file,'r'), all_consensus_seq) #Add each file to the final output
				all_consensus_seq.close()

				## Do final clustering and chimera-killing ##
				vsearch_cline = "%s --sortbysize %s --output %s" %(Vsearch, taxon_folder + '_Cluster_Finalconsensus.fa', taxon_folder + '_Cluster_FinalconsensusSs.fa')
				process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
				(out, err) = process.communicate() #the stdout and stderr

				vsearch_cline = "%s --cluster_fast %s --id %f --gapopen 3I/1E --consout %s --uc %s --sizein --sizeout" % (Vsearch, taxon_folder + '_Cluster_FinalconsensusSs.fa', clustID4, taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + '.fa', taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + '.uc')
				process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
				(out, err) = process.communicate() #the stdout and stderr

				vsearch_cline = "%s --uchime_denovo %s --abskew %s --nonchimeras %s --uchimeout %s" % (Vsearch, taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + '.fa', abskew, taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.fa', taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.uchime')
				process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
				(out, err) = process.communicate() #the stdout and stderr

				## Count number of chimera seq found by parsing 'outFile_uchime' ##
				chimera_count5 = 0
				uchime_out = open(taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.uchime', 'r')
				for line in uchime_out:
					line = line.strip('\n')
					if line.split('\t')[-1] == 'Y':
						chimera_count5 = chimera_count5 + 1
				uchime_out.close()
				LocusTaxonCountDict_chimera[taxon_folder, locus_folder] = [chimera_count1, chimera_count2, chimera_count3, chimera_count4, chimera_count5]

				## Count clustered seq and store in LocusTaxonCountDict_clustd as {('C_dia_5316', 'ApP'): 28} for example ##
				try:
					clustered_seq_file = parse_fasta(taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.fa')
					for each_seq in clustered_seq_file:
						try:
							LocusTaxonCountDict_clustd[taxon_folder, locus_folder] += 1  # {('C_dia_5316', 'ApP'): 28} for example
						except:
							LocusTaxonCountDict_clustd[taxon_folder, locus_folder] = 1
				except:
					if verbose_level in [1,2]:
						log.write(str(all_consensus_seq) + 'is an empty file\n')

				## Rename sequences in the final fasta: add taxon name ##
				with open(taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.fa', "r") as renamingFile:
					with open(taxon_folder + '_OTUs.fa', "w") as renamedFile:
						for line in renamingFile:
							if line.startswith(">"):
								line = re.sub(r"centroid=","",line)
								line = re.sub(r";seqs=\d+;", ";", line)
							renamedFile.write(line)
				#sed_cmd = "sed 's/>/>%s_/g' %s > %s" % (taxon_folder, taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh.fa', taxon_folder + '_Cluster_FinalconsensusSsC' + str(clustID4) + 'dCh_renamed.fa')
				#process = subprocess.Popen(sed_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
				#(out, err) = process.communicate() #the stdout and stderr

				## Remove intermediate files if requsted ##
				if args.clean:
					filt_to_remove = list(set(glob.glob('*')) - set(glob.glob('*_OTUs.fa')) - set(glob.glob(taxon_folder+'.fa')))
					for file in filt_to_remove:
						os.remove(file)

				os.chdir("..") # To get out of the current barcode folder and ready for the next one
		os.chdir("..") # To get out of the current locus folder and ready for the next one
	log.write('\t...done\n\n')

	## Put all the sequences together ##
	sys.stderr.write('\rPutting all the sequences together......\n\n')
	for locus_folder in all_folders_loci: # Looping through each of the locus folders
		#outputfile_name = str(locus_folder) + '_clustered.txt' # "locus_folder" is also the name of the locus
		#outputfile = open(outputfile_name, 'w')
		outputfile_name_reconsensus = str(locus_folder) + '_OTUs.fa'
		outputfile_reconsensus = open(outputfile_name_reconsensus, 'w')

		os.chdir(locus_folder)
		taxonForThisLocus = glob.glob("*")
		for taxon_folder in taxonForThisLocus: # have to go into each barcode folder in each locus folder
			if os.path.isdir(taxon_folder): # the glob might have found some files as well as folders
				os.chdir(taxon_folder)
				#files_to_add = glob.glob("*Ss4.fa")
				#for file in files_to_add:
				#	shutil.copyfileobj(open(file,'rb'), outputfile) #Add each file to the final output
				files_to_add_reconsensus = glob.glob("*_OTUs.fa")
				for file in files_to_add_reconsensus:
					shutil.copyfileobj(open(file,'rb'), outputfile_reconsensus) #Add each file to the final output

				os.chdir('..')
		os.chdir('..')
		outputfile_reconsensus.close()

	return LocusTaxonCountDict_clustd, LocusTaxonCountDict_unclustd, LocusTaxonCountDict_chimera

def mafftIt(file, verbose_level=0):
	"""Aligns the sequences using MUSCLE"""
	outFile = re.sub(r"(.*)\..*", r"\1.aligned.fa", file) # The rs indicate "raw" and thus python's escaping gets turned off
	mafft_cline = '%s --auto %s > %s' % (Mafft, file, outFile)
	process = subprocess.Popen(mafft_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(out, err) = process.communicate() #the stdout and stderr
	savestdout = sys.stdout
	if verbose_level == 2:
		#print '\n**Mafft output on', file, '**\n'
		#print err
		log.write('\n**Mafft output on' + str(file) + '**\n')
		log.write(str(err))
	return outFile

def dada(annotd_seqs_file, raw_fastq_sequences, Forward_primer, Reverse_primer, minLen, maxLen, maxEE, RscriptPath, verbose_level = 1):
	log.write("DADA2\n")
	## Split sequences into separate files/folders for each locus ##
	sys.stderr.write('Splitting sequences into a folder/file for each locus...\n')
	if locus_to_cluster is not None: # specific target is given
		locusCounts, LocusTaxonCountDict_unclustd = SplitBy(annotd_seqs_file = annotd_seqs_file, split_by = "locus", Clust_particular = True, locus_to_clust = locus_to_cluster)
	else: # cluster all
		locusCounts, LocusTaxonCountDict_unclustd = SplitBy(annotd_seqs_file = annotd_seqs_file, split_by = "locus", Clust_particular = False, locus_to_clust = '')

	sys.stderr.write('Sorting sequences into ASVs...\n')
	log.write('#Sorting sequences into ASVs#\n')
	all_folders_loci = list(locusCounts.keys()) # SplitBy makes a dictionary where the keys are the subcategories (and thus also the
		# folders) and they correspond to the counts for each.
	LocusTaxonCountDict_clustd = {} # {('C_dia_5316', 'ApP'): 28} for example, to store clusted seq count
	LocusTaxonCountDict_chimera = {} # {('C_dia_5316', 'ApP'): [1,0,0,0,0]} for example, to store the chimerc seq count for each chimera-killing step

	## Go through each locus ##
	locusCount = 0
	for locus_folder in all_folders_loci: # locus_folder = locus name
		os.chdir(locus_folder)
		sys.stderr.write('\nWorking on: ' + locus_folder + '...\n')
		if verbose_level in [1,2]:
			log.write('\nWorking on ' + str(locus_folder) + ' ...\n')

		if not os.stat(locus_folder + ".fa").st_size == 0: # ie, the file is not empty
			## Split sequences into separate taxon folders ##
			if specimen_to_cluster is not None: # specific target is given
				taxonCounts = SplitBy(annotd_seqs_file = locus_folder + ".fa", split_by = "taxon", Clust_particular = True, specimen_to_clust = specimen_to_cluster)
			else: # cluster all
				taxonCounts = SplitBy(annotd_seqs_file = locus_folder + ".fa", split_by = "taxon", Clust_particular = False, specimen_to_clust = '')
			all_folders_taxon = list(taxonCounts.keys())

			for taxon_folder in all_folders_taxon:
				if verbose_level in [1,2]:
					log.write("Working on " + taxon_folder + '\n')
				os.chdir(taxon_folder)
				subset_fasta_seqs_from_fastq("%s.fa" % taxon_folder, raw_fastq_sequences)
				writeASV(taxon_folder, Forward_primer[locusCount], Reverse_primer[locusCount], minLen, maxLen, maxEE)
				with open("%s_DADA2.log" % taxon_folder, "w") as logfile:
					dadaCMD = "%s %s_DADA2.R" %(RscriptPath, taxon_folder)
					process = subprocess.Popen(dadaCMD, stdout=logfile, stderr=logfile, shell=True, text=True)
					process.communicate()
				## Count ASVs and store in LocusTaxonCountDict_clustd as {('C_dia_5316', 'ApP'): 28} for example ##
				try:
					clustered_seq_file = parse_fasta(taxon_folder + '_ASVs.fa')
					for each_seq in clustered_seq_file:
						try:
							LocusTaxonCountDict_clustd[taxon_folder, locus_folder] += 1  # {('C_dia_5316', 'ApP'): 28} for example
						except:
							LocusTaxonCountDict_clustd[taxon_folder, locus_folder] = 1
				except:
					print("WARNING: No ASVs found for %s" % taxon_folder)
					log.write("WARNING: No ASVs found for %s" % taxon_folder)

				## Count chimeras
				try:
					chimera_seq_file = parse_fasta(taxon_folder + '_chimeras.fa')
					for each_seq in chimera_seq_file:
						try:
							LocusTaxonCountDict_chimera[taxon_folder, locus_folder] += 1  # {('C_dia_5316', 'ApP'): 28} for example
						except:
							LocusTaxonCountDict_chimera[taxon_folder, locus_folder] = 1
				except:
					if verbose_level in [2]:
						log.write("No chimeras detected in %s\n" % taxon_folder)
					LocusTaxonCountDict_chimera[taxon_folder, locus_folder] = 0
				os.chdir("..")
		locusCount += 1
		os.chdir("..")
	log.write('\t...done\n\n')

	## Put all the sequences together ##
	sys.stderr.write('\n\nPutting all the sequences together...\n\n')
	for locus_folder in all_folders_loci: # Looping through each of the locus folders
		outputfile_name_reconsensus = Output_prefix + '_4_' + str(locus_folder) + '_ASVs.fa'
		outputfile_reconsensus = open(outputfile_name_reconsensus, 'w')

		os.chdir(locus_folder)
		taxonForThisLocus = glob.glob("*")
		for taxon_folder in taxonForThisLocus: # have to go into each barcode folder in each locus folder
			if os.path.isdir(taxon_folder): # the glob might have found some files as well as folders
				os.chdir(taxon_folder)
				files_to_add_reconsensus = glob.glob("*_ASVs.fa")
				for file in files_to_add_reconsensus:
					shutil.copyfileobj(open(file,'r'), outputfile_reconsensus) #Add each file to the final output

				os.chdir('..')
		os.chdir('..')
		outputfile_reconsensus.close()

	return LocusTaxonCountDict_clustd, LocusTaxonCountDict_chimera

###### RUN ######

# Parse the arguments #
parser = argparse.ArgumentParser(description=usage)

required = parser.add_argument_group('required arguments')
required.add_argument('-f','--annotated_file', type=str, action="store", metavar='\b', required=True,
               help='The annotated fasta file, generated by purc')
required.add_argument('-o','--output_folder', type=str, action="store", metavar='\b', required=True,
               help='The output folder')
optional.add_argument('-c','--clustering_identities', type=float, nargs=4, metavar='\b',
               help='The similarity criterion for the first, second, third and fourth USEARCH clustering',
               default=[0.997, 0.995, 0.99, 0.997])
optional.add_argument('-s','--size_threshold', type=int, nargs=2, metavar='\b',
                    help='The min. number of sequences/cluster necessary for that cluster to be retained (set to 2 to remove singletons, 3 to remove singletons and doubles, etc)',
                    default=[1, 4])
optional.add_argument('-m', '--method', type=str, action="store", metavar="\b", help='Method to use for clustering sequences (OTU or ASV). Reserved for future use',default='OTU')
optional.add_argument('--minimum_length', type=int, action="store",metavar="\b", help="Minimum read length for ASV inference (0 = autodetect)", default=0)
optional.add_argument('--maximum_length', type=int, action="store",metavar="\b", help="Maximum read length for ASV inference (0 = autodetect)", default=0)
optional.add_argument('--maximum_errors', type=int, action="store",metavar="\b", help="Maximum number of expected errors allows in a read for ASV inference", default=5)

optional = parser.add_argument_group('optional arguments')
optional.add_argument('-a','--abundance_skew', type=float, nargs=1, metavar='\b',
                    help='The parameter to control chimera-killing; the default is 1.9',
                    default=[1.9])
optional.add_argument('-l','--locus_to_cluster', type=str, action="store", metavar='\b',
               		help='Only do clustering on this particular locus; e.g. -l APP')
optional.add_argument('-i','--specimen_to_cluster', type=str, action="store", metavar='\b',
               		help='Only do clustering on this particular specimen; e.g. -i Cystopteris_sp_7635')
optional.add_argument('-t','--thread_number', type=int, nargs=1, metavar='\b',
                    help='The number of threads for clustering; the default is 1',
                    default=[1])
optional.add_argument('--clean', action="store_true",
                    help='Remove the intermediate files')

args = parser.parse_args()

annotated_file = args.annotated_file
masterFolder = args.output_folder
clusteringMethod = args.method
clustID = args.clustering_identities[0]
clustID2 = args.clustering_identities[1]
clustID3 = args.clustering_identities[2]
clustID4 = args.clustering_identities[3]
sizeThreshold = args.size_threshold[0]
sizeThreshold2 = args.size_threshold[1]
abskew = args.abundance_skew[0]
num_threads = args.thread_number[0]
specimen_to_cluster = args.specimen_to_cluster
locus_to_cluster = args.locus_to_cluster

#if clusteringMethod == "OTU" and clustID not in globals() or clustID2 not in globals() or clustID3 not in globals() or clustID4 not in globals() or sizeThreshold not in globals() or sizeThreshold2 not in globals():
#	print("ERROR: OTU clustering requires all entries in --clustering_identities and --size_threshold to be set")
#	sys.exit(usage)

#purc_location = os.path.dirname(os.path.abspath( __file__ ))
Vsearch = 'vsearch'
Mafft = 'mafft'

# Check if dependencies are in place
sys.stderr.write('Checking dependencies...\n')
#if not os.path.isfile(Vsearch):
#	sys.exit("Error: couldn't find the Vsearch executable")
#if not os.path.isfile(Mafft):
#	sys.exit("Error: couldn't find the Mafft executable")

# Check if mafft can be executed
#if not os.access(Mafft, os.R_OK):
	#print("Error: count not execute mafft")
	#sys.exit(1)
#mafft_cline = '%s -version' % (Mafft)
#process = subprocess.Popen(mafft_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#(out, err) = process.communicate() #the stdout and stderr
#if not str(out).startswith('MUSCLE'):
#	sys.exit("Error: could not execute Mafft")

# Check if Vsearch can be executed
#if not os.access(Vsearch, os.R_OK):
	#print("Error: count not execute vsearch")
	#sys.exit(1)
#vsearch_cline = '%s -version' % (Vsearch)
#process = subprocess.Popen(vsearch_cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#(out, err) = process.communicate() #the stdout and stderr
#if not str(out).startswith('vsearch'):
#	sys.exit("Error: could not execute Vsearch")

## Make output folder ##
if os.path.exists(masterFolder): # overwrite existing folder
	shutil.rmtree(masterFolder)
os.makedirs(masterFolder)

ts = time.time()
time_stamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')
log_file = 'purc_log_' + time_stamp + '.txt'
log = open(log_file, 'w')
log.write(logo + '\n')
log.write('purc_recluster.py ' + annotated_file + ' ' + masterFolder + ' ' + str(clustID) + ' ' + str(clustID2) + ' ' + str(clustID3) + ' ' + str(sizeThreshold) + ' ' + str(sizeThreshold2) + '\n\n')
if not os.path.isfile(annotated_file):
	sys.exit("Error: could not find " + annotated_file)

## Recluster and redechimera ##

#print "annotated_file is ", annotated_file
if clusteringMethod == "OTU":
	LocusTaxonCountDict_clustd, LocusTaxonCountDict_unclustd, LocusTaxonCountDict_chimera = IterativeClusterDechimera(annotated_file, clustID, clustID2, clustID3, sizeThreshold, sizeThreshold2)
elif clusteringMethod == "ASV":
	LocusTaxonCountDict_clustd, LocusTaxonCountDict_chimera = dada(annoFileName, fastq_sequences, Forward_primer, Reverse_primer, minLen, maxLen, maxEE, RscriptPath)
else:
	print("ERROR: Clustering method not recognized. Must be OTU or ASV")
	sys.exit(1)

taxon_list = []
locus_list = []
for taxon_locus in list(LocusTaxonCountDict_clustd.keys()):
	taxon_list.append(taxon_locus[0])
	locus_list.append(taxon_locus[1])
taxon_list = set(taxon_list)
locus_list = set(locus_list)

## Producing a summary ##
count_output = open('purc_cluster_counts.xls', 'w')

# Output read count per taxon per locus
log.write("\n**Raw reads per accession per locus**\n")
count_output.write('\n**Raw reads per accession per locus**\n')
count_output.write('\t' + '\t'.join(locus_list) + '\n')
log.write('\t' + '\t'.join(locus_list) + '\n')
for each_taxon in set(taxon_list):
	count_output.write(each_taxon + '\t')
	log.write(each_taxon + '\t')
	for each_locus in locus_list:
		try:
			count_output.write(str(LocusTaxonCountDict_unclustd[each_taxon, each_locus]) + '\t')
			log.write(str(LocusTaxonCountDict_unclustd[each_taxon, each_locus]) + '\t')
		except:
			count_output.write('0\t')
			log.write('0\t')
	count_output.write('\n')
	log.write('\n')

# Output clustered seq count
count_output.write('\n**Final clustered sequences per accession per locus**\n')
log.write('\n**Final clustered sequences per accession per locus**\n')
count_output.write('\t' + '\t'.join(locus_list) + '\n')
log.write('\t' + '\t'.join(locus_list) + '\n')
for each_taxon in set(taxon_list):
	count_output.write(each_taxon + '\t')
	log.write(each_taxon + '\t')
	for each_locus in locus_list:
		try:
			count_output.write(str(LocusTaxonCountDict_clustd[each_taxon, each_locus]) + '\t')
			log.write(str(LocusTaxonCountDict_clustd[each_taxon, each_locus]) + '\t')
		except:
			count_output.write('0\t')
			log.write('0\t')
	count_output.write('\n')
	log.write('\n')

# Output cluster count per locus
count_output.write('\n**Allele/copy/cluster/whatever count by locus**\n')
log.write('\n**Allele/copy/cluster/whatever count by locus**\n')

for each_locus in locus_list:
	file_name = str(each_locus) + '_clustered_reconsensus.fa'
	try:
		seq_no = len(parse_fasta(file_name))
		count_output.write(str(each_locus) + '\t' + str(seq_no) + '\n')
		log.write(str(each_locus) + '\t' + str(seq_no) + '\n')
	except:
		count_output.write(str(each_locus) + '\t0\n')
		log.write(str(each_locus) + '\t0\n')

# Output chimeric seq count
count_output.write('\n**Chimeric clusters/sequence count by locus**\n')
for each_locus in locus_list:
	count_output.write(each_locus + '\n')
	#log.write(each_locus + '\t')
	for each_taxon in set(taxon_list):
		try:
			count_output.write('\t' + each_taxon + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][0]) + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][1]) + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][2]) + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][3]) + '\t' + str(LocusTaxonCountDict_chimera[each_taxon, each_locus][4]))
		except:
			count_output.write('\t' + each_taxon)

		count_output.write('\n')
	#log.write('\n')

## Aligning the sequences ##
fastas = glob.glob("*_OTUs.fa")
for file in fastas:
	sys.stderr.write("Aligning " + file + "\n")
	log.write("Aligning " + file + "\n")
	outFile = mafftIt(file)

fastas = glob.glob("*_clustered_reconsensus.fa")
for file in fastas:
	sys.stderr.write("Aligning " + file + "\n")
	log.write("Aligning " + file + "\n")
	outFile = mafftIt(file)
