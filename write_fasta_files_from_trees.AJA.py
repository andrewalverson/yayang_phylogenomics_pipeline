"""
Read a directory of subtrees and write out the corresponding FASTA files
Modified by AJA to read individual FASTA files rather than a concatenated FASTA file (Feb 2019)
"""

import sys,os,newick3,phylo3
import tree_utils
import re
from seq import read_fasta_file

def main(alignDIR,treDIR,tree_file_ending,outDIR):
	if alignDIR[-1] != "/": alignDIR += "/"
	if treDIR[-1] != "/": treDIR += "/"
	if outDIR[-1] != "/": outDIR += "/"

	filecount = 0

	for i in os.listdir(treDIR):
		if i.endswith(tree_file_ending):
			# print i
			filecount += 1
			with open(treDIR+i,"r")as infile:
				intree = newick3.parse(infile.readline())
			clusterID = tree_utils.get_clusterID(i)

			# create FASTA file name
			fasta = re.sub("_\d+$", "", clusterID)
			fasta += ".fasta"
			fasta = alignDIR+fasta
			print i, clusterID, fasta

			seqDICT = {} #key is seqID, value is seq

			# open and parse FASTA file
			for s in read_fasta_file(fasta):
				seqDICT[s.name] = s.seq

			# name of output file
			outname = outDIR+clusterID+"_rr.fa"

			with open(outname,"w") as outfile:
				for label in tree_utils.get_front_labels(intree):
					# added replace to get rid of single quotes, which are in the tree labels but not the alignment labels
					outfile.write(">"+label.replace("'","")+"\n"+seqDICT[label.replace("'","")]+"\n")

	assert filecount > 0,\
		"No file ends with "+tree_file_ending+" found in "+treDIR

if __name__ =="__main__":
	if len(sys.argv) != 5:
		print "usage: python write_fasta_files_from_trees.py alignDIR treDIR tree_file_ending outDIR"
		sys.exit()
	
	alignDIR,treDIR,tree_file_ending,outDIR = sys.argv[1:]
	main(alignDIR,treDIR,tree_file_ending,outDIR)

