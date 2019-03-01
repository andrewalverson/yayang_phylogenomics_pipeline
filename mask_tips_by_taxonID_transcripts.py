"""
Input is a dir of trees that end with ".tt". Change INTREE_FILE_ENDING if needed 

Mask both mono- and paraphyletic tips that belong to the same taxon
If only mask monophyletic tips, comment out this line:
curroot = paraphyly_masking(curroot,unamb_chrDICT)
Keep the tip that has the most un-ambiguous, well-aligned charactors in the trimmed alignment
"""

import newick3,phylo3,os,sys
from tree_utils import get_name,remove_kink
from seq import read_fasta_file

INTREE_FILE_ENDING = ".tt"

def get_clusterID(filename):
	"""given a file name return the cluster id"""
	return filename.split(".")[0]
	
def mask_monophyletic_tips(curroot,unamb_chrDICT):
	going = True
	while going and len(curroot.leaves()) >= 4:
		going = False
		for node in curroot.iternodes(): # walk through nodes
			if not node.istip: continue  # only look at tips
			for sister in node.get_sisters():
				# parse taxon names from the tip labels the check whether taxon names are the same (i.e., sister tips from the same taxon)
				if sister.istip and get_name(node.label)==get_name(sister.label): # masking
					# print node.label,unamb_chrDICT[node.label],sister.label,unamb_chrDICT[sister.label]
					if unamb_chrDICT[node.label.replace("'","")] > unamb_chrDICT[sister.label.replace("'","")]: # added code to remove single quotes from labels
						node = sister.prune()			
					else: node = node.prune()
					if len(curroot.leaves()) >= 4:
						if (node==curroot and node.nchildren==2) or (node!=curroot and node.nchildren==1):
							node,curroot = remove_kink(node,curroot)
					going = True
					break
	return curroot
	
def mask_paraphyletic_tips(curroot,unamb_chrDICT):
	going = True
	while going and len(curroot.leaves()) >= 4:
		going = False
		for node in curroot.iternodes(): #walk through nodes
			if not node.istip: continue #only look at tips
			parent = node.parent
			if node == curroot or parent == curroot:
				continue #no paraphyletic tips for the root
			for para in parent.get_sisters():
				if para.istip and get_name(node.label)==get_name(para.label):
					if unamb_chrDICT[node.label.replace("'","")] > unamb_chrDICT[para.label.replace("'","")]: # added code to remove single quotes from labels
						node = para.prune()
					else: node = node.prune()
					if len(curroot.leaves()) >= 4:
						if (node==curroot and node.nchildren==2) or (node!=curroot and node.nchildren==1):
							node,curroot = remove_kink(node,curroot)
					going = True
					break
	return curroot
	
def main(treDIR,clnDIR,para,intree_file_ending=INTREE_FILE_ENDING):
	if treDIR[-1] != "/": treDIR += "/"
	if clnDIR[-1] != "/": clnDIR += "/"
	assert para=="y" or para=="n", "mask paraphyletic tips? (y/n)"
	mask_para = True if para == "y" else False
	filecount = 0
	
	filematch = {} # key = clusterID, value = the .aln-cln file
	for i in os.listdir(clnDIR):
		if i.endswith(".aln-cln"):
			clusterID = get_clusterID(i)
			assert clusterID not in filematch, \
				"The clusterID "+clusterID+" repeats in "+clnDIR
			filematch[clusterID] = i
			
	# loop over all the tree files
	for i in os.listdir(treDIR):
		# only deal with files that have the right file extension
		if i.endswith(intree_file_ending):
			# print filename for this tree
			print i

			# get clusterID, which is the base file name that connects the tree file with its alignment
			clusterID = get_clusterID(i)
			print clusterID

			# open and parse tree file
			with open(treDIR+i,"r") as infile:
				intree = newick3.parse(infile.readline())
			
			filecount += 1

			# key = seqid (the entire tip label = taxon ID + sequence info), value = number of unambiguous chrs
			chrDICT = {} 

			# read FASTA file for this tree
			for s in read_fasta_file(clnDIR+filematch[clusterID]):
				# remove non-sequence characters
				for ch in ['-','X',"x","?","*"]:
					s.seq = s.seq.replace(ch,"") #ignore gaps, xs and Xs
				chrDICT[s.name] = len(s.seq)
	# get rid of single quotes surrounding taxon name - AJA added
	

			curroot = mask_monophyletic_tips(intree,chrDICT)
			if mask_para: curroot = mask_paraphyletic_tips(curroot,chrDICT)
			with open(treDIR+i+".mm","w") as outfile:
				outfile.write(newick3.tostring(curroot)+";\n")

	assert filecount > 0, \
		"No file ends with "+intree_file_ending+" found in "+treDIR
	
if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "python mask_tips_by_taxonID_transcripts.py treDIR aln-clnDIR mask_paraphyletic(y/n)"
		sys.exit(0)
		
	treDIR,clnDIR,para = sys.argv[1:]
	main(treDIR,clnDIR,para)
