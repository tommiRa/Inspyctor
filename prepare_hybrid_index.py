#! /usr/bin/env python
import subprocess
import argparse
import os
import shutil 

parser = argparse.ArgumentParser(description='Generate index files for BWA')
parser.add_argument('-g','--genome', dest="genome", help="The genomic sequence in fasta format")
parser.add_argument('-i','--insert', dest="insert", help="The insert sequence in fasta format")
parser.add_argument('-h','--insert_header', dest="insert_header", help="The header used to identify the inserted sequence")
parser.add_argument('-o','--outputdir', dest="outputdir", help="The output directory")
args = parser.parse_args()

def main():
	# Output directory 
	outdir = args.outputdir
	# Merge the fasta files 
	if os.path.isdir(outdir) == False:
		os.mkdir(outdir)
	else:
		pass
	print("Changing the fasta header for insert fasta file")
	insert_fasta = open(args.plasmid)
	temp_fasta = open("temp_insert.fa","w")
	for i,ln in enumerate(insert_fasta.readlines()):
		if i == 0:
			temp_fasta.write(">{0}".format(args.insert_header))
		else:
			temp_fasta.write(ln)
	insert_fasta.close()
	temp_fasta.close()
	print("Merging fasta files")
	proc = subprocess.Popen(["cat", args.genome, "temp_insert.fa"], stdout = open(outdir + "/hybrid_genome.fa", "w"))
	proc.wait()
	print("Finished merging fasta files")
	# Call BWA to create the index 
	print("Creating index")
	proc2 = subprocess.Popen(["bwa", "index", "-a", "bwtsw", outdir + "/hybrid_genome.fa"])
	proc2.wait()
	print("Cleaning up")
	os.remove("temp_insert.fa")
	print("Finished!")

if __name__ == "__main__":
	main()
