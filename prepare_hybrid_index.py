#! /usr/bin/env python
import subprocess
import argparse
import os
import shutil 

parser = argparse.ArgumentParser(description='Generate index files for BWA')
parser.add_argument('-g','--genome', dest="genome", help="The genomic sequence in fasta format")
parser.add_argument('-p','--plasmid', dest="plasmid", help="The plasmid sequence in fasta format")
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
	print("Merging fasta files")
	proc = subprocess.Popen(["cat", args.genome, args.plasmid], stdout = open(outdir + "/hybrid_genome.fa", "w"))
	proc.wait()
	print("Finished merging fasta files")
	# Call BWA to create the index 
	print("Creating index")
	proc2 = subprocess.Popen(["bwa", "index", "-a", "bwtsw", outdir + "/hybrid_genome.fa"])
	proc2.wait()
	print("Finished!")

if __name__ == "__main__":
	main()
