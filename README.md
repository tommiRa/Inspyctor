# Inspyctor toolkit for plasmid insertion discovery
The following instruction describe an example workflow for detection of plasmid insertion sites. 

## Prerequisites
Samtools, BWA and PICARD need to be installed and the paths to binaries have to be added to PATH variable. Samtools can be found from : https://github.com/samtools/samtools and BWA can be found from : http://bio-bwa.sourceforge.net/

## Preparation of hybrid index for read alignment 
Before the read alignment a hybrid index containing the chromosomes and the plasmid sequence must be prepared. Inspyctor toolkit provides a wrapper script for generating the hybrid index files for BWA. 
```
python prepare_hybrid_index.py -g genome.fa -p plasmid.fa -o bwa_index_dir
```
## Read alignment with BWA
Reads are aligned using BWA mem (without -M option). The following command aligns the reads using 4 cores and pipes the output to samtools wich converts the alignments to bam-format:
```
bwa mem -t 4 bwa_index_dir reads_pair1.fq reads_pair2.fq | samtools view -bS - > alignment.bam
```
## Removal of unmapped reads and secondary alignments 
```
samtools view -h -b -F 0x4 alignment.bam |samtools view -h -b -F 0x8 - > unmapped_removed.bam
samtools view -b -h -F 0x800 unmapped_removed.bam |samtools view -b -h -F 0x100 - > final.bam
```
## Sort alignments and marking duplicates 
Before marking the duplicates with PICARD the alignments need to be sorted according to the chromosomal coordinates
```
samtools sort final.bam > final.sorted.bam
```
Finally the duplicates are marked using PICARD
```
java -Xmx24g -jar picard.jar MarkDuplicates I=final_sorted.bam O=final_dp_markded.bam M=duplication_metrix
```
## Detection of plasmid insertion sites  

