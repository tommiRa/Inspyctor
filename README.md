# Inspyctor toolkit for plasmid insertion discovery
The following instruction describe an example workflow for detection of plasmid insertion sites. 

## Prerequisites
The following software need to be installed and added to path:
1. [Samtools](https://github.com/samtools/samtools)    
2. [BWA](http://bio-bwa.sourceforge.net/)
3. [PICARD](https://broadinstitute.github.io/picard/)

## Preparation of composite genome index for read alignment 
Before the read alignment a composite index containing the chromosomes and the insert sequence must be prepared. Inspyctor toolkit provides a wrapper script for generating the hybrid index files for BWA. 
```
python prepare_hybrid_index.py -g genome.fa -p insert.fa -h insert_header -o bwa_index_dir
```
## Read alignment with BWA
Reads are aligned using BWA mem (without -M option). The following command aligns the reads using 4 cores and pipes the output to samtools wich converts the alignments to bam-format:
```
bwa mem -t 4 bwa_index_dir reads_pair1.fq reads_pair2.fq | samtools view -bS - > alignment.bam
```
## Removal of unmapped reads, supplementary and secondary alignments 
```
samtools view -h -b -F 2316 alignment.bam > filtered.bam
```
## Sort alignments and marking duplicates 
Before marking the duplicates with PICARD the alignments need to be sorted according to the chromosomal coordinates
```
samtools sort filtered.bam > filtered.sorted.bam
```
Finally the duplicates are marked using PICARD
```
java -Xmx24g -jar picard.jar MarkDuplicates I=filtered.sorted.bam O=final_dp_marked.bam M=duplication_metrix
```
## Detection of plasmid insertion sites  
We need to indicate the insert header (-h flag) which refers to the same identifier we set earlier when preparing the composite genome index. In this example we want to omit the duplicates by adding the -d flag. The script also requires the length (bases) of the inserted sequence.
```
python find_insertion_sites.py -h insert_header -i filtered.sorted.bam -d -l 6201 -o output_directory
```
