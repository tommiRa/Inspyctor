#! /usr/bin/env python
import subprocess
from operator import itemgetter
import argparse
import numpy as np
import re
import os

parser = argparse.ArgumentParser(description='Find potential insertion sites')
parser.add_argument('-i','--inputBAM', dest="inputBAM")
parser.add_argument('-a','--insert_header', dest="insertID")
parser.add_argument('-l','--insert_length', dest="insertLength")
parser.add_argument('-d','--ignore_duplicates', dest="iduplicates", action="store_true")
parser.add_argument('-m','--minimum_read_count',type=int, dest="min_readcount", default = 3)
parser.add_argument('-o','--outdir', dest="outdir")
args = parser.parse_args()

# Sort based on the read name such that the paired reads will be
# in adjacent rows
def sortN(inputbam, outdir):
    proc = subprocess.Popen(['samtools', 'sort','-O','sam', '-n', inputbam], stdout=open(outdir + "/ntest.sam","w"))
    proc.wait()
    proc2 =  subprocess.Popen(['samtools', 'view', outdir + "/ntest.sam"], stdout=open(outdir + "/ntest_nohead.sam","w"))
    proc2.wait()
    return("finished")

# Determines if read is a PCR or optical duplicate
def isDuplicate(flag):
    binarized_flag = bin(int(flag)).replace("0b","")
    if len(binarized_flag) >= 11:
        if binarized_flag[-11] == "1":
            return(True)
        else:
            return(False)
    else:
        return(False)

def parseChimericAlignments(cols, insertID, insert_match):

    # Find if SA-tag is present
    sa_col_index = ""
    for i,col in enumerate(cols):
        if "SA:Z" in col:
            sa_col_index = i
        else:
            pass

    if sa_col_index != "":

        # SA-tag present
        sa_col = cols[sa_col_index].replace("SA:Z:","")
        split_alignments = sa_col.split(";")[:-1]

        # Iterate over the chimeric alignments
        # Based on insert_match look for genomic or insert mapped alignments
        chim_alignments = []
        if insert_match == True:
            # Look for insert matches
            for split_align in split_alignments:
                if split_align.split(",")[0] == insertID:
                    chim_alignments.append(split_align)
                else:
                    pass
        else:
            # Look for genomic matches
            for split_align in split_alignments:
                if split_align.split(",")[0] != insertID:
                    chim_alignments.append(split_align)
                else:
                    pass

        # Check how many acceptd chimeric matches are present
        if len(chim_alignments) == 0:
            return(None)
        elif len(chim_alignments) == 1:
            return(chim_alignments)
        else:
            return("multihits")
    else:
        # No SA-tag
        return(None)

# Find alignment length
def get_match_length(cigar):
    # Number of matching and deleted bases
    cigar_MD = re.findall("\d+(?=M)", cigar) + re.findall("\d+(?=D)", cigar)
    return(sum(map(int,cigar_MD)))

# Find the clip-sites used from inference of the insertion sites
def get_clipsites(start, end, cigar):
    regex = re.compile("\d+H|\d+S")
    matches = re.finditer(regex,cigar)
    if matches != None:
        starts = [m.start() for m in matches]
        if len(starts) == 1:
            mstart = starts[0]
            if mstart == 0:
                return([int(start)])
            else:
                return([int(end) - 1])
        else:
            return([int(start), int(end) - 1])
    else:
        return([])

def countMean(ls):
    return(np.sum(np.array(ls, dtype = float)) / len(ls))

# Summarize the clipsite information
def summarizeClipsites(sorted_alignment_records, ignore_duplicates):

    current_chrom, current_pos = "",""
    summary_table = []
    summarized_dc = {"chrom":"", "pos":"", "mapq_genome" : [], "mapq_insert": [], "span_genome": [], "span_insert": [], "insert_distance": []}
    for record in sorted_alignment_records:
            chrom = record[2]
            pos = record[8]
            flag = record[1]
            if ((ignore_duplicates == True) & (isDuplicate(flag) == True)):
                pass
            else:
                if ((chrom == current_chrom) & (pos == current_pos)):
                    summarized_dc["mapq_genome"].append(record[6])
                    summarized_dc["mapq_insert"].append(record[12])
                    summarized_dc["span_genome"].append(record[7])
                    summarized_dc["span_insert"].append(record[13])
                    summarized_dc["insert_distance"].append(record[14])

                else:
                    if current_chrom != "" :
                        # Update the final table
                        # Calculate the average mapping quality and
                        # span for genomic and insert sites
                        count = len(summarized_dc["mapq_genome"])
                        summary_table.append([summarized_dc["chrom"],
                                summarized_dc["pos"],
                                count,
                                countMean(summarized_dc["mapq_genome"]),
                                countMean(summarized_dc["mapq_insert"]),
                                countMean(summarized_dc["span_genome"]),
                                countMean(summarized_dc["span_insert"]),
                                countMean(summarized_dc["insert_distance"])])

                    else:
                        pass

                    # Initialize a new summarized dictionary
                    summarized_dc = {"chrom":"", "pos":"", "mapq_genome" : [], "mapq_insert": [], "span_genome": [], "span_insert": [], "insert_distance": []}
                    summarized_dc["chrom"] = chrom
                    summarized_dc["pos"] = pos
                    summarized_dc["mapq_genome"].append(record[6])
                    summarized_dc["mapq_insert"].append(record[12])
                    summarized_dc["span_genome"].append(record[7])
                    summarized_dc["span_insert"].append(record[13])
                    summarized_dc["insert_distance"].append(record[14])
                    current_chrom, current_pos = chrom, pos

    # The final update
    count = len(summarized_dc["mapq_genome"])
    summary_table.append([summarized_dc["chrom"],
                          summarized_dc["pos"],
                          count,
                          countMean(summarized_dc["mapq_genome"]),
                          countMean(summarized_dc["mapq_insert"]),
                          countMean(summarized_dc["span_genome"]),
                          countMean(summarized_dc["span_insert"]),
                          countMean(summarized_dc["insert_distance"])])

    return(summary_table)

# Function for creating a bed-file
# The bed-file includes all the primary and secondary alignments
# The primary alignments are indicated in red and the secondary in blue
def writeAlignmentBed(outdir, sorted_output_entries, insertID, ignore_duplicates):
    out = open(outdir + "/clipsite_reads.bed", "w")
    out.write("trackName=\"clipsite_reads\" itemRgb=\"On\"\n")
    for item in sorted_output_entries:
        name = item[0]
        primary_alignment = item[-1]
        flag = item[1]
        duplicate = isDuplicate(flag)
        if ((ignore_duplicates == False) | (duplicate == False)):
            if primary_alignment == "primary_insert":
                # Genomic part is considered secondary
                prim_chr = insertID
                prim_start = str(int(item[9]) - 1)
                prim_end = str(item[10])
                prim_score = item[12]
                prim_strand = "."
                prim_thick_start = prim_start
                prim_thick_end = prim_start
                prim_rgb = "255,0,0"

                out.write("\t".join([prim_chr, prim_start,
                                    prim_end, name, prim_score,
                                    prim_strand, prim_thick_start,
                                    prim_thick_end, prim_rgb]) + "\n")

                sec_chr = item[2]
                sec_start = str(item[3] - 1)
                sec_end = str(item[4])
                sec_score = item[6]
                sec_strand = "."
                sec_thick_start = sec_start
                sec_thick_end = sec_start
                sec_rgb = "0,0,255"

                out.write("\t".join([sec_chr, sec_start,
                                    sec_end, name, sec_score,
                                    sec_strand, sec_thick_start,
                                    sec_thick_end, sec_rgb]) + "\n")

            else:
                # Genomic part is considered primary
                prim_chr = item[2]
                prim_start = str(item[3] - 1)
                prim_end = str(item[4])
                prim_score = item[6]
                prim_strand = "."
                prim_thick_start = prim_start
                prim_thick_end = prim_start
                prim_rgb = "255,0,0"

                out.write("\t".join([prim_chr, prim_start,
                                    prim_end, name, prim_score,
                                    prim_strand, prim_thick_start,
                                    prim_thick_end, prim_rgb]) + "\n")

                sec_chr = insertID
                sec_start = str(int(item[9]) - 1)
                sec_end = str(item[10])
                sec_score = item[12]
                sec_strand = "."
                sec_thick_start = sec_start
                sec_thick_end = sec_start
                sec_rgb = "0,0,255"

                out.write("\t".join([sec_chr, sec_start,
                                    sec_end,name, sec_score,
                                    sec_strand, sec_thick_start,
                                    sec_thick_end, sec_rgb]) + "\n")
        else:
            pass

    out.close()

# Function sweeps through windows genomic and selects the
# position with highest clipsite count as the most likely insert site
def determineInsertSites(outdir, summary_table, min_readcount):

    #Initialize the first region
    prev_chrom = summary_table[0][0]
    prev_start = int(summary_table[0][1])
    region_entries = [summary_table[0]]
    insert_sites = []
    for item in summary_table[1:]:
        if ((item[0] == prev_chrom) & (int(item[1]) < prev_start + 100)):
            region_entries.append(item)
            prev_chrom = item[0]
            prev_start = item[1]
        else:
            # The position does not belong to neighbourhood of the previous
            # positions
            # Find the most likely insert site based on maximum clipsite count
            max_entry = region_entries[np.argmax([int(reg[2]) for reg in region_entries])]
            if int(max_entry[2]) >= min_readcount:
                insert_sites.append(max_entry)
            else:
                pass

            # Initialize the new region
            region_entries = []
            region_entries.append(item)

            # Initialize a new region
            prev_chrom = item[0]
            prev_start = item[1]


    max_entry = region_entries[np.argmax([int(reg[2]) for reg in region_entries])]
    if int(max_entry[2]) >= min_readcount:
        insert_sites.append(max_entry)
    else:
        pass

    output = open(outdir + "/insertionsite_candidates.txt","w")
    output.write("\t".join(["chrom", "pos", "count" , "average_mapq_genome", "average_mapq_insert", "average_span_genome", "average_span_insert", "average_insert_end_distance"]) + "\n")
    for item in insert_sites:
        output.write("\t".join(map(str,item)) + "\n")
    output.close()

# Main function
def main():

    # The contig name for insert
    insertID = args.insertID

    # Create outdir if does not exist
    outdir = args.outdir

    # The lenght of the insert sequence
    insertLength = int(args.insertLength)

    # The minimum required reads to consider a position to be
    # an insertion site candidate
    min_readcount = args.min_readcount

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Preprosess the bam file (Sort by read name)
    procfin = sortN(args.inputBAM, args.outdir)

    # Start processing the prerared input file
    output_entries = []

    with open(outdir + "/ntest_nohead.sam") as infile:
        while True:
            ln1 = infile.readline()
            ln2 = infile.readline()
            if ln2 != "":
                if ln1[0] == "@":
                     # Skip the header lines
                     pass
                else:

                        # Parse the paired entries together

                        # The first entry
                        cols1 = ln1[:-1].split("\t")
                        name1 = cols1[0]
                        flag1 = cols1[1]
                        chrom1 = cols1[2]
                        pair_chrom1 = cols1[6]

                        # The second entry
                        cols2 = ln2[:-1].split("\t")
                        name2 = cols2[0]
                        flag2 = cols2[1]
                        chrom2 = cols2[2]
                        pair_chrom2 = cols2[6]

                        # Check that the entries are truly paired reads
                        if name1 == name2:

                            if ((chrom1 != insertID) & (pair_chrom1 == insertID)):

                                # The first read of the pair maps to genome and the
                                # pair maps to insert

                                # The genomic alignment is primary and the insert
                                # alignment is secondary

                                start_genome = int(cols1[3])
                                mapq_genome = cols1[4]
                                cigar_genome = cols1[5]
                                alignment_length_genome = get_match_length(cigar_genome)
                                end_genome = start_genome + alignment_length_genome - 1
                                clipsites_genome = get_clipsites(start_genome, end_genome, cigar_genome)

                                # Fetch the rest of the data from the pair
                                start_insert = int(cols2[3])
                                mapq_insert = cols2[4]
                                cigar_insert = cols2[5]
                                alignment_length_insert = get_match_length(cigar_insert)
                                end_insert = start_insert + alignment_length_insert - 1

                                insert_end_distance = min([int(start_insert)-1, insertLength - int(end_insert)])

                                # Append entry
                                # Make each clipsite as separate entry
                                for clipsite in clipsites_genome:
                                    output_entries.append([name1, flag1, chrom1, int(start_genome),
                                                      int(end_genome),
                                                      cigar_genome, mapq_genome,
                                                      str(alignment_length_genome),
                                                      clipsite,
                                                      start_insert, end_insert,
                                                      cigar_insert, mapq_insert,
                                                      str(alignment_length_insert),insert_end_distance,"primary_genome"])

                            elif ((chrom1 != insertID) & (pair_chrom1 != insertID)):

                                # Neither reads of the pair are mapped to insert

                                # Try to find chimeric alignments
                                # First test the first entry
                                chim_alignments_1 = parseChimericAlignments(cols1, insertID, insert_match = True)
                                if ((chim_alignments_1 == None) | (chim_alignments_1 == "multihits")):
                                    # Test the second entry
                                    chim_alignments_2 = parseChimericAlignments(cols2, insertID, insert_match = True)
                                    if ((chim_alignments_2 == None) | (chim_alignments_2 == "multihits")):
                                        # No chimeric alignments or multihits
                                        pass
                                    else:
                                        start_genome = int(cols2[3])
                                        mapq_genome = cols2[4]
                                        cigar_genome = cols2[5]
                                        alignment_length_genome = get_match_length(cigar_genome)
                                        end_genome = start_genome + alignment_length_genome
                                        clipsites_genome = get_clipsites(start_genome, end_genome, cigar_genome)

                                        start_insert, strand_insert, cigar_insert, mapq_insert = chim_alignments_2[0].split(",")[1:-1]
                                        alignment_length_insert = get_match_length(cigar_insert)
                                        end_insert = int(start_insert) + alignment_length_insert - 1

                                        insert_end_distance = min([int(start_insert)-1, insertLength - int(end_insert)])

                                        # Append entry
                                        # Make each clipsite as separate entry
                                        for clipsite in clipsites_genome:
                                            output_entries.append([name2, flag2, chrom2, int(start_genome),
                                                                   int(end_genome),
                                                                   cigar_genome, mapq_genome,
                                                                   str(alignment_length_genome),
                                                                   clipsite,
                                                                   start_insert, end_insert,
                                                                   cigar_insert, mapq_insert,
                                                                   str(alignment_length_insert),insert_end_distance,"primary_genome"])

                                else:
                                    start_genome = int(cols1[3])
                                    mapq_genome = cols1[4]
                                    cigar_genome = cols1[5]
                                    alignment_length_genome = get_match_length(cigar_genome)
                                    end_genome = start_genome + alignment_length_genome
                                    clipsites_genome = get_clipsites(start_genome, end_genome, cigar_genome)

                                    start_insert, strand_insert, cigar_insert, mapq_insert = chim_alignments_1[0].split(",")[1:-1]
                                    alignment_length_insert = get_match_length(cigar_insert)
                                    end_insert = int(start_insert) + alignment_length_insert - 1

                                    # Append entry
                                    # Make each clipsite as separate entry
                                    insert_end_distance = min([int(start_insert)- 1, insertLength - int(end_insert)])

                                    for clipsite in clipsites_genome:
                                        output_entries.append([name1, flag1, chrom1, int(start_genome),
                                                           int(end_genome),
                                                           cigar_genome, mapq_genome,
                                                           str(alignment_length_genome),
                                                           clipsite,
                                                           start_insert, end_insert,
                                                           cigar_insert, mapq_insert,
                                                           str(alignment_length_insert),insert_end_distance,"primary_genome"])


                            elif ((chrom1 == insertID) & (pair_chrom1 == "=")):
                                # Both reads map to insert

                                # Try to find chimeric alignments
                                # First test the first entry
                                chim_alignments_1 = parseChimericAlignments(cols1, insertID, insert_match = False)
                                if ((chim_alignments_1 == None) | (chim_alignments_1 == "multihits")):
                                    # Test the second entry
                                    chim_alignments_2 = parseChimericAlignments(cols2, insertID, insert_match = False)
                                    if ((chim_alignments_2 == None) | (chim_alignments_2 == "multihits")):
                                        # No chimeric alignments or multihits
                                        pass
                                    else:
                                        start_insert = int(cols2[3])
                                        mapq_insert = cols2[4]
                                        cigar_insert = cols2[5]
                                        alignment_length_insert = get_match_length(cigar_insert)
                                        end_insert = start_insert + alignment_length_insert - 1

                                        # Extract information from the SA tag
                                        chrom, start_genome, strand_genome, cigar_genome, mapq_genome = chim_alignments_2[0].split(",")[:-1]
                                        alignment_length_genome = get_match_length(cigar_genome)
                                        end_genome = int(start_genome) + alignment_length_genome
                                        clipsites_genome = get_clipsites(int(start_genome), end_genome, cigar_genome)

                                        insert_end_distance = min([int(start_insert)-1, insertLength - int(end_insert)])

                                        # Append entry
                                        # Make each clipsite as separate entry
                                        for clipsite in clipsites_genome:
                                            output_entries.append([name1, flag1, chrom, int(start_genome),
                                                           int(end_genome),
                                                           cigar_genome, mapq_genome,
                                                           str(alignment_length_genome),
                                                           clipsite,
                                                           start_insert, end_insert,
                                                           cigar_insert, mapq_insert,
                                                           str(alignment_length_insert),insert_end_distance,"primary_insert"])

                                else:
                                    start_insert = int(cols1[3])
                                    mapq_insert = cols1[4]
                                    cigar_insert = cols1[5]
                                    alignment_length_insert = get_match_length(cigar_insert)
                                    end_insert = start_insert + alignment_length_insert - 1

                                    # Extract information from the SA tag
                                    chrom, start_genome, strand_genome, cigar_genome, mapq_genome = chim_alignments_1[0].split(",")[:-1]
                                    alignment_length_genome = get_match_length(cigar_genome)
                                    end_genome = int(start_genome) + alignment_length_genome
                                    clipsites_genome = get_clipsites(int(start_genome), end_genome, cigar_genome)

                                    insert_end_distance = min([int(start_insert)-1, insertLength - int(end_insert)])

                                    # Append entry
                                    # Make each clipsite as separate entry
                                    for clipsite in clipsites_genome:
                                        output_entries.append([name1, flag1, chrom, int(start_genome),
                                                           int(end_genome),
                                                           cigar_genome, mapq_genome,
                                                           str(alignment_length_genome),
                                                           clipsite,
                                                           start_insert, end_insert,
                                                           cigar_insert, mapq_insert,
                                                           str(alignment_length_insert),insert_end_distance,"primary_insert"])

                            else:

                                # The 1. pair maps to insert but the mate maps to genome
                                start_insert = int(cols1[3])
                                mapq_insert = cols1[4]
                                cigar_insert = cols1[5]
                                alignment_length_insert = get_match_length(cigar_insert)
                                end_insert = start_insert + alignment_length_insert - 1

                                start_genome = int(cols2[3])
                                mapq_genome = cols2[4]
                                cigar_genome = cols2[5]
                                alignment_length_genome = get_match_length(cigar_genome)
                                end_genome = start_genome + alignment_length_genome
                                clipsites_genome = get_clipsites(start_genome, end_genome, cigar_genome)

                                insert_end_distance = min([int(start_insert)-1, insertLength - int(end_insert)])

                                # Append entry
                                # Make each clipsite as separate entry
                                for clipsite in clipsites_genome:
                                    output_entries.append([name2, flag2, chrom2, int(start_genome),
                                                           int(end_genome),
                                                           cigar_genome, mapq_genome,
                                                           str(alignment_length_genome),
                                                           clipsite,
                                                           start_insert, end_insert,
                                                           cigar_insert, mapq_insert,
                                                           str(alignment_length_insert),insert_end_distance,"primary_insert"])


            else:
                break;

    print(output_entries[0])

    # Sort entries based on the clipsites
    sorted_output_entries = sorted(output_entries, key = itemgetter(2,8))

    # Write to bed
    writeAlignmentBed(outdir, sorted_output_entries, insertID, args.iduplicates)

    # Write sorted entries to an intermediate file
    output_intermediate = open(outdir + "/selected_insert_boundary_pairs.txt","w")
    output_intermediate.write("\t".join(["name","flag","chrom","startg","endg","cigarg","mapqg","alengthg","clipsite","starti","endi","cigari","mapqi","alengthi","insert_end_distance","primary_aligment"]) + "\n")
    for item in sorted_output_entries:
        output_intermediate.write("\t".join(map(str,item)) + "\n")
    output_intermediate.close()

    # Summarize the clipsite information (aggregate the clipsites in the same
    # positions)
    summary_table = summarizeClipsites(sorted_output_entries, args.iduplicates)

    # Write results to output
    output = open(outdir + "/summarized_clipsites.txt","w")
    output.write("\t".join(["chrom", "pos", "count" , "average_mapq_genome", "average_mapq_insert", "average_span_genome", "average_span_insert", "average_insert_end_distance"]) + "\n")
    for item in summary_table:
        output.write("\t".join(map(str,item)) + "\n")
    output.close()

    # Select most likely insert sites based on the summarized clipsite data
    insertsites_table = determineInsertSites(outdir, summary_table, min_readcount)

if __name__=="__main__":
    main()
