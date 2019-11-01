#!/usr/bin/env python

# this program takes genome assembly of two species
# their two transcripts in FASTA file
# their two GFF3 files
# if two species are A and B, the gene coordinates of A and B are obtained from
# respective gff3 files.
# Using the gene start and end coordinates, the sequences are extracted from
# the genome assembly and creates gene sequence FASTA files for A and B

# A Transcripts are blastn aligned to B genes and B Transcripts are blastn are
# aligned to A genes, generating two blastn outputs

# Blastn outputs and gff3 files are parsed to python dictionary objects

import parse_blastn_output
from parse_gff3 import gff3
from merge_overlapping_intervals import mergeIntervals

def parse_blast_result(blastfile):

    return parse_blastn_output.get_blast_data(blastfile)

def parse_gff_data(gff_file):
    annotation_data = gff3(gff_file)
    annotation_data.get_positions()
    return annotation_data.get_transcript_data(), annotation_data.get_gene_data()

def total_length_of_allfeatures(list_of_positions_start_end):
    # this return the total length of exons/cds
    total=0
    for a, b in list_of_positions_start_end:
        total += abs(b-a + 1)  # get absolute value, works if strand is + or -

    return total

def is_center_left(a1,b1,a2,b2):
    if a1>a2:
        return True

    return None

def is_a_in_b(a, b):
    # a and b are range
    # return is range a is in range b
    if a[0] >= b[0] and a[1] <=b[1]:
        return True

def right_overlap(a,b):
    if a[0] < b[0] and a[1] <=b[1]:
        return True
def left_overlap(a,b):
    if a[0] >=b[0] and a[0] < b[1] and a[1] > b[1]:
        return True
def a_completely_overlaps_b(a,b):
    if a[0] <= b[0]  and a[1] >= b[1]:
        return True
    else:
        return False

def all_jxn_match(transcript_a_aln_positions, gene_b_aln_positions):
    if len(transcript_a_aln_positions) == len(gene_b_aln_positions):
        a1, b1 =  transcript_a_aln_positions[0]
        a2, b2 = gene_b_aln_positions[0]

        if a1 <= b2 and b1 >=a2:
            return True

def source_contained(transcript_a_aln_positions, gene_b_aln_positions):
    # if not all_jxn_match, lets look other
    total_contained=0
    for a1, b1 in transcript_a_aln_positions:
        for a2, b2 in gene_b_aln_positions:
            if a1 > b2 :
                continue
            elif is_a_in_b([a1,b1], [a2,b2]):
                total_contained+=1
            else:
                pass
    if total_contained == len(transcript_a_aln_positions):
        return True

def target_contained(transcript_a_aln_positions, gene_b_aln_positions):
    total_contained=0
    for a1, b1 in transcript_a_aln_positions:
        for a2, b2 in gene_b_aln_positions:
            if is_a_in_b([a1, b1], [a2,b2]) or right_overlap([a1, b1], [a2,b2]) or left_overlap([a1, b1], [a2,b2]):
                total_contained+=1
            else:
                pass
    if len(transcript_a_aln_positions) >= total_contained and total_contained >=1:
        return True


def is_unique_transcript(transcript_a_aln_positions, gene_b_aln_positions):
    if transcript_a_aln_positions == gene_b_aln_positions:
        return 'exact_match'
    elif all_jxn_match(transcript_a_aln_positions, gene_b_aln_positions):
        return 'all_jxn_match'

    elif source_contained(transcript_a_aln_positions, gene_b_aln_positions):
        return 'source_contained'
    elif target_contained(transcript_a_aln_positions, gene_b_aln_positions):
        return 'target_contained'



def is_multiple_transcript(transcript_a, transcript_a_all_aligned_genes):
    pass


def get_transcript_details_from_gene(gene):
    pass

def get_exon_positions_starting_1(chr_exon_positions, chr_transcript_start_pos):
    # return the exon positions starting from 1, not in chromosome
    positions = []
    for start, end in chr_exon_positions:
        positions.append((start - chr_transcript_start_pos+1, end-chr_transcript_start_pos+1))

    return positions

def homology_analysis(blastAB, blastBA, gffA, gffB):
    print("Parsing A->B align blast data")
    speciesA_blast_speciesBgenes = parse_blast_result(blastAB)
    print("Parsing B->a align blast data")
    speciesB_blast_speciesAgenes = parse_blast_result(blastBA)

    #print(speciesA_blast_speciesBgenes)
    #print(speciesB_blast_speciesAgenes)
    A_gff_transcript_data, A_gff_gene_data = parse_gff_data(gffA)
    B_gff_transcript_data, B_gff_gene_data = parse_gff_data(gffB)

    # checking A transcripts homology to B transcript

    for transcript_a in speciesA_blast_speciesBgenes.keys():
        #print(transcript_a)
        if transcript_a in A_gff_transcript_data.keys():
            pass
        else:
            continue

        transcript_a_gene = A_gff_transcript_data[transcript_a]['gene']
        transcript_a_gene_start = A_gff_gene_data[transcript_a_gene]['start']
        transcript_a_gene_end = A_gff_gene_data[transcript_a_gene]['end']
        transcript_a_allfeatures_positions = sorted(set( get_exon_positions_starting_1(A_gff_transcript_data[transcript_a]['allfeatures'], transcript_a_gene_start ) ))

        # list of B genes that transcript_a aligned to
        list_of_b_genes = speciesA_blast_speciesBgenes[transcript_a].keys()
        for b_gene_aligned_to in speciesA_blast_speciesBgenes[transcript_a].keys():

            print(transcript_a, " aligned to ", b_gene_aligned_to)


            b_gene_aln_start_pos = speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['start']
            b_gene_aln_end_pos = speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['end']

            b_gene_aln_pos = sorted(zip(b_gene_aln_start_pos, b_gene_aln_end_pos) )

            transcript_aln_pos = sorted( speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][transcript_a]['aln_pos'] )

            total_aligned_positions = len(b_gene_aln_start_pos )

            #print("transcript ", transcript_a, transcript_aln_pos, " gene ", b_gene_aligned_to, b_gene_aln_pos)
            b_gene_transcript_list = B_gff_gene_data[b_gene_aligned_to]['transcript']
            for b_gene_transcript in b_gene_transcript_list:
                b_gene_transcript_exon_positions = B_gff_gene_data[b_gene_aligned_to][b_gene_transcript]['exon']
                b_gene_transcript_allfeatures = B_gff_gene_data[b_gene_aligned_to][b_gene_transcript]['allfeatures']

                b_gene_chr_start_pos, b_gene_chr_end_pos = B_gff_gene_data[b_gene_aligned_to]['start'],B_gff_gene_data[b_gene_aligned_to]['end']
                # positions are paired [start, end]
                # these positions are with respective to the chromosome position

                b_gene_transcript_exon_positions_starting_1 = sorted( set( get_exon_positions_starting_1(b_gene_transcript_exon_positions, b_gene_chr_start_pos) ) )
                b_gene_transcript_allfeatures_starting_1 = sorted( set( get_exon_positions_starting_1(b_gene_transcript_allfeatures, b_gene_chr_start_pos) ) )


                b_gene_chromosome = B_gff_gene_data[b_gene_aligned_to]['chr']
                b_gene_strand = B_gff_gene_data[b_gene_aligned_to]['strand']

                print(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1)
                unique_transcript_call = is_unique_transcript(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1)
                if unique_transcript_call:
                    print(transcript_a, b_gene_transcript, 'unique_transcript', transcript_a_gene, b_gene_aligned_to, unique_transcript_call)
                # if b_gene_aln_pos == b_gene_transcript_allfeatures_starting_1:
                #     print("B gene chr ",b_gene_chromosome, " transcript ", b_gene_transcript, " allfeatures ", b_gene_transcript_allfeatures, b_gene_transcript_allfeatures_starting_1)
                # for align_position in range(total_aligned_positions):
                #     a_start = speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['start'][align_position]
                #     a_end = speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['end'][align_position]
                #     b_start = speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][transcript_a]['start'][align_position]
                #     b_end = speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][transcript_a]['end'][align_position]
                #
                #     print("Transcript " , transcript_a, " aligned ", a_start, a_end, b_start, b_end )



if __name__=='__main__':

    import sys
    import argparse
    Description = 'Python program to find homology between transcripts from two closely related species'
    parser=argparse.ArgumentParser(description=Description)

    parser.add_argument('-A', '--species1assembly',  action='store', dest='species1assembly', help='Species 1 Assembly FASTA file')
    parser.add_argument('-a', '--species2assembly', action='store', dest='species2assembly', help='Species 2 Assembly FASTA file')
    parser.add_argument('-T','--species1transcript', action='store', dest='species1transcript', help='Speices 1 Transcript FASTA file')
    parser.add_argument('-t', '--species2transcript', action='store', dest='species2transcript', help='Species 2 Transcript FASTA file')
    parser.add_argument('-G','--species1gff',  action='store', dest='species1gff', help='Species 1 GFF3 file')
    parser.add_argument('-g','--species2gff',  action='store', dest='species2gff', help='Species 2 GFF3 file')


    options=parser.parse_args()
    #print(options)
    #homology_analysis("results/test_A_transcripts_align_Bgenes.txt","results/test_B_transcripts_align_Agenes.txt","data/testA.gff3","data/testB.gff3")
    homology_analysis("results/A_transcripts_align_Bgenes.txt", "results/B_transcripts_align_Agenes.txt", "data/A.gff3", "data/B.gff3")


    exit(0)
