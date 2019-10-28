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
    return annotation_data.get_positions()

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
    A_gff_data = parse_gff_data(gffA)
    B_gff_data = parse_gff_data(gffB)

    print(A_gff_data)
    print(B_gff_data)

    # checking A transcripts homology to B transcript

    for transcript_a in speciesA_blast_speciesBgenes.keys():
        print(transcript_a)

        for b_gene_aligned_to in speciesA_blast_speciesBgenes[transcript_a].keys():

            print(transcript_a, " aligned to ", b_gene_aligned_to)


            b_gene_aln_start_pos = speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['start']
            b_gene_aln_end_pos = speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['end']

            b_gene_aln_pos = sorted(zip(b_gene_aln_start_pos, b_gene_aln_end_pos) )

            transcript_aln_pos = sorted( speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][transcript_a]['aln_pos'] )

            total_aligned_positions = len(b_gene_aln_start_pos )

            print("transcript ", transcript_a, transcript_aln_pos, " gene ", b_gene_aligned_to, b_gene_aln_pos)
            b_gene_transcript = B_gff_data[b_gene_aligned_to]['transcript']
            b_gene_transcript_exon_positions = B_gff_data[b_gene_aligned_to][b_gene_transcript]['exon']
            b_gene_transcript_allfeatures = B_gff_data[b_gene_aligned_to][b_gene_transcript]['allfeatures']

            b_gene_chr_start_pos, b_gene_chr_end_pos = B_gff_data[b_gene_aligned_to]['start'],B_gff_data[b_gene_aligned_to]['end']
            # positions are paired [start, end]
            # these positions are with respective to the chromosome position

            b_gene_transcript_exon_positions_starting_1 = sorted( set( get_exon_positions_starting_1(b_gene_transcript_exon_positions, b_gene_chr_start_pos) ) )
            b_gene_transcript_allfeatures_starting_1 = sorted( set( get_exon_positions_starting_1(b_gene_transcript_allfeatures, b_gene_chr_start_pos) ) )


            b_gene_chromosome = B_gff_data[b_gene_aligned_to]['chr']
            b_gene_strand = B_gff_data[b_gene_aligned_to]['strand']

            print("B gene chr ",b_gene_chromosome, " transcript ", b_gene_transcript, " allfeatures ", b_gene_transcript_allfeatures, b_gene_transcript_allfeatures_starting_1)
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
    homology_analysis("results/test_A_transcripts_align_Bgenes.txt","results/test_B_transcripts_align_Agenes.txt","data/testA.gff3","data/testB.gff3")


    exit(0)
