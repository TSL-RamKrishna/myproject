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
from changed_exons import is_changed_exon, is_changed_exon_incl_kept_intron
from novel_retained_intron import is_novel_retained_intron
from unique_transcripts import unique_transcript
import compare_source_target_exons
from gene_start_end_fragments import is_gene_start_overlap, is_gene_end_overlap
from antisense import is_antisense

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

            b_gene_aln_start_pos = speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['start']
            b_gene_aln_end_pos = speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['end']
            transcript_a_strand = speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['qstrand']
            transcript_b_strand = speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['sstrand']

            b_gene_aln_pos = sorted(zip(b_gene_aln_start_pos, b_gene_aln_end_pos) )

            transcript_aln_pos = sorted( speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][transcript_a]['aln_pos'] )

            total_aligned_positions = len(b_gene_aln_start_pos )

            #print("transcript ", transcript_a, transcript_aln_pos, " gene ", b_gene_aligned_to, b_gene_aln_pos)
            b_gene_transcript_list = B_gff_gene_data[b_gene_aligned_to]['transcript']
            aligned_gene_found=False

            for b_gene_transcript in b_gene_transcript_list:
                b_gene_transcript_exon_positions = B_gff_gene_data[b_gene_aligned_to][b_gene_transcript]['exon']
                b_gene_transcript_allfeatures = B_gff_gene_data[b_gene_aligned_to][b_gene_transcript]['allfeatures']

                b_gene_chr_start_pos, b_gene_chr_end_pos = B_gff_gene_data[b_gene_aligned_to]['start'],B_gff_gene_data[b_gene_aligned_to]['end']
                # positions are paired [start, end]
                # these positions are with respective to the chromosome position

                b_gene_transcript_exon_positions_starting_1 = sorted( set( get_exon_positions_starting_1(b_gene_transcript_exon_positions, b_gene_chr_start_pos) ) )
                b_gene_transcript_allfeatures_starting_1 = sorted( set( get_exon_positions_starting_1(b_gene_transcript_allfeatures, b_gene_chr_start_pos) ) )


                b_gene_chromosome = B_gff_gene_data[b_gene_aligned_to]['chr']

                #print(b_gene_transcript, transcript_a_strand, transcript_b_strand)

                #print(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1)
                unique_transcript_call = unique_transcript(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1)
                if unique_transcript_call:
                    speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['call'].append('unique_transcript')
                    speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['category']=unique_transcript_call
                    print(transcript_a, " aligned to ", b_gene_aligned_to, " call ", ' unique_transcript ', unique_transcript_call )
                    aligned_gene_found = True
                    continue

                if is_novel_retained_intron(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1):
                    speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['call'].append('novel_retained_intron')
                    print(transcript_a, " aligned to ", b_gene_aligned_to, " call ", 'novel_retained_intron')
                    aligned_gene_found = True
                    continue

                if is_changed_exon_incl_kept_intron(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1):
                    speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['call'].append('changed_exon_incl_kept_intron')
                    print(transcript_a, " aligned to ", b_gene_aligned_to, " call ", 'changed_exon_incl_kept_intron')
                    aligned_gene_found = True
                    continue
                if is_changed_exon(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1):
                    speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['call'].append('changed_exons')
                    print(transcript_a, " aligned to ", b_gene_aligned_to, ' call ', ' changed_exons')
                    aligned_gene_found = True
                    continue
                if is_gene_start_overlap(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1):
                    speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['call'].append('gene_start_overlap')
                    print(transcript_a, " aligned to ", b_gene_aligned_to, ' call ', ' gene_start_overlap')
                    aligned_gene_found = True
                    continue
                if is_gene_end_overlap(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1):
                    speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['call'].append('gene_end_overlap')
                    print(transcript_a, " aligned to ", b_gene_aligned_to, ' call ', ' gene_end_overlap')
                    aligned_gene_found = True
                    continue
                if is_antisense(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1, transcript_a_strand, transcript_b_strand):
                    speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['call'].append('antisense')
                    print(transcript_a, " aligned to ", b_gene_aligned_to, ' call ', 'antisense')
                    aligned_gene_found = True
                    continue
            if aligned_gene_found == False:
                print(transcript_a, " aligned to ", b_gene_aligned_to, ' call ', ' not found')



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
