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

def is_exon_exact_match(source_exon, target_exon):
    a1,b1=target_exon
    a2,b2=source_exon
    if a1 == a2 and b1 == b2:
        return

def is_exon_right_justified(source_exon, target_exon):
    a1,b1=target_exon
    a2,b2=source_exon
    if b1 == b2 and a1 < a2:
        return True

def is_exon_left_justified(source_exon, target_exon):
    a1,b1=target_exon
    a2,b2=source_exon
    if a1 == a2 and b1 > b2:
        return True

def is_source_exon_longer(source_exon, target_exon):
    a1,b1=target_exon
    a2,b2=source_exon
    if b2-a2 > b1-a1:
        return True

def is_target_exon_longer(source_exon, target_exon):
    a1,b1=target_exon
    a2,b2=source_exon
    if b1-a1 > b2-a2:
        return True

def is_source_exon_overlap_longer_than_target_exon(source_exon, target_exon):
    a1,b1=target_exon
    a2,b2=source_exon
    if(a2 < a1 and b2 > b1):
        return True

def is_novel_retained_intron(target_exon1, target_exon2, source_exon):
    # this is true if there is intron in between two target exons or source mapped across all target exons
    a1,b1=target_exon1
    a2,b2=target_exon2
    a3,b3=source_exon

    # test if there is intron in between target exon1 and target exon2
    if (b3 == b2 and a3 < b1) or (a3==a1 and b3 < b2):
        return True
    # if (a1==a3 or b2==b3) and (b3-a3 > b1-a1 or b3-a3 > b2-a2):
    #     return True
    # test if source exon mapped across all first and last exons
    elif a3 < a1 and b3 > b2:
        return True

def is_changed_exon_incl_kept_intron(target_exon1, target_exon2, source_exon):
    # test if there is intron between target exons and source exon overlaps more than target exons
    a1,b1 = target_exon1
    a2,b2 = target_exon2
    a3,b3 = source_exon

    if (a3 < a1 or b3 > b2):
        return True
    elif (a3 < b1 and b3 > a2):
        return True


def is_changed_exon(target_exon, source_exon):
    # test if all source exons overlap target exons
    a1,b1 = target_exon
    a2,b2 = source_exon

    pass

def is_center_left(a1,b1,a2,b2):
    if a1>a2:
        return True

    return None

def is_a_in_b(source, target):
    # a and b are range
    # return is range a is in range b
    if source[0] >= target[0] and source[1] <= target[1]:
        return True

def right_overlap(source, target):
    if source[0] < target[0] and source[1] <= target[1]:
        return True
def left_overlap(source, target):
    if source[0] >=target[0] and source[0] < target[1] and source[1] >= target[1]:
        return True
def a_completely_overlaps_b(source, target):
    if source[0] <= target[0]  and source[1] >= target[1]:
        return True
    else:
        return False

def exact_match(source_positions, target_positions):
    total = 0
    for i in range(len(source_positions)):
        if is_exon_exact_match(source_positions[i], target_positions[i]):
            total +=1
    if total == len(source_positions):
        return True

def all_jxn_match(source_positions, target_positions):

    if len(source_positions) == len(target_positions):
        counter=0
        for i in range(len(source_positions)):

            a1, b1 =  source_positions[i]
            a2, b2 = target_positions[i]

            if a1 <= b2 and b1 >=a2:
                counter+=1
        if counter==len(source_positions):
            return True

def source_contained(source_positions, target_positions):
    # if not all_jxn_match, lets look other
    total_contained=0
    if len(source_positions) < len(target_positions):
        for i in range(len(source_positions)):
            if is_exon_exact_match(target_positions[i], source_positions[i]) or is_exon_left_justified(target_positions[i], source_positions[i]) or is_exon_right_justified(target_positions[i], source_positions[i]):
                total_contained=0

        if total_contained == len(source_positions):
            return True


def target_contained(source_positions, target_positions):
    total_contained=0

    for target_exon in target_positions:
        for source_exon in source_positions:
            if is_a_in_b(target_exon, source_exon):
                total_contained+=1
            else:
                pass
    if len(source_positions) >= total_contained and total_contained >=1:
        return True


def unique_transcript(source_positions, target_positions):
    if len(source_positions) ==1 and  len(target_positions) ==1  and source_positions == target_positions:
        return 'mono_exonic'
    elif exact_match(source_positions, target_positions):
        return 'exact_match'
    elif all_jxn_match(source_positions, target_positions):
        return 'all_jxn_match'
    elif source_contained(source_positions, target_positions):
        return 'source_contained'
    elif target_contained(source_positions, target_positions):
        return 'target_contained'


def novel_retained_intron(source_positions, target_positions):
    # test if transcript a has aligned to intronic region of transcript b,
    # getting pair of exonic positions from transcript b

    transcript_b_aln_positions_paired = [(target_positions[i-1], target_positions[i]) for i in range(1, len(target_positions))]

    for exon1, exon2 in transcript_b_aln_positions_paired:
        for source_exon in source_positions:
            if is_novel_retained_intron(exon1, exon2, source_exon):
                return True

    # test if source might align across all exons and introns of target
    if len(source_positions)>=2 and is_novel_retained_intron(target_positions[0], target_positions[-1], source_positions):
        return True

def changed_exon(source_positions, target_positions):
    number_of_source_exons = len(source_positions)

    for source_exon in source_positions:

    return

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
                unique_transcript_call = unique_transcript(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1)
                if unique_transcript_call:
                    speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['call']='unique_transcript'
                    speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['category']=unique_transcript_call

                    continue

                if novel_retained_intron(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1):
                    speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['call']='novel_retained_intron'
                    continue

                if is_changed_exon_incl_kept_intron(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1):
                    speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['call']='is_changed_exon_incl_kept_intron'
                    continue
                if changed_exon(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1):
                    speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to]['call']='changed_exon'
                    continue



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
