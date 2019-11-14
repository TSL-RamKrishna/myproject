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
from target_intronic import is_target_intronic

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

def homology_analysis(blastAB, blastBA, gffA, gffB, species1name, species2name):
    print("Parsing A->B align blast data")
    speciesA_blast_speciesBgenes = parse_blast_result(blastAB)
    print("Parsing B->a align blast data")
    speciesB_blast_speciesAgenes = parse_blast_result(blastBA)

    #print(speciesA_blast_speciesBgenes)
    #print(speciesB_blast_speciesAgenes)
    A_gff_transcript_data, A_gff_gene_data = parse_gff_data(gffA)
    B_gff_transcript_data, B_gff_gene_data = parse_gff_data(gffB)

    # checking A transcripts homology to B transcript

    for transcript_a in A_gff_transcript_data.keys():

        if transcript_a in speciesA_blast_speciesBgenes.keys():
            #print(transcript_a)
            if transcript_a in A_gff_transcript_data.keys():
                pass
            else:
                continue
            transcript_aligned_called=False
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
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript]= {'call':'unique_transcript'}
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript].update({'category':unique_transcript_call})
                        #print(transcript_a, b_gene_transcript, transcript_a_gene, b_gene_aligned_to, " call ", ' unique_transcript ', unique_transcript_call )
                        transcript_aligned_called = True
                        continue

                    if is_novel_retained_intron(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1):
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript] = {'call':'absent_transcript'}
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript].update({'category':'novel_retained_intron'})
                        #print(transcript_a, b_gene_transcript, transcript_a_gene,b_gene_aligned_to, " call ", 'novel_retained_intron')
                        transcript_aligned_called = True
                        continue

                    if is_changed_exon_incl_kept_intron(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1):
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript]={'call':'absent_transcript'}
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript].update({'category':'changed_exon_incl_kept_intron'})
                        #print(transcript_a, b_gene_transcript, transcript_a_gene,b_gene_aligned_to, " call ", 'changed_exon_incl_kept_intron')
                        transcript_aligned_called = True
                        continue
                    if is_changed_exon(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1):
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript]={'call':'absent_transcript'}
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript].update({'category':'changed_exons'})
                        #print(transcript_a, b_gene_transcript, transcript_a_gene,b_gene_aligned_to, ' call ', ' changed_exons')
                        transcript_aligned_called = True
                        continue
                    if is_gene_start_overlap(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1):
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript]= {'call':'absent_transcript'}
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript].update({'category':'genic_start_overlap'})
                        #print(transcript_a, b_gene_transcript, transcript_a_gene,b_gene_aligned_to, ' call ', ' gene_start_overlap')
                        transcript_aligned_called = True
                        continue
                    if is_gene_end_overlap(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1):
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript]={'call':'absent_transcript'}
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript].update({'category':'genic_end_overlap'})
                        #print(transcript_a, b_gene_transcript, transcript_a_gene,b_gene_aligned_to, ' call ', ' gene_end_overlap')
                        transcript_aligned_called = True
                        continue
                    if is_antisense(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1, transcript_a_strand, transcript_b_strand):
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript]={'call':'absent_gene'}
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript].update({'category':'antisense'})
                        #print(transcript_a, b_gene_transcript, transcript_a_gene,b_gene_aligned_to, ' call ', 'antisense')
                        transcript_aligned_called = True
                        continue
                    if is_target_intronic(transcript_a_allfeatures_positions, b_gene_transcript_allfeatures_starting_1):
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript]={'call':'absent_gene'}
                        speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript].update({'category':'target_intronic'})
                        #print(transcript_a, b_gene_transcript, transcript_a_gene,b_gene_aligned_to, ' call ', 'target_intronic')
                        transcript_aligned_called = True
                        continue
                #if transcript_aligned_called == False:
                #print(transcript_a, " ", " ", b_gene_aligned_to, ' call ', ' not found')
        else:
            pass


    for transcript_a in A_gff_transcript_data.keys():
        transcript_a_gene = A_gff_transcript_data[transcript_a]['gene']
        if transcript_a in speciesA_blast_speciesBgenes.keys():
            #print(transcript_a)

            for b_gene_aligned_to in speciesA_blast_speciesBgenes[transcript_a].keys():
                b_gene_transcript_list = B_gff_gene_data[b_gene_aligned_to]['transcript']
                call_transcript_a_to_b={}
                for b_gene_transcript in b_gene_transcript_list:
                    if b_gene_transcript in speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to].keys():
                        call = speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript]['call']
                        category = speciesA_blast_speciesBgenes[transcript_a][b_gene_aligned_to][b_gene_transcript]['category']
                        if call in call_transcript_a_to_b.keys():
                            call_transcript_a_to_b[call].update({category:[b_gene_aligned_to, b_gene_transcript]})
                        else:
                            call_transcript_a_to_b={call:{category:[b_gene_aligned_to, b_gene_transcript]}}

                for key in call_transcript_a_to_b.keys():
                    if key == 'unique_transcript':
                        category_keys = call_transcript_a_to_b['unique_transcript'].keys()
                        if len(category_keys) ==1:
                            # this is unique transcript, exact match
                            for category_key in call_transcript_a_to_b['unique_transcript'].keys():
                                print(species1name + "\t" + species2name + "\t" + transcript_a + "\t" + call_transcript_a_to_b['unique_transcript'][category_key][1] + "\t" + 'unique_transcript' + "\t" + transcript_a_gene + "\t"  + call_transcript_a_to_b['unique_transcript'][category_key][0] + "\t" + category_key)
                        elif len(category_keys) > 1:
                            for category_key in call_transcript_a_to_b['unique_transcript'].keys():
                                print(transcript_a + "\t" + call_transcript_a_to_b['unique_transcript'][category_key][1] + "\t" + 'unique_transcript' + "\t" + transcript_a_gene + "\t"  + call_transcript_a_to_b['unique_transcript'][category_key][0] + "\t" + category_key)
                    else:
                        for category_key in call_transcript_a_to_b[key].keys():
                            print(species1name + "\t" + species2name + "\t" + transcript_a + "\t" + call_transcript_a_to_b[key][category_key][1] + "\t" + key + "\t" + transcript_a_gene + "\t"  + call_transcript_a_to_b[key][category_key][0] + "\t" + category_key)
        else:
            print(species1name + "\t" + species2name + "\t" + transcript_a + "" + "\t" + transcript_a_gene + "\t" + "\t" + "\t")



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
    parser.add_argument('-N', '--species1name', action='store', dest='species1name', default='A', help='Name of species1')
    parser.add_argument('-n', '--species2name', action='store', dest='species2name', default='B', help='Name of species2')


    options=parser.parse_args()
    #print(options)
    #homology_analysis("results/test_A_transcripts_align_Bgenes.txt","results/test_B_transcripts_align_Agenes.txt","data/testA.gff3","data/testB.gff3", options.species1name, options.species2name)
    homology_analysis("results/A_transcripts_align_Bgenes.txt", "results/B_transcripts_align_Agenes.txt", "data/A.gff3", "data/B.gff3", options.species1name, options.species2name)


    exit(0)
