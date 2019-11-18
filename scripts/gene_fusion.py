# gene fusion - Matches more than one gene in one alignment

from merge_overlapping_intervals import mergeIntervals
from statistics import mean

from changed_exons import is_changed_exon, is_changed_exon_incl_kept_intron
from novel_retained_intron import is_novel_retained_intron
from unique_transcripts import unique_transcript
import compare_source_target_exons
from gene_start_end_fragments import is_gene_start_overlap, is_gene_end_overlap
from antisense import is_antisense
from target_intronic import is_target_intronic


def add_intervals(array_values):
    total=0
    for a, b in array_values:
        total=b-a+1

    return total



def multiple_gene_merge_intervals(transcript_a, transcript_a_allfeatures_positions, transcript_a_aligned_b_genes):

    #print(transcript_a, transcript_a_allfeatures_positions, transcript_a_aligned_b_genes)
    merged_interval = []
    perc_identity = []
    for gene in transcript_a_aligned_b_genes.keys():
        interval_list=[[a,b] for a, b in transcript_a_aligned_b_genes[gene][transcript_a]['aln_pos']]  # b+1 to check if adding one base overlaps with nearest interval
        #print('interval list ', interval_list)
        merged_interval = mergeIntervals(sorted(interval_list) + merged_interval)
        #print('after merge ', merged_interval)
        perc_identity=transcript_a_aligned_b_genes[gene]['perc_identity'] + perc_identity

    return merged_interval, str(mean(perc_identity)), ";".join(transcript_a_aligned_b_genes.keys())


def call_gene_fusion(source_positions, transcript_positions):

    unique_transcript_call = unique_transcript(source_positions, transcript_positions)
    if unique_transcript_call:
        return unique_transcript_call
    if is_novel_retained_intron(source_positions, transcript_positions):
        return 'novel_retained_intron'
    if is_changed_exon_incl_kept_intron(source_positions, transcript_positions):
        return 'changed_exon_incl_kept_intron'
    if is_changed_exon(source_positions, transcript_positions):
        return 'changed_exons'
    if is_gene_start_overlap(source_positions, transcript_positions):
        return 'genic_start_overlap'
    if is_gene_end_overlap(source_positions, transcript_positions):
        return 'genic_end_overlap'
    #if is_antisense(source_positions, transcript_positions, transcript_a_strand, transcript_b_strand):
    #    return 'antisense'
    if is_target_intronic(source_positions, transcript_positions):
        return 'target_intronic'
