# multiple transcript
from unique_transcripts import unique_transcript
from statistics import mean
from merge_overlapping_intervals import mergeIntervals

def is_multiple_transcript(transcript_a, transcript_a_allfeatures_positions, transcript_a_aligned_b_genes):
    unique_count=0
    calls = []
    genes = []
    scores = []
    for gene in transcript_a_aligned_b_genes.keys():
        #print(transcript_a,  transcript_a_aligned_b_genes[gene]['perc_identity'])
        interval_list=[[a,b] for a, b in transcript_a_aligned_b_genes[gene][transcript_a]['aln_pos']]  # b+1 to check if adding one base overlaps with nearest interval
        score = mean(transcript_a_aligned_b_genes[gene]['perc_identity'])
        transcript_b_positions = mergeIntervals(sorted(interval_list) )

        call = unique_transcript(transcript_a_allfeatures_positions, transcript_b_positions)
        if call == 'exact_match' or call == 'all_jxn_match' or call == 'target_contained':
            unique_count +=1
            calls.append(call)
            genes.append(gene)
            scores.append(score)

    return unique_count, calls, genes, scores
