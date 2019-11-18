# get call category

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
    if is_antisense(source_positions, transcript_positions, transcript_a_strand, transcript_b_strand):
        return 'antisense'
    if is_target_intronic(source_positions, transcript_positions):
        return 'target_intronic'
