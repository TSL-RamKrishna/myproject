# target intronic - A transcript is fully contained in other genome gene but overlaps no transcript exons of that gene

def is_target_intronic(source_positions, target_positions):
    # source is wholly in the intron of the target
    if len(source_positions) < 1:
        return
    elif len(source_positions) == 1:
        first_source_exon = source_positions[0]
        last_source_exon = source_positions[0]
    else:

        first_source_exon = source_positions[0]
        last_source_exon = source_positions[-1]

    for target_exon_index in range(len(target_positions)-1):

        target_exon_left = target_positions[target_exon_index]
        target_exon_right = target_positions[target_exon_index + 1]
        if target_exon_left[1] < first_source_exon[0] and last_source_exon[1] < target_exon_right[0]:
            return True
