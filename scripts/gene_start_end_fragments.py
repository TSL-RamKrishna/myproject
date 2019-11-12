# gene_start_overlap or gene_end_overlap - The source transcript partially overlaps the beginning or end of a target transcript

def is_gene_start_overlap(source_positions, target_positions):

    #            ---------     --------      -------- target
    # -------   ------- source
    if len(source_positions) < 1:
        return
    source_last_exon = source_positions[-1]
    target_first_exon = target_positions[0]
    if source_last_exon[0] < target_first_exon[0] < source_last_exon[1] < target_first_exon[1]:
        return True


def is_gene_end_overlap(source_positions, target_positions):
    # ----------   -----------  ----------  target
    #                               --------- ----- source
    if len(source_positions) < 1:
        return
    target_last_exon = target_positions[-1]
    source_first_exon = source_positions[0]

    if target_last_exon[0] <= source_first_exon[0] < target_last_exon[1] <= source_first_exon[1]:
        return True
