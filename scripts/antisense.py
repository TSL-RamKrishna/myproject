# antisense - Source transcript aligns to opposite strand of the only target gene in that location

def is_antisense(source_positions, target_positions, source_strand, target_strand):
    # check if source transcript is antisense to target transcript
    if source_strand != target_strand:  # meaning they are in different strands of dna
        total = 0
        for target_exon in target_positions:
            for source_exon in source_positions:
                if source_exon[0] < target_exon[0] < source_exon[1] <= target_exon[1]:
                    total +=1
                elif source_exon[0] <= target_exon[0] < source_exon[1]:
                    total +=1

        if total == len(source_positions):
            return True
