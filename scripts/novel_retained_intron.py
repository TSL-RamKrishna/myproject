# novel_retained_intron - Source transcript exons match compatible target junctions but include sequence from a target intron

import compare_source_target_exons
def is_novel_retained_intron(source_positions, target_positions):
    # test if transcript a has aligned to intronic region of transcript b,
    # getting pair of exonic positions from transcript b

    # test if source might align across all exons and introns of target
    if len(source_positions)==1 and compare_source_target_exons.is_source_exon_in_target_exon([target_positions[0][0], target_positions[-1][1] ], source_positions[0]):
        return True

    else:
        total = 0
        novel_retained_intron = 0
        total_source_positions=len(source_positions)
        for target_index_counter in range(len(target_positions)):
            if target_index_counter + 1 < len(target_positions):
                target_exon, target_exon1=target_positions[target_index_counter], target_positions[target_index_counter+1]
            else:
                target_exon = target_positions[target_index_counter]
            for source_exon in source_positions:
                # in this, source exon is across intron between target exons and source exon has boundary with target exon/s
                if target_index_counter + 1 < len(target_positions) and compare_source_target_exons.is_source_exon_in_target_exon([target_exon[1], target_exon1[0]], source_exon) and (target_exon1[1] == source_exon[1] or target_exon[0] == source_exon[0]):
                    #print('novel', source_exon, target_exon)
                    novel_retained_intron+=1
                    source_positions.remove(source_exon)
                    target_index_counter +=2
                    break
                else:
                    if compare_source_target_exons.is_source_exon_in_target_exon(source_exon, target_exon):
                        #print('equal', source_exon, target_exon)
                        total+=1
                        source_positions.remove(source_exon)

                    elif compare_source_target_exons.is_source_exon_left_equal_right_longer(source_exon, target_exon):
                        #print('left', source_exon, target_exon)
                        total+=1
                        source_positions.remove(source_exon)

                    elif compare_source_target_exons.is_source_exon_left_equal_right_shorter(source_exon, target_exon):
                        #print('right',source_exon, target_exon)
                        total+=1
                        source_positions.remove(source_exon)

        if novel_retained_intron > 0 and (total + novel_retained_intron * 2) == len(target_positions):
            return True
