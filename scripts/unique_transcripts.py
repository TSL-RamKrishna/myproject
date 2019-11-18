
import compare_source_target_exons


def exact_match(source_positions, target_positions):

    if len(source_positions) != len(target_positions):
        return

    total_source_exons = len(source_positions)
    total = 0

    # check first source exon
    source_exon = source_positions[0]
    target_exon = target_positions[0]
    if source_exon[0] <= target_exon[0] < source_exon[1] == target_exon[1]:
        total+=1

    #check last source exon
    source_exon = source_positions[-1]
    target_exon = target_positions[-1]
    if source_exon[0] == target_exon[0] < target_exon[1] <= source_exon[1]:
        total+=1

    for target_exon in target_positions[1:-2]:
        for source_exon in source_positions[1:-2]:

            if source_exon[0] == target_exon[0] and source_exon[1] == target_exon[1]:
                total+=1
            else:
                pass
    if total == total_source_exons:
        return True



def unique_transcript_call(source_positions, target_positions):
    # if not all_jxn_match, lets look other
    total = 0
    total_source_exons = len(source_positions)
    first_source_exon_found=False

    for source_exon_index in range(len(source_positions)):
        source_exon = source_positions[source_exon_index]
        for target_exon in target_positions:
            if first_source_exon_found==False and (source_exon[1] == target_exon[1]):
                total +=1
                first_source_exon_found=True
                break
            elif source_exon_index > 0 and source_exon_index < len(source_positions)-1 and (source_exon[0] == target_exon[0] and source_exon[1] == target_exon[1]):
                total +=1
                break
            elif source_exon_index == len(source_positions) -1 and (source_exon[0] == target_exon[0]):
                total +=1
                break

    return total



def unique_transcript(source_positions, target_positions):
    if len(source_positions) ==1 and  len(target_positions) ==1  and source_positions == target_positions:
        return 'mono_exonic'
    elif exact_match(source_positions, target_positions):
        return 'exact_match'
    elif unique_transcript_call(source_positions, target_positions) == len(source_positions) == len(target_positions):
        return 'all_jxn_match'
    elif unique_transcript_call(source_positions, target_positions) == len(source_positions) < len(target_positions):
        return 'source_contained'
    elif unique_transcript_call(source_positions, target_positions) == len(target_positions) < len(source_positions):
        return 'target_contained'
