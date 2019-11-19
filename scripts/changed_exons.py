# changed exon - The source transcript contains exons which overlap target transcript exons but junctions differ

import compare_source_target_exons

def is_changed_exon_incl_kept_intron(source_positions, target_positions):
    transcript_b_aln_positions_paired = [(target_positions[i-1], target_positions[i]) for i in range(1, len(target_positions))]

    for exon1, exon2 in transcript_b_aln_positions_paired:
        for source_exon in source_positions:
            # in this, source exon is across intron between two target exon but source exon boundary not same as target exon boundary
            if compare_source_target_exons.is_source_exon_in_target_exon([exon1[1], exon2[0]], source_exon) and (exon1[0] != source_exon[0] and exon2[1] != source_exon[1]):
                return True


def is_changed_exon(source_positions, target_positions):
    # changed_exon : The source transcript contains exons which overlap target transcript exons with at least one junction match
    total_source_positions = len(source_positions)
    non_junction_match=0
    junction_match=0
    for target_exon in target_positions:
        for source_exon in source_positions:

            if source_exon[0] < target_exon[0] < source_exon[1] < target_exon[1]:
                #      ------------ target
                # -----------   source
                non_junction_match +=1
                source_positions.remove(source_exon)
                break
            elif target_exon[0] < source_exon[0] < target_exon[1] < source_exon[1]:
                # ------------ target
                #     ------------ source
                non_junction_match +=1
                source_positions.remove(source_exon)
                break
            elif source_exon[0] < target_exon[0] < source_exon[1] == target_exon[1]:
                #     --------------  target
                # ------------------  source
                junction_match +=1
                source_positions.remove(source_exon)
                break
            elif target_exon[0] == source_exon[0] < target_exon[1] > source_exon[1]:
                # ---------- target
                # -------------- source
                junction_match +=1
                source_positions.remove(source_exon)
                break
            elif source_exon[0] < target_exon[0] < target_exon[1] < source_exon[1]:
                #      ------------ target
                #  ------------------ source
                non_junction_match +=1
                source_positions.remove(source_exon)
                break
            else:
                pass
    #print("junciton match ", junction_match, "non junction match" , non_junction_match)
    if junction_match >= 1 and junction_match + non_junction_match == total_source_positions:
        return True

def is_new_exons(source_positions, target_positions):
    # changed_exon : The source transcript contains exons which overlap target transcript exons with at least one junction match
    total_source_positions = len(source_positions)
    non_junction_match=0
    junction_match=0
    for target_exon in target_positions:
        for source_exon in source_positions:

            if source_exon[0] == target_exon[0] or  source_exon[1] ==target_exon[1]:
                junction_match +=1
                source_positions.remove(source_exon)
                break
            else:
                pass
    #print("junciton match ", junction_match, "non junction match" , non_junction_match)
    if junction_match == 0:
        return True
