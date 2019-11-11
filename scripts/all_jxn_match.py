import sys

def is_all_jxn_match(source_positions, target_positions):
    # test for all_jxn_match
    # criteria:
    # number of exons in both source and target should be same
    # no intron matching, exon right and/or left end should be same

    if len(source_positions) == len(target_positions):
        counter=0
        for i in range(len(source_positions)):

            a1, b1 =  target_positions[i]
            a2, b2 = source_positions[i]

            if a2 <= a1 and b2>=b1:
                counter+=1
        if counter==len(source_positions):
            return True

if __name__ == '__main__':
    target_exon=[(4,20),(30,50), (60, 100),(150, 250)]
    source_exon = [(1,20), (30,50), (60,100), (150,350)]
    print(is_all_jxn_match(source_exon, target_exon) )
