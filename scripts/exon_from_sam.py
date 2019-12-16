import os, sys

def get_exon_positions(start, cigar):
    """
    from the start position, use cigar to extract exon positions
    """
    total=0 + int(start)
    exon_lengths = []
    number=''
    for char in cigar:
        if char in '0123456789':
            number+=char
        else:
            total+=int(number)
            number=''
        # elif char in 'MIDNSHP=X':
        #     if char == "M" or char == "I" or char == "D" or char == "X":
        #         match=int(number)
        #
        #     elif char=="N":
        #         nucleotide=int(number)
        #     elif char=="S":
        #         softclip=int(number)
        #     elif char=="H":
        #         hardclip=int(number)
        #     elif char=="P":
        #         padding=int(number)
        #     elif char=="X":
        #         mismatch=int(number)
    return (total )


def get_exon_positions_from_sam(samfile):
    """
    get exon positions from sam file
    """

    with open(samfile) as sam:
        for line in sam:
            if line.startswith("@"):    continue
            linearray=line.split()
            query=linearray[0]
            flag=linearray[1]
            target=linearray[2]
            align_start_position = linearray[3]
            quality=linearray[4]
            cigar=linearray[5]
            get_exon_positions(align_start_position, cigar)

print(get_exon_positions(sys.argv[1], sys.argv[2]) )
