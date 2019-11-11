

def is_exon_exact_match(source_exon, target_exon):
    # ---------
    # ---------
    a1,b1=target_exon
    a2,b2=source_exon
    if a1 == a2 and b1 == b2:
        return True

def is_source_exon_right_equal_left_longer(source_exon, target_exon):
    #  ----------------------
    #       -----------------

    a1,b1=target_exon
    a2,b2=source_exon
    if a1 < a2 < b1 == b2:
        return True

def is_source_exon_right_equal_left_shorter(source_exon, target_exon):
    #     -------------------
    #  ----------------------

    a1,b1=target_exon
    a2,b2=source_exon
    if a1 > a2 < b1 == b2:
        return True

def is_source_exon_right_longer_left_longer(source_exon, target_exon):
    #  ------------------
    #      -----------------
    a1,b1=target_exon
    a2,b2=source_exon
    if a1 < a2 < b1 < b2:
        return True

def is_source_exon_left_equal_right_shorter(source_exon, target_exon):
    # --------------------
    # ----------------
    a1,b1=target_exon
    a2,b2=source_exon
    if a1 == a2 and b1 > b2:
        return True

def is_source_exon_left_equal_right_longer(source_exon, target_exon):
    # --------------------
    # -----------------------
    a1,b1=target_exon
    a2,b2=source_exon
    if a1 == a2 and b1 < b2:
        return True

def is_source_exon_left_shorter_right_shorter(source_exon, target_exon):
    #     ------------------
    # ------------------
    a1,b1=target_exon
    a2,b2=source_exon

    if a2 < a1 < b2 < b1:
        return True

def is_source_exon_longer_both_sides(source_exon, target_exon):
    #     a1 --------- b1
    # a2------------------ b2
    a1,b1=target_exon
    a2,b2=source_exon
    if a2 < a1 < b1 < b2:
        return True

def is_target_exon_longer(source_exon, target_exon):
    #  a1------------------b1
    #     a2------------b2
    a1,b1=target_exon
    a2,b2=source_exon
    if a1 < a2 <b2 < b1:
        return True

def is_source_exon_in_target_exon(source, target):
    #  -----------   ------------
    #    -------     ------------
    if source[0] >= target[0] < source[1] <= target[1]:
        return True
    elif source[0] < target[0] < source[1] == target[1]:
        return True
    elif source[0] < target[0] < source[1] == target[1]:
        return True
    else:
        pass

def is_source_exon_one_end_equal(source_exon, target_exon):
    # --------------
    # ----------
    # or
    # ----------------
    #     ------------

    a1,b1=target_exon
    a2,b2=source_exon
    if (a1 == a2 and b2 <= b1 ) or (a2 >= a1 and b1==b2):
        return True

def right_overlap(source, target):
    if source[0] < target[0] and source[1] <= target[1]:
        return True
def left_overlap(source, target):
    if source[0] >=target[0] and source[0] < target[1] and source[1] >= target[1]:
        return True
def a_completely_overlaps_b(source, target):
    if source[0] <= target[0]  and source[1] >= target[1]:
        return True
