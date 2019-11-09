import analysis
from analysis import all_jxn_match, source_contained, target_contained

if analysis.is_source_exon_left_equal_right_longer((5,12), (5,10)):
    print("True")
else:
    print("False")

if analysis.is_source_exon_right_equal_left_longer((10, 20), (8,20)):
    print("True")
else:
    print("False")

target_exon=[(1,20),(30,50), (60, 100),(150, 250)]
source_exon = [(4,20), (30,50), (60,100), (150,350)]

if all_jxn_match(source_exon, target_exon) == True:
    print('all_jxn_match passed')
else:
    print('all_jxn_match failed')

target_exon=[(1,20),(30,50), (60, 100),(150, 250), (270, 340)]
source_exon = [(35,50), (60,100), (270, 340)]

if source_contained(source_exon, target_exon) == True:
    print('source_contained passed')
else:
    print('source_contained failed')

source_exon=[(1,20),(25,50), (60, 100),(150, 250), (270, 345)]
target_exon = [(35,50), (60,100), (150, 250),(270, 340)]

if target_contained(source_exon, target_exon) == True:
    print('target_contained passed')
else:
    print('target_contained failed')

target_exon=[(1, 25),(35, 50), (60,100), (150, 250), (270, 340)]
source_exon=[(20,50), (60, 100),(150, 250), (270, 345)]

if analysis.novel_retained_intron(source_exon, target_exon):
    print('novel_retained_intron passed')
else:
    print('novel_retained_intron failed')

target_exon=[(1, 25), (35, 50), (60,100), (150, 250), (270, 340)]
source_exon=[(1, 25), (40, 100),(150, 250), (270, 345)]

if analysis.novel_retained_intron(source_exon, target_exon):
    print('novel_retained_intron passed')
else:
    print('novel_retained_intron failed')

target_exon=[(1, 25), (35, 50), (60,100), (150, 250), (270, 340),(450, 500)]
source_exon=[(1, 25), (40, 100),(150, 250), (270, 500)]

if analysis.novel_retained_intron(source_exon, target_exon):
    print('novel_retained_intron passed')
else:
    print('novel_retained_intron failed')

target_exon=[(1, 25), (35, 50), (60,100), (150, 250), (270, 340)]
source_exon=[(1, 345)]

if analysis.novel_retained_intron(source_exon, target_exon):
    print('novel_retained_intron passed')
else:
    print('novel_retained_intron failed')
