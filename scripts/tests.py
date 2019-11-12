import analysis
import compare_source_target_exons
from novel_retained_intron import is_novel_retained_intron
from unique_transcripts import unique_transcript
from changed_exons import is_changed_exon, is_changed_exon_incl_kept_intron
from gene_start_end_fragments import is_gene_start_overlap, is_gene_end_overlap
from target_intronic import is_target_intronic

if compare_source_target_exons.is_source_exon_left_equal_right_longer((5,12), (5,10)):
    print("True")
else:
    print("False")

if compare_source_target_exons.is_source_exon_right_equal_left_longer((10, 20), (8,20)):
    print("True")
else:
    print("False")

target_exon=[(1,20),(30,50), (60, 100),(150, 250)]
source_exon = [(4,20), (30,50), (60,100), (150,350)]

if unique_transcript(source_exon, target_exon) == 'all_jxn_match':
    print('all_jxn_match passed')
else:
    print('all_jxn_match failed')

target_exon=[(1,20),(30,50), (60, 100),(150, 250), (270, 340)]
source_exon = [(35,50), (60,100), (270, 340)]

if unique_transcript(source_exon, target_exon) == 'source_contained':
    print('source_contained passed')
else:
    print('source_contained failed')

source_exon=[(1,20),(25,50), (60, 100),(150, 250), (270, 345)]
target_exon = [(35,50), (60,100), (150, 250),(270, 340)]

if unique_transcript(source_exon, target_exon) == 'target_contained':
    print('target_contained passed')
else:
    print('target_contained failed')

target_exon=[(1, 25),(35, 50), (60,100), (150, 250), (270, 340)]
source_exon=[(20,50), (60, 100),(150, 250), (270, 345)]

if is_novel_retained_intron(source_exon, target_exon):
    print('novel_retained_intron passed')
else:
    print('novel_retained_intron failed')

target_exon=[(1, 25), (35, 50), (60,100), (150, 250), (270, 340)]
source_exon=[(1, 25), (40, 100),(150, 250), (270, 345)]

if is_novel_retained_intron(source_exon, target_exon):
    print('novel_retained_intron passed')
else:
    print('novel_retained_intron failed')

target_exon=[(1, 25), (35, 50), (60,100), (150, 250), (270, 340),(450, 500)]
source_exon=[(1, 25), (40, 100),(150, 250), (270, 500)]

if is_novel_retained_intron(source_exon, target_exon):
    print('novel_retained_intron passed')
else:
    print('novel_retained_intron failed')

target_exon=[(1, 25), (35, 50), (60,100), (150, 250), (270, 340)]
source_exon=[(1, 345)]

if is_novel_retained_intron(source_exon, target_exon):
    print('novel_retained_intron passed')
else:
    print('novel_retained_intron failed')

target_exon=[(1, 25),(35, 50), (60,100), (150, 250), (270, 340)]
source_exon=[(20,55), (60, 100)]

if is_changed_exon_incl_kept_intron(source_exon, target_exon):
    print('changed_exon_incl_kept_intron passed')
else:
    print('changed_exon_incl_kept_intron failed')

target_exon=[(1, 25),(35, 50), (60,100), (150, 250), (270, 340)]
source_exon=[(20,45)]

if is_changed_exon_incl_kept_intron(source_exon, target_exon):
    print('changed_exon_incl_kept_intron passed')
else:
    print('changed_exon_incl_kept_intron failed')

target_exon=[(1, 25),(35, 50), (60,100), (150, 250), (270, 340)]
source_exon=[(230, 320),(500, 600)]

if is_changed_exon_incl_kept_intron(source_exon, target_exon):
    print('changed_exon_incl_kept_intron passed')
else:
    print('changed_exon_incl_kept_intron failed')

target_exon=[(1, 25),(35, 50), (60,100), (150, 250), (270, 340)]
source_exon=[(1, 20),(30, 55), (57,103), (265, 345)]

if is_changed_exon(source_exon, target_exon):
    print('changed_exon  passed')
else:
    print('changed_exon failed')

target_exon=[(1, 25),(35, 50), (60,100), (150, 250), (270, 340)]
source_exon=[(30, 50), (57,103), (145, 254), (265, 345)]

if is_changed_exon(source_exon, target_exon):
    print('changed_exon  passed')
else:
    print('changed_exon failed')

target_exon=[(35, 50), (60,100), (150, 250), (270, 340)]
source_exon=[(10, 40)]

if is_gene_start_overlap(source_exon, target_exon):
    print('gene_start_overlap  passed')
else:
    print('gene_start_overlap failed')

target_exon=[(35, 50), (60,100), (150, 250), (270, 340)]
source_exon=[(290, 400)]

if is_gene_end_overlap(source_exon, target_exon):
    print('gene_end_overlap  passed')
else:
    print('gene_end_overlap failed')

target_exon=[(35, 50), (60,100), (150, 250), (320, 340)]
source_exon=[(270, 290), (295,310)]
if is_target_intronic(source_exon, target_exon):
    print('target_intronic passed')
else:
    print('target_intronic failed')
