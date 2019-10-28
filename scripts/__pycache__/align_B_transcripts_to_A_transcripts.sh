#!/bin/bash
source exonerate-2.2.0
exonerate --model affine:local --score 1 --query data/B.transcripts.fasta --target data/A.transcripts.fasta --showalignment no --showvulgar no --refine full --ryo "%qi %ql %qS %qt %qab %qae %qal %ti %tl %tS %tt %tab %tae %tal %tab %tae %tal %s %r %ps %g\n" > results/B_transcript_align_to_A_transcripts.txt 
