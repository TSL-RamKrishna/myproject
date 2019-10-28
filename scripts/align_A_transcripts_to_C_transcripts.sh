#!/bin/bash

source exonerate-2.2.0
exonerate --model affine:local --score 1 --query data/A.transcripts.fasta --target data/C.transcripts.fasta --showalignment no --showvulgar no --refine full --ryo "%qi %ql %qS %qt %qab %qae %qal %ti %tl %tS %tt %tab %tae %tal %tab %tae %tal %s %r %ps %g\n" >  results/A_transcript_align_to_C_transcripts.txt
