for x in A B C; do cmd="transeq -sequence data/${x}.transcripts.fasta -outseq data/${x}.transcripts_pep.fasta -frame 1 -table 0";
	echo $cmd
	$cmd	
 done
