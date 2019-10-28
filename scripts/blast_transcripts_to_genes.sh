sbatch -o ../log/A_transcripts_align_Bgenes.log -J A_B --mem 10G --cpus 4 --wrap  "blastn -query A.transcripts.fasta -db B_gene_sequences.fasta -outfmt 6 -perc_identity 98.0  -out ../results/A_transcripts_align_Bgenes.txt "
sbatch -o ../log/A_transcripts_align_Cgenes.log -J A_C --mem 10G --cpus 4 --wrap  "blastn -query A.transcripts.fasta -db C_gene_sequences.fasta -outfmt 6 -perc_identity 90.0  -out ../results/A_transcripts_align_Cgenes.txt "

