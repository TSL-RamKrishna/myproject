#!/usr/bin/env python

import sys

def get_blast_data(blast_outfile):
    '''
    get blast result like which transcript aligned to which gene and their
    start and end positions
    e.g. transcriptA1={'geneB1':{'start':[25,175], 'end':[125, 225], 'transcriptA1':{'start':[1,150], 'end':[100,200]}, 'qstrand'='+', 'sstrand'='-', 'perc_identity':[100.0, 99.9]}}
    In the example above, transcriptA1 is aligned to geneB1.
    transcriptA1 1 to 100 is aligned to geneB1 25 to 125 and percent identity is 100.0
    transcriptA1 150 to 200 is aligned to geneB1 175 to 225 and percent identity is 99.9
    transcriptA1 is query, and its in positive '+' strand
    geneB1 is subject  and its in negative '-' strand

    '''

    transcript_gene_alignment = dict()
    with open(blast_outfile) as blastfile:
        for line in blastfile:
            line=line.rstrip()
            if line == "":
                continue

            linearray = line.split()
            query = linearray[0]
            subject = linearray[1]
            perc_identity = float(linearray[2])
            query_start = int(linearray[6])
            query_end = int(linearray[7])
            subject_start = int(linearray[8])
            subject_end = int(linearray[9])
            query_strand = '+'
            subject_strand = '+'

            # swap the position if the alignment is on the negative strand
            if subject_start > subject_end:
                subject_start, subject_end = subject_end, subject_start
                subject_strand = '-'
            if query_start > query_end:
                query_start, query_end = query_end, query_start
                query_strand = '-'

            if query in transcript_gene_alignment.keys():
                if subject in transcript_gene_alignment[query].keys():
                    transcript_gene_alignment[query][subject]['start'].append(subject_start)
                    transcript_gene_alignment[query][subject]['end'].append(subject_end)
                    transcript_gene_alignment[query][subject][query]['start'].append(query_start)
                    transcript_gene_alignment[query][subject][query]['end'].append(query_end)
                    transcript_gene_alignment[query][subject][query]['aln_pos'].append((query_start, query_end))
                    transcript_gene_alignment[query][subject]['perc_identity'].append(perc_identity)
                else:
                    transcript_gene_alignment[query].update({subject:{'start':[subject_start], 'end':[subject_end], query:{'aln_pos':[(query_start, query_end)], 'start':[query_start], 'end':[query_end]}, 'qstrand':query_strand, 'sstrand':subject_strand, 'perc_identity':[perc_identity]}})
            else:
                transcript_gene_alignment[query] = {subject:{'start':[subject_start], 'end':[subject_end], query:{'aln_pos':[(query_start, query_end)], 'start':[query_start], 'end':[query_end]}, 'qstrand':query_strand, 'sstrand':subject_strand, 'perc_identity':[perc_identity]}}


    return transcript_gene_alignment


if __name__=='__main__':

    result = get_blast_data(sys.argv[1])
    for transcript in result.keys():
        print(transcript, result[transcript])
