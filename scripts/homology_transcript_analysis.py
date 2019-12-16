import os, sys
import gzip
import parse_gff3

sample_A_transcripts=parse_gff3.gff3(sys.argv[1])
sample_B_transcripts=parse_gff3.gff3(sys.argv[2])
alignment_A_to_B=sys.argv[3]




# get sample A transcripts gene Ids

sample_A_transcripts.get_transcript_gene()
sample_A_transcripts_gene=sample_A_transcripts.transcript_data


sample_B_transcripts.get_transcript_gene()
sample_B_transcripts_gene=sample_B_transcripts.transcript_data
#print(sample_A_transcripts_gene)
#print(sample_B_transcripts_gene)
#print("gff3 parse done")


def transcripts_A_to_B_unique_multiple_analysis(samfile):
    '''
    return unique or multiple transcript relation from transcripts in A to B
    '''

    transcript_alignments=dict()
    openfunction=gzip.open if samfile.endswith(".gz") else open
    with openfunction(samfile) as samInput:
            for record in samInput:
                    if record.rstrip().startswith("@"):     continue
                    if record.rstrip() == "":       continue
                    recordArray=record.rstrip().split("\t")
                    try:
                            Qtranscript=recordArray[0]
                            Flag=recordArray[1]
                            Rtranscript=recordArray[2]
                            Position=recordArray[3]
                            Quality=recordArray[4]
                            Cigar=recordArray[5]
                    except IndexError as ierr:
                            print("Record: ", recordArray[:6])
                            print("See the error record above")
                            exit(0)

                    if Qtranscript in transcript_alignments.keys():
                            transcript_alignments[Qtranscript].update({Rtranscript:{"F":Flag, "T":Rtranscript, "P":Position, "Q":Quality, "C":Cigar}})
                    else:
                            transcript_alignments[Qtranscript]={Rtranscript:{"F":Flag, "T":Rtranscript, "P":Position, "Q":Quality, "C":Cigar}}

    #print(transcript_alignments)

    for transcriptA in transcript_alignments.keys():
            list_of_transcriptB=list(transcript_alignments[transcriptA].keys() )
            flags=[ transcript_alignments[transcriptA][transcriptB]['F'] for transcriptB in list_of_transcriptB]
            #print(transcriptA, list_of_transcriptB, flags)


            if '0' in flags:
                if '256' in flags:
                    # first print the transcrptA perfectly aligned transcriptB
                    print("A\tB\t" + transcriptA + "\t" + transcriptB + "\tmultiple_transcript\t0\t" + sample_A_transcripts_gene[transcriptA] + "\t" + sample_B_transcripts_gene[list_of_transcriptB[0]]+ "\t" + "" )
                    for index in range(len(flags)):
                        if flags[index] == '256':
                            print("A\tB\t" + transcriptA + "\t" + list_of_transcriptB[index] + "\tmultiple_transcript\t0\t" + sample_A_transcripts_gene[transcriptA] + "\t" + sample_B_transcripts_gene[list_of_transcriptB[index]]+ "\t" + "" )
                        else:
                            print("Other flag found ", flags[index], list_of_transcriptB[index])
                else:
                    transcriptB=list_of_transcriptB[0]
                    if transcriptB in sample_A_transcripts_gene.keys():
                        print("A\tB\t" + transcriptA + "\t" + list_of_transcriptB[0] + "\tunique_transcript\t0\t" + sample_A_transcripts_gene[transcriptA] + "\t" + sample_B_transcripts_gene[transcriptB]+ "\t" + "" )
            else:
                pass

transcripts_A_to_B_unique_multiple_analysis(alignment_A_to_B)
