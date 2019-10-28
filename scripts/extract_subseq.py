import sys
import csv
from Bio import SeqIO


def get_coordinates(txtfile):

    data = {}
    with open(txtfile, 'r') as genepos:
        for line in genepos:
            line=line.rstrip()
            if line.startswith("#") or line=="":
                continue
            else:
                linearray=line.split()
                for position in range(len(linearray)):
                    chromosome=linearray[0]
                    gene = linearray[1]
                    start = int(linearray[2])
                    end = int(linearray[3])
                    if chromosome in data.keys():
                        data[chromosome].update({gene:[start, end]})
                    else:
                        data[chromosome] = {gene:[start, end]}
    return data
def get_subseq(gene_positions, fasta):
    with open(fasta) as fastahandle:
        for record in SeqIO.parse(fastahandle, 'fasta'):
            recordid=record.id
            ntseq = str(record.seq)
            if recordid in gene_positions.keys():
                for gene in sorted(gene_positions[recordid].keys()):

                    print(">" + gene + " chromosome " + recordid + " start " + str(gene_positions[recordid][gene][0]) + " end " + str(gene_positions[recordid][gene][1]))
                    print(ntseq[int(gene_positions[recordid][gene][0] ) -1 : int(gene_positions[recordid][gene][1])] )

gene_positions = get_coordinates(sys.argv[1])
get_subseq(gene_positions, sys.argv[2])
