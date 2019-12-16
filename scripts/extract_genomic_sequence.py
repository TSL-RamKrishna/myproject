# extract gene sequence from chromosomes
import sys
import get_coordinates
import gzip
from Bio import SeqIO

class extract_genomic_sequence():
    # extract genomic sequence from start to end positions

    def __init__(self, fasta, gff, seqfeature):
        # get fasta file
        self.fasta = fasta
        self.gff = gff
        self.seqfeature=seqfeature
    def read_gene_coordinates_from_gff(self):
        # get chromosome for define sequence feature
        # default sequence feature is gene
        # user can pass sequence feature as mRNA or exon or intron to get the
        # coordinates for those sequence features only
        gene_data = {}

        with open(self.gff) as gfffile:
            for line in gfffile:
                line=line.rstrip()
                if line.startswith('#') or line=="":
                    continue
                else:

                    chromosome, feature, feature_id, start, end, strand = get_coordinates.coordinates(line)
                    if self.seqfeature == feature:
                        if chromosome in gene_data.keys():
                            gene_data[chromosome].update({feature_id:[start, end]} )
                        else:
                            gene_data[chromosome] = {feature_id:[start, end]}

        self.gene_data = gene_data

    def extract_gene_sequence(self):
        # extract from start to end
        self.read_gene_coordinates_from_gff()

        openfunc = gzip.open if self.fasta.endswith(".gz") else open
        fastaseq = openfunc(self.fasta)
        for record in SeqIO.parse(fastaseq, 'fasta'):
            recordid = record.id
            recordseq = str(record.seq)
            
            if recordid in self.gene_data.keys():

                for rec in self.gene_data[recordid].keys():
                    # print(recordid, rec, self.gene_data[recordid][rec][0], self.gene_data[recordid][rec][1]  )
                    start = int(self.gene_data[recordid][rec][0])
                    end = int(self.gene_data[recordid][rec][1])
                    print(">" + rec)
                    print(recordseq[start-1:end])   # as python string index starts at 0 and nt seq start from 1
                    exit(0)


if __name__ == '__main__':
    '''
    pass assembly fasta and gff files as positional arguments
    e.g. extract_genomic_sequence.py A.chromosome.fasta A.gff gene/mRNA/exon
    '''
    a = extract_genomic_sequence(sys.argv[1], sys.argv[2], sys.argv[3])
    a.extract_gene_sequence()
