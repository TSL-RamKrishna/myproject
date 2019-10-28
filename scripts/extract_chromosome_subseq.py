import os, sys
import get_coordinates
import json

class chromosome_sequence():
    '''
    gets transcript DNA seqeunce from chromosome/assembly file
    '''
    def __init__(self, reference, gff):
        # this takes reference as chromosome level assembly, transcript file and a gff file

        self.reference = reference
        self.gff =gff

    def open_gff(self):
        self.opengff = open(self.gff, 'r')
    def next(self):
        return self.opengff.readline()

    def get_gene_coordinates_from_gff(self):
        # reads gff file to get gene coordinates
        gene_data = dict()

        with open(self.gff) as filegff:
            for line in filegff:
                if line.startswith("#"):
                    continue
                # print(line)
                chromosome, feature, feature_id, start, end, strand = get_coordinates.coordinates(line)
                if chromosome in gene_data.keys():
                    if feature_id in gene_data[chromosome].keys():
                        gene_data[chromosome][feature_id].update()
                    gene_data[chromosome].update()
                if 'gene' in feature:
                    gene_data[feature_id] = {'feature': feature, 'start': start, 'end': end, 'strand': strand}

                # if feature in gene_data.keys():
                #     gene_data[feature_id].update()
                # else:
                #     gene_data[feature_id] = {'feature': feature, 'start': start, 'end': end, 'strand': strand}

        return gene_data

if __name__=='__main__':

    a  = chromosome_sequence(sys.argv[1], sys.argv[2])
    print(json.dumps(a.get_gene_coordinates_from_gff(), sort_keys = False, indent=4 ) )
