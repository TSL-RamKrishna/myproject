#!/usr/bin/env python
import sys
import gzip
import urllib.parse
import json

'''
A simple gff3 parser to return a dictionary object with start and end
position for all types e.g. gene, transcript, exon, intron, UTRs

Version  0.0.1
'''

class gff3():
    '''
    class to parse GFF3

    '''

    def __init__(self, gff3):
        '''
        gff3 - gff3 filename
        feature e.g. gene/mRNA/exon/intron/UTRs/
        '''
        self.in_file_gff3 = gff3
        self.gene_data=dict()
        self.transcript_data=dict()

    def parseGFFAttributes(self, attributeString):
        """Parse the GFF3 attribute column and return a dict"""#
        if attributeString == ".": return {}
        ret = {}
        for attribute in attributeString.split(";"):
            key, value = attribute.split("=")
            ret[urllib.parse.unquote(key)] = urllib.parse.unquote(value)
        return ret

    def coordinates(self, line):
        '''
        returns feature and its start, end coordinates
        '''
        line=line.rstrip()
        linearray=line.split()
        chromosome = linearray[0]
        feature = linearray[2]
        attributes = self.parseGFFAttributes(linearray[8])
        feature_id = attributes["ID"]
        parent_id = attributes["Parent"] if "Parent" in attributes.keys() else ""
        start = int(linearray[3])
        end = int(linearray[4])
        strand = linearray[6]
        #print("returning ", feature, feature_id, start, end, strand)
        return chromosome, feature, feature_id, parent_id, start, end, strand

    def get_positions(self):
        '''
        get start and end positions for all attributes
        '''
        #print("Getting positions")
        openFunc = gzip.open if self.in_file_gff3.endswith('.gz') else open

        infile = openFunc(self.in_file_gff3)
        for line in infile:

            if line.startswith("#") or line.rstrip() == "":
                continue

            chromosome, feature, feature_id, parent_id, start, end, strand = self.coordinates(line)
            #print(chromosome, feature, feature_id, parent_id, start, end, strand)
            if 'gene' in feature.lower():
                ## no two gene have same id, so just add to dict
                self.gene_data[feature_id] = {'start':start, 'end':end, 'strand':strand, 'chr':chromosome}
                last_geneid = feature_id
            elif 'rna' in feature.lower():
                self.transcript_data[feature_id]={'start':start, 'end': end, 'gene':parent_id   , 'allfeatures': [], 'strand':strand}
                if 'transcript' in self.gene_data[last_geneid].keys():
                    self.gene_data[last_geneid]['transcript'].append(feature_id)
                else:
                    self.gene_data[last_geneid]['transcript'] = [feature_id]

                self.gene_data[last_geneid].update({feature_id:{'start':start, 'end':end, 'allfeatures':[]}})
                last_transcript = feature_id
            elif 'exon' in feature.lower():
                # add all other features start and end to one array
                self.transcript_data[last_transcript]['allfeatures'].append((start, end))
                self.gene_data[last_geneid][last_transcript]['allfeatures'].append((start, end))

                if feature in self.gene_data[last_geneid][last_transcript].keys():
                    self.gene_data[last_geneid][last_transcript][feature].append((start,end))
                else:
                    self.gene_data[last_geneid][last_transcript].update({feature:[(start,end)]})

    def get_transcript_gene(self):
        '''
        get start and end positions for all attributes
        '''
        #print("Getting positions")
        openFunc = gzip.open if self.in_file_gff3.endswith('.gz') else open

        infile = openFunc(self.in_file_gff3)
        for line in infile:

            if line.startswith("#") or line.rstrip() == "":
                continue

            chromosome, feature, feature_id, parent_id, start, end, strand = self.coordinates(line)
            #print(chromosome, feature, feature_id, parent_id, start, end, strand)
            if 'rna' in feature.lower():
                self.transcript_data[feature_id]=parent_id

    def save_coordinates(self, outfilename):
        '''
        save the dict object containing attribute positions
        in json format in the filename provided
        '''
        print("Saving dict object in json")
        if len(self.transcript_data) == 0:
            self.get_positions()
        with open(outfilename, 'w') as outfile:
            json.dump(self.transcript_data, outfile)
    def print_coordinates(self):
        '''
        prints the dict object containing attribute positions
        in json format on screen
        '''
        if len(self.transcript_data) == 0:
            self.get_positions()
        print(json.dumps(self.transcript_data, indent=4, sort_keys=True))

    def get_transcript_data(self):
        return self.transcript_data

    def get_gene_data(self):
        return self.gene_data


if __name__=='__main__':

    sample=gff3(sys.argv[1])
    #sample.save_coordinates('output.json')
    sample.get_positions()
    print(sample.get_transcript_data() )

    print(sample.get_gene_data())
