# get chromosome gene coordinates

def coordinates(line):
    # returns feature and its start, end coordinates
    line=line.rstrip()
    linearray=line.split()
    chromosome = linearray[0]
    feature = linearray[2]
    feature_id = linearray[8].split(";")[0].replace("ID=", "")
    start = int(linearray[3])
    end = int(linearray[4])
    strand = linearray[6]
    #print("returning ", feature, feature_id, start, end, strand)
    return chromosome, feature, feature_id, start, end, strand
