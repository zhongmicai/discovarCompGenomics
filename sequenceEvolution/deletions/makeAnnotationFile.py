import csv
import re

def read_csv(file, header=None):
    '''
    Takes a csv file, 
    outputs a list, where each element is a list which contains the cells of a single row as strings.
    '''
    data = []    
    reader = csv.reader(open(file, 'rU'))
    for row in reader:
        data.append(row)
    if header == True:
        return data[1:]
    else:
        return data

def makeAnnotationFile(geneTable, outFile):
    table=read_csv(geneTable,header=True)
    output=open(outFile,"a+")
    for line in table:
        ID=line[0]
        GOs=line[6]+line[7]+line[8]
        terms = re.findall("GO:[0-9]+", GOs)
        for term in terms:
            output.write('''%s = %s\n''' % (ID,term[3:]))
    output.close()
        