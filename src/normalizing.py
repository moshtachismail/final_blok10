"""
Author: Moshtach Ismail
This script converts the raw counts to normalized counts and removes
the versions out of the gene name
"""
import csv
import pandas as pd
import sys

def main(input, output):
    filter(input, output)
    print("normalizing the raw_counts is done. ")

def filter(input, output):
    """
    make a new file where only the first and last gene names are in
    the file, because the last column of the file already has the
    values that we need.
    :return: the new file with 1st and last column
    """
    print("normalizing.py started running")
    print(pd.__version__)
    print(csv.__version__)
    gene_values= [] # gene names
    values = []
    #df = pd.read_table(input, sep='\t', usecols=[0,3], skiprows=4)
    df = pd.read_table(input, sep='\t', skiprows=4)
    q = df.iloc[:,0] # only gets first value of each row
    r = df.iloc[:,1] # gets the second value

    for a in r:
        values.append(a)
    for x in q:
        gene_values.append(x)

    print("amount of values found", len(values))
    print("b", len(gene_values))

    #normalize values
    amin, amax = min(values), max(values)
    for i, val in enumerate(values):
        values[i] = (val - amin) / (amax - amin)
    #print(values) # prints normalized values

    # put the new values in a txt file
    with open(output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(gene_values, values))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
    #1. is the input file
    #2. is the output file

