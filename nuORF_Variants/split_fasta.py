'''
# Author: Tamara Ouspenskaia
# Date: 01/14/2020
# Objective: This script generates all possible peptides given a fasta file
'''


import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--fasta_file", required = True, help = "Path to the FASTA file of proteins")
parser.add_argument("--shortest_peptide", required = True,
                    help = "Shortest length of peptides to be predicted")
parser.add_argument("--longest_peptide", required = True,
                    help = "Longest length of peptides to be predicted")

args = parser.parse_args()

def find_peptides(seq,shortest,longest):
    peptides = []
    longestplus1 = int(longest) + 1
    for i in range(int(shortest),longestplus1):
        unit = '.{' + str(i) + ',' + str(i) + '}'
        for j in range(0,i):
            out = re.findall(unit, seq[j:])
            peptides = peptides + out
    return(peptides)

for record in SeqIO.parse(args.fasta_file, 'fasta'):
    peptides=find_peptides(str(record.seq), args.shortest_peptide, args.longest_peptide)
    print("\n".join(i for i in peptides))
