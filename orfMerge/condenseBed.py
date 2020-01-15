# Python 3.6
'''
# Author: Travis Law
# Date: 08/20/2018
# Objective: Condense Bed Files
'''
FASTA = '/ahg/regevdata/projects/Ribo-seq/MHC-I/ref/orfeome/B721.RibORF.Price/B721.RibORF.Price.canonical.prot.fasta'  # noqa
BED = '/ahg/regevdata/projects/Ribo-seq/MHC-I/ref/orfeome/B721.RibORF.Price/B721.RibORF.Price.canonical.bed'  # noqa
OUT = '/ahg/regevdata/projects/Ribo-seq/MHC-I/ref/orfeome/B721.RibORF.Price/B721.RibORF.Price.dedup.bed'  # noqa
fasta = open(FASTA, 'r')
bed = open(BED, 'r')
out = open(OUT, 'w')
keep = set()
for line in fasta:
    if line[0] == '>':
        ID = line[1:line.find('::')]
        keep.add(ID)
fasta.close()
duplicates = set()
for line in bed:
    name = line.split('\t')[3]
    if name in keep and name not in duplicates:
        out.write(line)
        duplicates.add(name)
bed.close()
out.close()
