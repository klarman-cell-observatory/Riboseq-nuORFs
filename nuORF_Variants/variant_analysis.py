'''
# Author: Tamara Ouspenskaia
# Date: 01/14/2020
# Objective: This script takes the output from variant_to_protein script, filters it and outputs a new fasta of filtered proteins, variant counts per type, and an info file.
'''


import argparse
import pandas as pd
import numpy as np
from functools import reduce
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


parser = argparse.ArgumentParser()

parser.add_argument("--wt_orfs", required = True, help = "Path to fasta of wild-type proteins")
parser.add_argument("--mut_orfs", required = True, help = "Path to fasta of mutant proteins")
parser.add_argument("--mut_annot", required = True, help = "Path to the annotation file produced by variant_to_protein.sh script")
parser.add_argument("--filtered_fasta", required = True, help = "Path to output filtered fasta with variants")
parser.add_argument("--all_variants_counts", required = True, help = "Path to output file with variant counts per ORF type")
parser.add_argument("--high_qual_variants_counts", required = True, help = "Path to output file with high quality variant counts per ORF type")
parser.add_argument("--indel_analysis", required = True, help = "Path to output info file")

args = parser.parse_args()



#### Load wild-type and mutant protein fastas to be compared and the annotation output from variant_to_protein script:

def ids_seqs(fasta_file):
    ids = []
    seqs = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        ids.append(str(record.id))
        seqs.append(str(record.seq))
    df = pd.DataFrame({'header': ids,'seq': seqs})
    return(df)

def parse_mut_id(ORF_ID):
    ORF_ID = ORF_ID.split('|')[0:3]
    ORF_ID = '|'.join(ORF_ID)
    return ORF_ID

wt = ids_seqs(args.wt_orfs)
cols = ['ORF_ID','wt_seq']
wt.columns = cols

mut = ids_seqs(args.mut_orfs)
mut['ORF_ID'] = mut['header'].apply(parse_mut_id)
cols = ['header', 'mut_seq', 'ORF_ID']
mut.columns = cols

mut_vs_wt = mut.merge(wt, on='ORF_ID', how='left')

different = mut_vs_wt[mut_vs_wt['mut_seq'] != mut_vs_wt['wt_seq']]

expanded = different['ORF_ID'].str.split('|', expand=True)
cols = ['gene_name','ORF_ID','ORF_type']
expanded.columns = cols

different['ORF_type'] = expanded['ORF_type']

annot = pd.read_csv(args.mut_annot, header=0, sep='\t')
annot_to_keep = annot[['mut_location', 'variant_qc',
       'ORF_genomic_start', 'ORF_genomic_end', 'ORF_ID', 'trans_strand',
       'variant_id', 'trans_id','variant_ORF_id']]


#### Compare lengths and start codons of wild-type vs. mutant proteins:

def seq_len(row):
    row['wt_len'] = len(row['wt_seq'])
    row['mut_len'] = len(row['mut_seq'])
    conditions = [row['mut_len'] == row['wt_len'], row['mut_len'] > row['wt_len'], row['mut_len'] < row['wt_len']]
    choices = ['equal','mut_longer','mut_shorter']
    row['comp_lengths'] = np.select(conditions, choices, default=np.nan)
    return(row)

different = different.apply(seq_len, axis=1)

def start_comp(row):
    try:
        wt_start = row['wt_seq'][0]
        mut_start = row['mut_seq'][0]
        conditions = [
            wt_start == mut_start, wt_start != mut_start
        ]
        choices = ['same','different']
        row['comp_starts'] = np.select(conditions, choices, default=np.nan)
    except:
        print('This ORF had its start codon mutated to stop codon:')
        print(row['ORF_ID'])
    return(row)

different = different.apply(start_comp, axis=1)

different = different.merge(annot_to_keep, left_on='header', right_on='variant_ORF_id', how='left')


# #### Find variants unique to canonical_extended vs. canonical:

# Sort ORFs by ENST, stop, variant, wild-type length (strand-aware). Remove duplicates by ENST. Keep canonical.

types = ['canonical','canonical_extended','CDS']
can_ext = different[different['ORF_type'].isin(types)]

def split(df):
    plus_strand = df[df['trans_strand']=='+']
    minus_strand = df[df['trans_strand'] == '-']
    return plus_strand, minus_strand

def filter_plus_str(df):
    df = df.sort_values(['trans_id', 'ORF_genomic_end', 'mut_location', 'wt_len'], ascending=[True, True, True, True])
    df['group'] = df["trans_id"].map(str) + '_' + df["ORF_genomic_end"].map(str) + '_' + df["mut_location"].map(str)
    unique_plus_strand = df.drop_duplicates(subset=['group'], keep='first')
    unique_plus_strand.drop(['group'], axis =1, inplace=True)
    return unique_plus_strand

def filter_minus_str(df):
    df = df.sort_values(['trans_id', 'ORF_genomic_start', 'mut_location', 'wt_len'], ascending=[True, True, True, True])
    df['group'] = df["trans_id"].map(str) + '_' + df["ORF_genomic_start"].map(str) + '_' + df["mut_location"].map(str)
    unique_minus_strand = df.drop_duplicates(subset=['group'], keep='first')
    unique_minus_strand.drop(['group'], axis = 1, inplace=True)
    return unique_minus_strand

def filter_df(df):
    plus_strand, minus_strand = split(df)
    unique_plus = filter_plus_str(plus_strand)
    unique_minus = filter_minus_str(minus_strand)
    frames = [unique_plus, unique_minus]
    unique = pd.concat(frames)
    unique = unique.sort_index()
    return unique

filtered = filter_df(can_ext)

no_can = different[~different['ORF_type'].isin(types)]
frames = [no_can, filtered]
new_diff = pd.concat(frames)
new_diff = new_diff.sort_index()


# #### Remove duplicates between canonical_truncated and canonical, in case some canonical_truncated do not have a matching partner.
types = ['canonical','canonical_truncated','CDS','Trunc']
can_trunc = new_diff[new_diff['ORF_type'].isin(types)]

def split(df):
    plus_strand = df[df['trans_strand']=='+']
    minus_strand = df[df['trans_strand'] == '-']
    return plus_strand, minus_strand

def filter_plus_str(df):
    df = df.sort_values(['trans_id', 'ORF_genomic_end', 'mut_location', 'wt_len'], ascending=[True, True, True, False])
    df['group'] = df["trans_id"].map(str) + '_' + df["ORF_genomic_end"].map(str) + '_' + df["mut_location"].map(str)
    unique_plus_strand = df.drop_duplicates(subset=['group'], keep='first')
    unique_plus_strand.drop(['group'], axis =1, inplace=True)
    return unique_plus_strand

def filter_minus_str(df):
    df = df.sort_values(['trans_id', 'ORF_genomic_start', 'mut_location', 'wt_len'], ascending=[True, True, True, False])
    df['group'] = df["trans_id"].map(str) + '_' + df["ORF_genomic_start"].map(str) + '_' + df["mut_location"].map(str)
    unique_minus_strand = df.drop_duplicates(subset=['group'], keep='first')
    unique_minus_strand.drop(['group'], axis = 1, inplace=True)
    return unique_minus_strand

def filter_df(df):
    plus_strand, minus_strand = split(df)
    unique_plus = filter_plus_str(plus_strand)
    unique_minus = filter_minus_str(minus_strand)
    frames = [unique_plus, unique_minus]
    unique = pd.concat(frames)
    unique = unique.sort_index()
    return unique

filtered = filter_df(can_trunc)

no_can = new_diff[~new_diff['ORF_type'].isin(types)]

frames = [no_can, filtered]
new_diff = pd.concat(frames)
new_diff = new_diff.sort_index()


# #### Prepare a new filtered fasta for mass spec search:

def subset_fasta(df, input_fasta, output_fasta):
    record_dict = SeqIO.index(input_fasta, 'fasta')
    with open(output_fasta, "w") as output:
        for ORF_id in df['header']:
            record = SeqRecord(seq = record_dict[ORF_id].seq, id = ORF_id, description = '')
            SeqIO.write(record, output, "fasta")

subset_fasta(new_diff, args.mut_orfs, args.filtered_fasta)


# #### Save output files:

# All indels:

mut_counts = new_diff.groupby('ORF_type').count()
mut_counts['ORF_ID_x'].to_csv(args.all_variants_counts, sep='\t')

# Indels that pass QC metrics:

high_qual = new_diff[new_diff['variant_qc'] == 'PASS']
mut_counts_high_qual = high_qual.groupby('ORF_type').count()['ORF_ID_x']
mut_counts_high_qual.to_csv(args.high_qual_variants_counts, sep='\t')

# Variants annotation table:

new_diff.to_csv(args.indel_analysis, header=True, index=None, sep='\t')


#### Print out information to stdout:

print('Number of proteins with indels: ' + str(len(new_diff)))
print('Number of proteins with indels passing QC: ' + str(len(high_qual)))
