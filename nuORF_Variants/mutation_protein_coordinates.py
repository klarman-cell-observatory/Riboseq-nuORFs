'''
# Author: Tamara Ouspenskaia
# Date: 01/14/2020
# Objective: This script generates the protein coordinates of the variant
'''

import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("--snv_annot", required = True, help = "Path to the annotation file produced by variant_to_protein.sh script")
parser.add_argument("--snv_analysis", required = True, help = "Path to output info file")
parser.add_argument("--output", required = True, help = "Path to output file")

args = parser.parse_args()

annot = pd.read_csv(args.snv_annot, sep='\t', header=0)
analysis = pd.read_csv(args.snv_analysis, sep='\t', header=0)

keep = analysis[['ORF_type','mut_seq', 'wt_seq','mut_len', 'wt_len', 'variant_ORF_id']]
annot = annot.merge(keep, on='variant_ORF_id', how='inner')
annot = annot[annot['variant_qc'] =='PASS']


def mut_protein_loc(row):
    if row['mut_len'] == row['wt_len']:
        mut_prot_loc = int((row['mut_trans_loc'] - row['ORF_start_trans_loc']) / 3)
        wt_aa = row['wt_seq'][mut_prot_loc]
        mut_aa = row['mut_seq'][mut_prot_loc]
        row['prot_mut'] = str(wt_aa) + str(mut_prot_loc + 1) + str(mut_aa)
    if row['mut_len'] < row['wt_len']:
        row['prot_mut'] = 'premature_stop'
    if row['mut_len'] > row['wt_len']:
        row['prot_mut'] = 'mutated_stop'
    return(row)

annot = annot.apply(mut_protein_loc, axis=1)
annot = annot.drop_duplicates(subset='mut_seq')

annot[['ORF_ID','variant_ORF_id','trans_variant_id', 'wt', 'mut', 'variant_qc', 'ORF_type', 'prot_mut']].to_csv(args.output, sep='\t', header=True, index=None)
