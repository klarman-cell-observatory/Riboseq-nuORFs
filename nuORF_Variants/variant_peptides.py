'''
# Author: Tamara Ouspenskaia
# Date: 01/14/2020
# Objective: This script will identify cancer-specific peptides, and couple them with coverage and translation information.
'''




# ### Goal:

# Generate unique peptides from a fasta of mutant ORFs. Search them against peptides from the original database. Re-assign them to their ORFs.
# First, generate all unique peptides for the fasta of mutant ORFs and PanSample using the previous script and load into this script.


import pandas as pd
import re
import gzip
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--mut_peptides", required = True,
                    help = "Path to fastq.gz file with unique peptides from the ORF fasta with variants")
parser.add_argument("--reference_peptides", required = True, help = "Path to fastq.gz file with unique peptides from the ORF fasta with variants")
parser.add_argument("--variant_annot", required = True, help = "Path to the analysis.tab file for variants")
parser.add_argument("--variant_coverage", required = True, help = "Path to the coverage file for variants")
parser.add_argument("--peptide_coverage", required = True, help = "Path to the peptide coverage OUTPUT file")


args = parser.parse_args()

# #### Extract peptides unique to cancer:

def load_peptides(gzip_file):
    with gzip.open(gzip_file) as f:
        peptide_list = f.readlines()
        peptide_list = [x.strip() for x in peptide_list]
        return(peptide_list)

snvs = load_peptides(args.mut_peptides)
pansample = load_peptides(args.reference_peptides)

mut_peptides = list(set(snvs) - set(pansample))
mut_peptides = [x.decode("utf-8") for x in mut_peptides]

# #### Assign peptides to their ORFs:

# Filter for the variants that pass GATK QC:

annot = pd.read_csv(args.variant_annot, header=0, sep='\t')

annot_pass = annot[annot['variant_qc'] == 'PASS']
annot_pass.dropna(axis=0, how='any', inplace=True)

index_series = []
for i in mut_peptides:
    annot_contains_index = annot_pass[annot_pass['mut_seq'].str.contains(i)].index
    index_series.append(annot_contains_index)

annot_pass_index = pd.DataFrame({'annot_index_mut': index_series,
                                    'peptides': mut_peptides,
                                    })

def nonzero_index(col):
    return(len(col))

annot_pass_index['mut_index_len'] = annot_pass_index['annot_index_mut'].apply(nonzero_index)


# Only keep peptides that were mapped to ORFs with with PASS filter variants:

annot_pass_index = annot_pass_index[annot_pass_index['mut_index_len'] != 0]


# #### For each mutant peptide, get the ORF types that it can be found in.

def get_orf_type(col):
    return(list(annot_pass.loc[col,'ORF_type']))

annot_pass_index['orf_types'] = annot_pass_index['annot_index_mut'].apply(get_orf_type)


# If a mutant peptide can be found in a canonical or CDS ORF, count it as canonical. Otherwise, could it as nuORF.

def test_canonical(col):
    if (('canonical' in col) | ('CDS' in col)):
        orf_cat = 'annotated'
    else:
        orf_cat = 'nuORF'
    return(orf_cat)

annot_pass_index['peptide_type'] = annot_pass_index['orf_types'].apply(test_canonical)


# ### Get the read coverage supporting each of these peptides.

coverage = pd.read_csv(args.variant_coverage, sep='\t', header=0)

# #### For each peptide, get the variant_ORF_ids from the annot table. For each variant_ORF_id, get the read coverage of the variant.

# Split each peptide into multiple rows, each row for a different variant_ORF_id that the peptide can be found in. Then, each variant_ORF_id can have its associated coverages and TPM.

def get_variant_orf_id(col):
    return(list(annot_pass.loc[col,'variant_ORF_id']))

annot_pass_index['variant_ORF_ids'] = annot_pass_index['annot_index_mut'].apply(get_variant_orf_id)


# #### Create a new dataframe, where for each peptide, there are as many lines as variant_ORF_ids, with each variant_ORF_id associated with SNV coverage.

def create_df(peptides_df):
    df = pd.DataFrame(columns=['peptide', 'mut_index_len', 'orf_types','peptide_type','variant_ORF_ids','variant_ORF_id'])
    for i in list(peptides_df.index):
        listofSeries = []
        for j in peptides_df.loc[i, 'variant_ORF_ids']:
            listofSeries.append(pd.Series([peptides_df.loc[i, 'peptides'],peptides_df.loc[i, 'mut_index_len'],peptides_df.loc[i, 'orf_types'],peptides_df.loc[i, 'peptide_type'],peptides_df.loc[i, 'variant_ORF_ids'],j],
                                          index=df.columns))
        df = df.append(listofSeries, ignore_index=True)
    return(df)

peptides_to_orfs = create_df(annot_pass_index)


# #### Merge with the coverage table to get the read coverage for each peptide:

pep_coverage = peptides_to_orfs.merge(coverage, on='variant_ORF_id', how='left')

pep_coverage.to_csv(args.peptide_coverage, sep='\t', header=True, index=False)
