'''
# Author: Tamara Ouspenskaia
# Date: 01/14/2020
# Objective: This script incorporates variants into ORFs
'''




#This script requires the following files:
# - Trucated genepred of transcripts
# - BED12 of ORFs
# - VCF/BED12 intersect between variants and transcripts BED12
# - Fasta of transcripts where headers match ID's in the transcripts BED12 file

import argparse
import pandas as pd
import numpy as np
from functools import reduce
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


parser = argparse.ArgumentParser()

parser.add_argument("--mutect_variants", required = True, help = "Path to VCF/BED12 intersect output from mutect2")
parser.add_argument("--orf_bed", required = True, help = "Path to the ORF BED12")
parser.add_argument("--trans_genepred", required = True, help = "Path to the truncated genepred of transcripts")
parser.add_argument("--trans_fasta", required = True, help = "Path to fasta of transcripts")
parser.add_argument("--mut_trans_fasta_snv", required = True, help = "Path to output fasta of mutated transcripts with SNVs")
parser.add_argument("--mut_trans_fasta_indel", required = True, help = "Path to output fasta of mutated transcripts with InDels")
parser.add_argument("--splice_site_indels", required = True, help = "Path to output file with ORFs where indels span exon-intron juctions")
parser.add_argument("--mut_orf_fasta_snv", required = True, help = "Path to output fasta of mutated ORFs with SNVs")
parser.add_argument("--mut_orf_fasta_indel", required = True, help = "Path to output fasta of mutated ORFs with InDels")
parser.add_argument("--mut_annot_snv", required = True, help = "Path to output info of mutated ORFs with SNVs")
parser.add_argument("--mut_annot_indel", required = True, help = "Path to output info of mutated ORFs with InDels")


args = parser.parse_args()

mutect_variants_unfilt = pd.read_csv(args.mutect_variants, header=None, sep='\t')

def variant_id(df):
    df['variant_id'] = df[14] + '|' + df[0].map(str) + '|' + df[1].map(str) + '|' + df[3] + '|' + df[4]
    return(df)

dfs = [mutect_variants_unfilt]
for df in dfs:
    df = variant_id(df)

unfilt_variants = mutect_variants_unfilt
unfilt_variants = unfilt_variants.drop_duplicates(subset='variant_id')
unfilt_variants = unfilt_variants.drop([2,5,7,8,9,10,15,17,18,19,23], axis=1)

cols = ['mut_chr','mut_location','wt','mut','variant_qc','trans_chr','trans_genomic_start','trans_genomic_end','trans_id','trans_strand','trans_num_exons','trans_exon_lengths','exon_start_transcript','variant_id']
unfilt_variants.columns = cols


# #### Load the BED12 for the ORFs:

cols = ['ORF_chr', 'ORF_genomic_start', 'ORF_genomic_end', 'ORF_ID', 'Score', 'ORF_strand',
       'thick_start', 'thick_end', 'color', 'num_exons',
       'exon_lengths', 'exon_genomic_relative_starts']
ORFs_bed = pd.read_csv(args.orf_bed, header=None, names=cols, sep='\t')
ORFs_bed = ORFs_bed[['ORF_genomic_start', 'ORF_genomic_end', 'ORF_ID']]


# #### Load the Genepred for transcripts:

cols = ['trans_id', 'trans_exon_starts_genomic', 'trans_exon_ends_genomic']
trans_genepred = pd.read_csv(args.trans_genepred, header=None, names=cols, sep='\t')


# #### Extract transcript id from ORF_ID in the ORF BED12 file:

def extract_trans_id(ORF_ID):
    ENST = ORF_ID.split('|')[1]
    trans_id = ENST.split(':')[0]
    trans_id = trans_id.split('_')[0:-1]
    trans_id = '_'.join(trans_id)
    return trans_id

ORFs_bed['trans_id'] = ORFs_bed['ORF_ID'].apply(extract_trans_id)


# #### Combine the necessary columns from variants, trans_bed, and trans_genepred.

unfilt_variants = unfilt_variants[['mut_chr', 'mut_location', 'wt', 'mut', 'variant_qc', 'trans_id',
       'trans_strand', 'trans_num_exons','trans_exon_lengths', 'variant_id']]

variants_genepred = unfilt_variants.merge(trans_genepred, on='trans_id', how='left')
merge = variants_genepred.merge(ORFs_bed, on='trans_id', how = 'left')
merge.dropna(subset=['ORF_ID'], inplace=True)

def fix_col(col):
    col = col.str.rstrip(',').str.split(',', expand = False)
    return(col)

merge['trans_exon_lengths'] = fix_col(merge['trans_exon_lengths'])
merge['trans_exon_starts_genomic'] = fix_col(merge['trans_exon_starts_genomic'])
merge['trans_exon_ends_genomic'] = fix_col(merge['trans_exon_ends_genomic'])

merge['trans_num_exons'] = merge['trans_num_exons'].astype(int)


# Not all of the variants fall within an ORF. Can further filter the table to enrich it for variants that fall within ORFs. I will just check if the mutation is between the genomic start and stop.

merge['variant_in_orf'] = np.where((merge['ORF_genomic_start'] <= merge['mut_location']) &
                               (merge['mut_location'] <= merge['ORF_genomic_end']) , 'in', 'not_in')

variant_in_orf = merge[merge['variant_in_orf'] == 'in']

# Determine if the variant is an SNV or an InDel:

def indel_or_snv(row):
    if len(row['wt']) == len(row['mut']):
        row['variant_type'] = 'SNV'
    else:
        row['variant_type'] = 'InDel'
    return(row)

variant_in_orf = variant_in_orf.apply(indel_or_snv, axis=1)

# #### Determine the position of the mutation within the transcript

def mut_trans_loc(row):
    if row['trans_strand'] == '+':
        trans_pos = 0
        for i in range(0,row['trans_num_exons']):
            if int(row['trans_exon_starts_genomic'][i]) <= int(row['mut_location']) <= int(row['trans_exon_ends_genomic'][i]):
                trans_pos = trans_pos + (int(row['mut_location']) - int(row['trans_exon_starts_genomic'][i]))
                break
            else:
                trans_pos = trans_pos + int(row['trans_exon_lengths'][i])
        row['mut_trans_loc'] = trans_pos - 1
    if row['trans_strand'] == '-':
        trans_pos = 0
        for i in range(0,row['trans_num_exons']):
            exon_start = int(list(reversed(row['trans_exon_ends_genomic']))[i])
            exon_end = int(list(reversed(row['trans_exon_starts_genomic']))[i])
            if exon_start >= int(row['mut_location']) >= exon_end:
                trans_pos = trans_pos + (exon_start - int(row['mut_location']))
                break
            else:
                trans_pos = trans_pos + int(list(reversed(row['trans_exon_lengths']))[i])
        row['mut_trans_loc'] = trans_pos
    return(row)

variants_mut_trans_loc = variant_in_orf.apply(mut_trans_loc, axis=1)

def trans_variant_id(df):
    df['trans_variant_id'] = df['trans_id'] + '_' + df['mut_chr'].map(str) + df['trans_strand'] + ':' + df['mut_location'].map(str) + ':' + df['wt'] + '>' + df['mut']
    return(df)

variants_mut_trans_loc = trans_variant_id(variants_mut_trans_loc)
variants_mut_trans_loc = variants_mut_trans_loc.reset_index(drop=True)

# Split the annotation into InDels and SNVs:

snv_annot = variants_mut_trans_loc[variants_mut_trans_loc['variant_type'] == 'SNV']
snv_annot = snv_annot.reset_index(drop=True)
indel_annot = variants_mut_trans_loc[variants_mut_trans_loc['variant_type'] == 'InDel']
indel_annot = indel_annot.reset_index(drop=True)

# #### Apply the mutation to the transcript in the fasta file:

# This generates a fasta of unique transcript+variant combinations.

def correct_fasta(df, input_fasta, output_fasta, splice_mutants):
    record_dict = SeqIO.index(input_fasta, 'fasta')
    mut_trans_dict = {}
    with open(splice_mutants, 'a') as splice_file:
        for i in range(0, len(df)):
            trans_id = df.loc[i,'trans_id']
            trans_variant_id = df.loc[i,'trans_variant_id']
            seq = str(record_dict[trans_id].seq)
            mut_pos = int(df.loc[i,'mut_trans_loc'])
            wt = df.loc[i,'wt']
            mut = df.loc[i,'mut']
            if df.loc[i,'trans_strand'] == '+':
                if wt == seq[mut_pos : mut_pos+len(wt)]:
                    new_seq = seq[:mut_pos] + mut + seq[mut_pos+len(wt):]
                else:
                    splice_file.write("Error: the WT sequence on the + strand does not match" + '\n')
                    splice_file.write(trans_variant_id + '\n')
                    splice_file.write(df.loc[i,'ORF_ID'] + '\n')
            if df.loc[i,'trans_strand'] == '-':
                if Seq(wt).reverse_complement() == seq[mut_pos-len(wt)+1 : mut_pos+1]:
                    new_seq = seq[:mut_pos-len(wt)+1] + Seq(mut).reverse_complement() + seq[mut_pos+1:]
                else:
                    splice_file.write("Error: the WT sequence on the - strand does not match" + '\n')
                    splice_file.write(trans_variant_id + '\n')
                    splice_file.write(df.loc[i,'ORF_ID'] + '\n')
            if trans_variant_id not in mut_trans_dict:
                mut_trans_dict[trans_variant_id] = new_seq
            if trans_variant_id in mut_trans_dict:
                pass
    with open(output_fasta, "w") as f:
            for trans_variant_id in mut_trans_dict:
                f.write('>' + trans_variant_id + '\n')
                f.write(str(mut_trans_dict[trans_variant_id]) + '\n')

correct_fasta(snv_annot, args.trans_fasta, args.mut_trans_fasta_snv, args.splice_site_indels)
correct_fasta(indel_annot, args.trans_fasta, args.mut_trans_fasta_indel, args.splice_site_indels)

# #### Re-derive the ORFs:
# ### Convert ORF start from genomic to transcript coordinates:

def ORF_start_loc(row):
    if row['trans_strand'] == '+':
        trans_pos = 0
        for i in range(0,row['trans_num_exons']):
            if int(row['trans_exon_starts_genomic'][i]) <= int(row['ORF_genomic_start']) <= int(row['trans_exon_ends_genomic'][i]):
                trans_pos = trans_pos + (int(row['ORF_genomic_start']) - int(row['trans_exon_starts_genomic'][i]))
                break
            else:
                trans_pos = trans_pos + int(row['trans_exon_lengths'][i])
        row['ORF_start_trans_loc'] = trans_pos
    if row['trans_strand'] == '-':
        trans_pos = 0
        for i in range(0,row['trans_num_exons']):
            exon_start = int(list(reversed(row['trans_exon_ends_genomic']))[i])
            exon_end = int(list(reversed(row['trans_exon_starts_genomic']))[i])
            if exon_start >= int(row['ORF_genomic_end']) >= exon_end:
                trans_pos = trans_pos + (exon_start - int(row['ORF_genomic_end']))
                break
            else:
                trans_pos = trans_pos + int(list(reversed(row['trans_exon_lengths']))[i])
        row['ORF_start_trans_loc'] = trans_pos
    return(row)

snv_annot = snv_annot.apply(ORF_start_loc, axis=1)
indel_annot = indel_annot.apply(ORF_start_loc, axis=1)

# #### Get the ORF starting with the defined start, but looking for a new stop.
# Given a transcript id, subset for ORFs that are on that ID. Then, for each ORF, derive its sequence and print to a new fasta.

def variant_ORF_id(df):
    df['variant_ORF_id'] = df['ORF_ID'] + '|' + df['mut_chr'].map(str) + '|' + df['mut_location'].map(str) + '|' + df['wt'] + '|' + df['mut']
    return(df)

snv_annot = variant_ORF_id(snv_annot)
indel_annot = variant_ORF_id(indel_annot)

def find_ORF(df, input_fasta, output_fasta):
    with open(output_fasta, "w") as output:
        for rec in SeqIO.parse(input_fasta, 'fasta'):
            ORFs = df[df['trans_variant_id'] == rec.id]
            ORFs.reset_index(inplace=True, drop=True)
            for i in range(0, len(ORFs)):
                ORF_start = ORFs.loc[i,'ORF_start_trans_loc']
                record = SeqRecord(seq = rec.seq[ORF_start:].translate(to_stop=True), id = ORFs.loc[i,'variant_ORF_id'], description = '')
                SeqIO.write(record, output, "fasta")

find_ORF(snv_annot, args.mut_trans_fasta_snv, args.mut_orf_fasta_snv)
find_ORF(indel_annot, args.mut_trans_fasta_indel, args.mut_orf_fasta_indel)

#### Save the annotation file

snv_annot[['mut_chr','mut_location','wt','mut','variant_qc','trans_id','trans_strand','variant_id','ORF_ID','mut_trans_loc','trans_variant_id','ORF_genomic_start','ORF_genomic_end','ORF_start_trans_loc','variant_ORF_id']].to_csv(args.mut_annot_snv, header=True, index=None, sep='\t')
indel_annot[['mut_chr','mut_location','wt','mut','variant_qc','trans_id','trans_strand','variant_id','ORF_ID','mut_trans_loc','trans_variant_id','ORF_genomic_start','ORF_genomic_end','ORF_start_trans_loc','variant_ORF_id']].to_csv(args.mut_annot_indel, header=True, index=None, sep='\t')
