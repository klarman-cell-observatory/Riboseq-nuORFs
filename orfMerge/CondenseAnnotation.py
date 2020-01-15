# Python 3.6
'''
# Author: Travis Law
# Date: 07/23/2018
# Objective: Condense Annotation Files
'''
import pandas


def split(df):
    plus_strand = df[df['#strand'] == '+']
    minus_strand = df[df['#strand'] == '-']
    return(plus_strand, minus_strand)


def filter_plus_str(df):
    df = df.sort_values(['transcriptID', '#end', '#orf_len'],
                        ascending=[True, True, False])
    df['transcript_end'] = (df["transcriptID"].map(str) +
                            '_' +
                            df["#end"].map(str))
    unique_plus_strand = df.drop_duplicates(subset=['transcript_end'],
                                            keep='first').copy()
    unique_plus_strand.drop(['transcript_end'],
                            axis=1,
                            inplace=True)
    return(unique_plus_strand)


def filter_minus_str(df):
    df = df.sort_values(['transcriptID', '#start', '#orf_len'],
                        ascending=[True, False, False])
    df['transcript_start'] = (df["transcriptID"].map(str) +
                              '_' +
                              df["#start"].map(str))
    unique_minus_strand = df.drop_duplicates(subset=['transcript_start'],
                                             keep='first').copy()
    unique_minus_strand.drop(['transcript_start'],
                             axis=1,
                             inplace=True)
    return(unique_minus_strand)


def filter_riborf(riborf_output):
    plus_strand, minus_strand = split(riborf_output)
    unique_plus = filter_plus_str(plus_strand)
    unique_minus = filter_minus_str(minus_strand)
    frames = [unique_plus, unique_minus]
    unique = pandas.concat(frames)
    unique = unique.sort_index()
    return(unique)


def mergeCanonical(canORF, allORF):
    frames = [canORF,
              allORF[allORF['#orf_type'] != 'canonical']]
    df = pandas.concat(frames)
    df = df.sort_index()
    return(df)


def saveBed(df, out):
    cols = ['#seqname',
            '#start',
            '#end',
            'ORF_ID',
            '#score',
            '#strand',
            '#thick_start',
            '#thick_end',
            '#color',
            '#num_exons',
            '#exon_lengths',
            '#exon_genomic_relative_starts']
    BED12 = df[cols]
    cols = ['#chr',
            '#start',
            '#end',
            '#id',
            '#score',
            '#strand',
            '#thick_start',
            '#thick_end',
            '#color',
            '#num_exons',
            '#exon_lengths',
            '#exon_genomic_relative_starts']
    BED12.columns = cols
    BED12.to_csv(out,
                 index=None,
                 header=True,
                 sep='\t')


ANNOTATION = '/ahg/regevdata/projects/Ribo-seq/MHC-I/ref/orfeome/RNA.Unfiltered/B721.RNA.mit.annot'  # noqa
BED = '/ahg/regevdata/projects/Ribo-seq/MHC-I/ref/orfeome/RNA.Unfiltered/B721.RNA.mit.cond.bed'  # noqa
INFO = '/ahg/regevdata/projects/Ribo-seq/MHC-I/ref/orfeome/RNA.Unfiltered/B721.RNA.mit.cond.info'  # noqa
annotations = pandas.read_csv(
    ANNOTATION,
    sep='\t',
    index_col=False,
    header=0,
    low_memory=False)
annotations['ORF_ID'] = (
    annotations['gene_name'].map(str) + '|' +
    annotations['#id'].map(str) + '|' +
    annotations['#orf_type'].map(str))
canonical = annotations[annotations['#orf_type'] == 'canonical']
unique_all = filter_riborf(annotations)
unique_canonical = filter_riborf(canonical)
comb_unique = mergeCanonical(unique_canonical, unique_all)
saveBed(comb_unique, BED)
comb_unique.to_csv(
    INFO,
    index=None,
    header=True,
    sep='\t')
