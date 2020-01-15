'''
# Author: Travis Law
# Date: 02/11/2019
# Objective: Calculate frame aware TPM for RPF
# Docker Image: tlaw/tpm:latest
# Expects
	bed12 file of all ORFs
	bam file of RibORF Offset-correct alignments
	path to output file
	path for temporary files
'''
import argparse
import pysam
import subprocess
import pandas
from tqdm import tqdm
NAMES = [
    'chr',
    'start',
    'end',
    'ORF_ID',
    'score',
    'strand',
    'tStart',
    'tEnd',
    'rgb',
    'numE',
    'lenE',
    'startE'
]
bedTypes = {
    'chr': str,
    'start': int,
    'end': int,
    'ORF_ID': str,
    'score': int,
    'strand': str,
    'tStart': int,
    'tEnd': int,
    'rgb': int,
    'numE': int,
    'lenE': str,
    'startE': str
}


def getInputs():
    print('Get CLI')
    parser = argparse.ArgumentParser(
        description='Calculate FPKM from bedtools coverage bed',
        epilog='Expects as input the results of bedtools coverage on a bed12',
        add_help=True,
        allow_abbrev=True)
    parser.add_argument(
        '-bam',
        metavar='Bamfile',
        help='Path to Bam',
        required=True,
        type=str,
        dest='bam')
    parser.add_argument(
        '-bed',
        metavar='Bedfile',
        help='Path to ORF regions',
        required=True,
        type=str,
        dest='bed')
    parser.add_argument(
        '-out',
        metavar='Output',
        help='Path to output',
        default='',
        type=str,
        dest='out')
    parser.add_argument(
        '-temp',
        metavar='Temp',
        help='Temporary Bam File Location',
        default='',
        type=str,
        dest='temp')
    args = parser.parse_args()
    if args.out == '':
        args.out = args.bed
    if args.temp == '':
        args.temp = 'temp.bam'
    print(args)
    return(args)


def prepBam(bam, temp):
    print('Prepare Bam')
    pysam.sort('-o', temp, bam)
    pysam.index(temp)


def cleanup(temp):
    print('Cleanup Temp Files')
    subprocess.call(['rm', temp])
    subprocess.call(['rm', temp + '.bai'])


def calculateRPK(row):
    return((row.loc['inFrameReads'] * 1000.0) / row.loc['length'])


def calculateTPM(row, scalingFactor):
    return(row.loc['RPK'] / scalingFactor)


def calculateLength(row):
    return(len(row.loc['coverage'].split(',')))


def calculateTotalReads(coverage):
    reads = 0
    for i in list(map(int, coverage.split(','))):
        reads += i
    return(reads)


def calculateInFrameReads(coverage):
    reads = 0
    mod = 0
    for i in list(map(int, coverage.split(','))):
        mod += 1
        if mod == 1:
            reads += i
        if mod == 3:
            mod = 0
    return(reads)


def calculatePurity(row):
    if row.loc['totalReads'] == 0:
        return(0.0)
    return(row.loc['inFrameReads'] / row.loc['totalReads'])


def getCoverage(row, bam):
    chromosome = str(row.loc['chr'])
    start = str(row.loc['start'])
    end = str(row.loc['end'])
    region = chromosome + ':' + start + '-' + end
    output = pysam.depth(
        '-aa',
        '-d',
        '0',
        '-r',
        region,
        bam
    )
    coverage = dict()
    for line in output.split('\n')[:-1]:
        line = line.split('\t')
        coverage[int(line[1])] = int(line[2])
    return(coverage)


def getPositions(row):
    start = row.loc['start']
    numE = row.loc['numE']
    lenE = row['lenE'].split(',')
    if lenE[-1] == '':
        del lenE[-1]
    lenE = list(map(int, lenE))
    startE = row['startE'].split(',')
    if startE[-1] == '':
        del startE[-1]
    startE = list(map(int, startE))
    positions = []
    for e in range(0, numE):
        for i in range(0, lenE[e]):
            pos = start + startE[e] + i
            positions.append(pos)
    return(positions)


def getTranscriptCoverage(row, bam):
    coverage = getCoverage(row, bam)
    transcript = getPositions(row)
    transcriptCoverage = []
    for pos in transcript:
        transcriptCoverage.append(coverage[pos])
    transcriptCoverage = ','.join(list(map(str, transcriptCoverage)))
    return(transcriptCoverage)


def main():
    print('Calculate RPF-TPM')
    args = getInputs()
    prepBam(args.bam, args.temp)
    bed = pandas.read_csv(
        args.bed,
        sep='\t',
        header=None,
        index_col=False,
        names=NAMES,
        dtype=bedTypes)
    tqdm.pandas(desc='Measure Coverage')
    bed['coverage'] = bed.progress_apply(
        getTranscriptCoverage,
        axis=1,
        args=(args.temp,)
    )
    tqdm.pandas(desc='Calculate Length')
    bed['length'] = bed.progress_apply(
        calculateLength,
        axis=1
    )
    tqdm.pandas(desc='Calculate In Frame Reads')
    bed['inFrameReads'] = bed['coverage'].progress_apply(
        calculateInFrameReads
    )
    tqdm.pandas(desc='Calculate Reads')
    bed['totalReads'] = bed['coverage'].progress_apply(
        calculateTotalReads
    )
    tqdm.pandas(desc='Calculate Purity')
    bed['Purity'] = bed.progress_apply(
        calculatePurity,
        axis=1
    )
    tqdm.pandas(desc='Calculate RPK')
    bed['RPK'] = bed.progress_apply(
        calculateRPK,
        axis=1
    )
    scalingFactor = bed['RPK'].sum() / 1000000.0
    tqdm.pandas(desc='Calculate TPM')
    bed['TPM'] = bed.progress_apply(
        calculateTPM,
        axis=1,
        args=(scalingFactor, )
    )
    bed.drop(
        ['coverage', 'RPK'],
        axis=1,
        inplace=True
    )
    print('Writing to ' + args.out)
    bed.to_csv(
        args.out,
        sep='\t',
        header=True,
        index=False)
    cleanup(args.temp)
    print('Done')


main()
