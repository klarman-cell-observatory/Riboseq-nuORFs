'''
# Author: Tamara Ouspenskaia
# Date: 01/14/2020
# Objective: This script will identify cancer-specific peptides, and couple them with coverage and translation information.
'''

#! /bin/bash
#$ -cwd
#$ -l h_vmem=50g
#$ -l h_rt=06:00:00


# Inputs:

mut_peptides=../snv_peptides.txt.gz #output from variant_pipeline.sh
reference_peptides=normal_peptides.txt.gz # list of normal peptides from the reference database
variant_annot=../snv.analysis.tab #output from variant_pipeline.sh
variant_coverage=../snv_coverage_across_samples_PASS.txt #output from coverage_across_samples.ipynb

# Output:

peptide_coverage=../peptide_coverage.tsv

python ../variant_peptides.py \
--mut_peptides $mut_peptides \
--reference_peptides $reference_peptides \
--variant_annot $variant_annot \
--variant_coverage $variant_coverage \
--peptide_coverage $peptide_coverage
