
'''
# Author: Tamara Ouspenskaia
# Date: 01/14/2020
# Objective: This script combines a variant VCF file with an ORF database
'''



#! /bin/bash
#$ -cwd
#$ -l h_vmem=50g
#$ -l h_rt=02:00:00

# Inputs:

OUT_DIR=$1
mutect_vcf=$2

mkdir -p $OUT_DIR


# Outputs:

trans_mutect_variants=$OUT_DIR/"trans.mutect2.txt"

#Follow these steps to obtain the references needed for variant incorporation into ORFs:

#### Combine GENCODE + MiTranscriptome transcripts BED files.
#### Convert transcript BED12 to GenePred
#### Extract the columns containing the transcript ID and genomic coordinates of exon boundaries.

# cat ../Resources/RPBP/human_g1k_v37/human_g1k_v37.annotated.thickfixed.bed \
# ../RPBP/shortgen_mi_unannot/mi_unannot.annotated.bed \
# > ../Analysis/variants/gencode_mi_unannot_transcripts.bed
#
# ../Scripts/bedToGenePred \
# ../Analysis/variants/gencode_mi_unannot_transcripts.bed \
# ../Analysis/variants/gencode_mi_unannot_transcripts.genepred
#
# cut -f1,9,10 ../Analysis/variants/gencode_mi_unannot_transcripts.genepred > ../Analysis/variants/gencode_mi_unannot_transcripts.genomic_exons.txt

trans_bed=../Analysis/variants/gencode_mi_unannot_transcripts.bed
trans_genepred=../Analysis/variants/gencode_mi_unannot_transcripts.genomic_exons.txt



#### Generate a fasta file for all transcripts and make sure that IDs match the BED12 file exactly:
## Do not need to re-run this every time.

# bedtools getfasta -name -s -split -fi ../references/b37/human_g1k_v37.fasta -bed \
# ../Analysis/variants/gencode_mi_unannot_transcripts.bed > \
# ../Analysis/variants/gencode_mi_unannot_transcripts.fasta

trans_fasta=../gencode_mi_unannot_transcripts.stripped.fasta


#### Intersect VCF with BED12 to get transcripts with mutations:

##Mutect2 variants (unfiltered):

bedtools intersect -wo -split -a $mutect_vcf -b $trans_bed > $trans_mutect_variants


#### Now ready to incorporate variants into ORFs.
#This script requires the following files:
# - Subset of genepred of transcripts generated above containing genomic coordinates of exon boundaries.
# - BED12 of ORFs
# - VCF/BED12 intersect between variants and transcripts BED12 generated above.
# - Fasta of transcripts where headers match ID's in the transcripts BED12 file generated above.

orf_bed=../nuORFdb.bed

python ../variant_to_protein_gatk_only.py \
--mutect_variants $trans_mutect_variants \
--trans_genepred $trans_genepred \
--orf_bed $orf_bed \
--trans_fasta $trans_fasta \
--splice_site_indels $OUT_DIR/"splice_site_indels.txt" \
--mut_trans_fasta_snv $OUT_DIR/"transcripts.snv.fasta" \
--mut_trans_fasta_indel $OUT_DIR/"transcripts.indel.fasta" \
--mut_orf_fasta_snv $OUT_DIR/"ORFs.snv.fasta" \
--mut_orf_fasta_indel $OUT_DIR/"ORFs.indel.fasta" \
--mut_annot_snv $OUT_DIR/"ORFs.snv.annot" \
--mut_annot_indel $OUT_DIR/"ORFs.indel.annot"

#Variant analysis - SNV


orf_fasta=../nuORFdb.fasta

python ../variant_analysis.py \
--wt_orfs $orf_fasta \
--mut_orfs $OUT_DIR/"ORFs.snv.fasta" \
--mut_annot $OUT_DIR/"ORFs.snv.annot" \
--filtered_fasta $OUT_DIR/"filtered.snv.fasta" \
--all_variants_counts $OUT_DIR/"all_snv_counts.tab" \
--high_qual_variants_counts $OUT_DIR/"hiqh_qual_snv_counts.tab" \
--indel_analysis $OUT_DIR/"snv.analysis.tab"

#Variant analysis - Indels

python ../variant_analysis.py \
--wt_orfs $orf_fasta \
--mut_orfs $OUT_DIR/"ORFs.indel.fasta" \
--mut_annot $OUT_DIR/"ORFs.indel.annot" \
--filtered_fasta $OUT_DIR/"filtered.indel.fasta" \
--all_variants_counts $OUT_DIR/"all_indel_counts.tab" \
--high_qual_variants_counts $OUT_DIR/"hiqh_qual_indel_counts.tab" \
--indel_analysis $OUT_DIR/"indel.analysis.tab"

#Split FASTA into peptides - snv

OUT_DIR_PEP=$OUT_DIR/peptides/snv
mkdir -p $OUT_DIR_PEP


PEP_OUT_FILE=$OUT_DIR_PEP/"snv_peptides_all.txt.gz"
PEP_FINAL=$OUT_DIR_PEP/"snv_peptides.txt.gz"

python ../split_fasta.py --fasta_file $OUT_DIR/"filtered.snv.fasta" --shortest_peptide 9 --longest_peptide 10 | gzip - > $PEP_OUT_FILE

cd $OUT_DIR_PEP

zcat $PEP_OUT_FILE | sort | uniq | gzip - > $PEP_FINAL

#Split FASTA into peptides - indels

OUT_DIR_PEP=$OUT_DIR/pepides/indel
mkdir -p $OUT_DIR_PEP

PEP_OUT_FILE=$OUT_DIR_PEP/"indel_peptides_all.txt.gz"
PEP_FINAL=$OUT_DIR_PEP/"indel_peptides.txt.gz"

python ../split_fasta.py --fasta_file $OUT_DIR/"filtered.indel.fasta" --shortest_peptide 9 --longest_peptide 10 | gzip - > $PEP_OUT_FILE

cd $OUT_DIR_PEP

zcat $PEP_OUT_FILE | sort | uniq | gzip - > $PEP_FINAL

#Mutation protein coordinates

python ../mutation_protein_coordinates.py \
--snv_annot $OUT_DIR/"ORFs.snv.annot" \
--snv_analysis $OUT_DIR/"snv.analysis.tab" \
--output $OUT_DIR/"snvs.analysis"


#qsub ../variant_pipeline.sh /path/to/output /path/to/variants.vcf
