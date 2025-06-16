#!/usr/bin/env bash

# coding_CDS > coding_UTR > coding_intron > non-coding_exon > non-coding_intron > intergenic

set -e -u -o pipefail

### gene annotation ###
mkdir -p fixed_numts polymorphic_numts
for type in coding_CDS coding_UTR coding_intron non-coding_exon non-coding_intron
do
  bedtools intersect -a ${fixed_nupos} -b chm13_gene_annotation_regions/T2T-CHM13.${type}.merged.bed -wa | cut -f4 | sort -u > fixed_numts/${type}.numts_id.list
done

for type in coding_CDS coding_UTR coding_intron non-coding_exon non-coding_intron
do
  bedtools intersect -a ${polymorphic_nupos} -b chm13_gene_annotation_regions/T2T-CHM13.${type}.merged.bed -wa | cut -f4 | sort -u > polymorphic_numts/${type}.numts_id.list
done

# shuffle
mkdir -p fixed_numts/shuffle_results

parallel --jobs 10 "
for num in {0..9999}
do
  bedtools intersect -a \"${fixed_shuffled}/fixed_numts.nu_pos.shuffled.\${num}.bed\" -b \"chm13_gene_annotation_regions/T2T-CHM13.{1}.merged.bed\" -wa | cut -f4 | sort -u | wc -l >> fixed_numts/shuffle_results/{1}.numts.list
done
" ::: coding_CDS coding_UTR coding_intron non-coding_exon non-coding_intron

mkdir -p polymorphic_numts/shuffle_results

parallel --jobs 10 "
for num in {0..9999}
do
  bedtools intersect -a "${polymorphic_shuffled}/polymorphic_numts.nu_pos.shuffled.\${num}.bed\" -b \"chm13_gene_annotation_regions/T2T-CHM13.{1}.merged.bed\" -wa | cut -f4 | sort -u | wc -l >> polymorphic_numts/shuffle_results/{1}.numts.list
done
" ::: coding_CDS coding_UTR coding_intron non-coding_exon non-coding_intron


### repeats ###
mkdir -p repeats_classification
parallel --jobs 16 '
for i in {0..9999}
do
  bedtools slop -g ${chrom_length} -b 100 -i {2}_shuffled/{2}.nu_pos.shuffled.${i}.bed | \
  bedtools intersect -a - -b T2T-CHM13.{1}.repeatmasker.merged.bed -wa | \
  sort -u | cut -f4 | wc -l >> repeats_classification/{2}.{1}.simulation.each100bp.list
done
' ::: DNA LINE Low_complexity LTR RC Retroposon rRNA Satellite scRNA Simple_repeat SINE snRNA srpRNA tRNA ::: fixed_numts polymorphic_numts
