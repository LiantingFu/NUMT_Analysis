#!/usr/bin/env bash

set -e -u -o pipefail

# get mtregions
# ALL; fixed; polymorphic

all_info="/home/ltfu/APG_project/NUMT/c-NUMT_genomic_features/a-numts_classification/final_classification/human_numts.info.tsv"
fixed_info="/home/ltfu/APG_project/NUMT/c-NUMT_genomic_features/a-numts_classification/final_classification/excl_chrY/fixed_numts.info.excl_chrY.tsv"
polymorphic_info="/home/ltfu/APG_project/NUMT/c-NUMT_genomic_features/a-numts_classification/final_classification/polymorphic_numts.info.tsv"
uncomparable_info="/home/ltfu/APG_project/NUMT/c-NUMT_genomic_features/a-numts_classification/final_classification/excl_chrY/uncomparable_numts.info.excl_chrY.tsv"

### mtdna region ###
# ALL
cut -f6-8 ${all_info} | awk 'BEGIN{OFS=FS="\t"} NR>1 {if($3>16569){print $1,$2,16569"\n"$1,1,$3-16569}else{print $0}}' > all_numts.mtregions.bed

# fixed
cut -f6-8 ${fixed_info} | awk 'BEGIN{OFS=FS="\t"} NR>1 {if($3>16569){print $1,$2,16569"\n"$1,1,$3-16569}else{print $0}}' > fixed_numts_excl_chrY.mtregions.bed

# polymorphic
cut -f6-8 ${polymorphic_info} | awk 'BEGIN{OFS=FS="\t"} NR>1 {if($3>16569){print $1,$2,16569"\n"$1,1,$3-16569}else{print $0}}' > polymorphic_numts.mtregions.bed

# uncomparable
cut -f6-8 ${uncomparable_info} | awk 'BEGIN{OFS=FS="\t"} NR>1 {if($3>16569){print $1,$2,16569"\n"$1,1,$3-16569}else{print $0}}' > pericentric_numts_excl_chrY.mtregions.bed

### coverage ###
bedtools genomecov -i fixed_numts_excl_chrY.mtregions.bed -g primates_calc/human_mtdna_length.txt -d > human_mtdna_cov/fixed_numts_excl_chrY.mtdna_cov.tsv

bedtools genomecov -i polymorphic_numts.mtregions.bed -g primates_calc/human_mtdna_length.txt -d > human_mtdna_cov/polymorphic_numts.mtdna_cov.tsv

bedtools genomecov -i pericentric_numts_excl_chrY.mtregions.bed -g primates_calc/human_mtdna_length.txt -d > human_mtdna_cov/pericentric_numts_excl_chrY.mtdna_cov.tsv
