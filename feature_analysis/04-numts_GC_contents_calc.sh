#!/usr/bin/env bash

set -e -u -o pipefail

mkdir -p body_regions

### fixed numts ###
# observed
bedtools getfasta -nameOnly -fi ${chm13_fasta} -bed ${fixed_nu_pos_bed} | ../GC_content.py - | awk '{count+=1; sum+=$2}END{print sum/count}'

# nuclear genome shuffled
mkdir -p nugenome_shuffled
parallel --jobs 100 '
bedtools shuffle -i '"${fixed_nu_pos_bed}"' -g '"${chm13_len}"' -excl '"${chm13_censhort_bed}"' | bedtools getfasta -nameOnly -fi '"${chm13_fasta}"' -bed - | ./GC_content.py - > nugenome_shuffled/fixed_numts.nugenome_shuffled.body_regions.{1}.tsv
' ::: {1..10000}

for i in nugenome_shuffled/*
do
  awk '{count+=1; sum+=$2}END{print sum/count}' $i >> fixed_numts.GC_content.nugenome_shuffled.body_regions.results.list
done

# mitochondria genome shuffled
mkdir -p mtgenome_shuffled
parallel --jobs 100 '
bedtools shuffle -i '"${fixed_nu_pos_bed}"' -g <(echo -e "human_chrM_1_chrM_combined\t33138") | bedtools getfasta -nameOnly -fi human_chrM_1.chrM_combined.fasta -bed - | ./GC_content.py - > mtgenome_shuffled/fixed_numts.mtgenome_shuffled.body_regions.{1}.tsv
' ::: {1..10000}


### polymorphic numts ###
bedtools getfasta -nameOnly -fi polymorphic_numts.pannumts.ALT_sequences.fasta -bed polymorphic_numts.pannumts.numts_sequences.bed > polymorphic_numts.pannumts.numts_sequences.fasta

bedtools getfasta -nameOnly -fi ${chm13_fasta} -bed polymorphic_numts.refnumts.numts_sequences.bed > polymorphic_numts.refnumts.numts_sequences.fasta

cat polymorphic_numts.pannumts.numts_sequences.fasta polymorphic_numts.refnumts.numts_sequences.fasta > polymorphic_numts.observed.body_regions.numts_sequences.fasta

# observed
./GC_content.py polymorphic_numts.pannumts.numts_sequences.fasta > polymorphic_numts.GC_content.observed.body_regions.distribution.tsv

# nuclear genome shuffled
mkdir -p nugenome_shuffled
parallel --jobs 100 '
bedtools shuffle -i sequences/polymorphic_numts.observed.body_regions.numts_sequences.bed -g '"${chm13_len}"' -excl '"${chm13_censhort_bed}"' | bedtools getfasta -nameOnly -fi '"${chm13_fasta}"' -bed - | ./GC_content.py - > nugenome_shuffled/polymorphic_numts.GC_content.nugenome_shuffled.body_regions.{1}.tsv
' ::: {1..10000}

for i in nugenome_shuffled/*
do
  awk '{count+=1; sum+=$2}END{print sum/count}' $i >> polymorphic_numts.GC_content.nugenome_shuffled.body_regions.results.list
done

# mitochondria genome shuffled
mkdir -p mtgenome_shuffled
parallel --jobs 100 '
bedtools shuffle -i sequences/polymorphic_numts.observed.body_regions.numts_sequences.bed -g <(echo -e "human_chrM_10x\t165690") | bedtools getfasta -nameOnly -fi sequences/human_chrM_10x.fasta -bed - | ./GC_content.py - > mtgenome_shuffled/polymorphic_numts.GC_content.mtgenome_shuffled.body_regions.{1}.tsv
' ::: {1..10000}

for i in mtgenome_shuffled/*
do
  awk '{count+=1; sum+=$2}END{print sum/count}' $i >> polymorphic_numts.GC_content.mtgenome_shuffled.body_regions.results.list
done
