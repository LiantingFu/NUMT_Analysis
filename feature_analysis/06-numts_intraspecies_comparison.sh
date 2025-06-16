#!/usr/bin/env bash

#SBATCH --job-name=numts_comparison
#SBATCH --partition=debug64c512g
#SBATCH --cpus-per-task 2
#SBATCH --ntasks-per-node 1
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err

set -e -u -o pipefail

species=$1


##### SPECIES aln HUMAN #####
mkdir -p species_aln_human
mkdir -p species_aln_human/${species}

## liftOver
liftOver ${human_numts_bed} ${species_aln_human_chain} species_aln_human/${species}/human_numts.${species}.liftOver.bed species_aln_human/${species}/human_numts.${species}.unmapped.bed

## classification
./post_liftover.py ${species} query species_aln_human/${species}/human_numts.${species}.unmapped.bed ${species_aln_human_bam}

## duplication
liftOver -multiple species_aln_human/${species}/human_numts.${species}.duplication.bed ${species_aln_human_chain} species_aln_human/${species}/human_numts.${species}.duplication.reliftOver.bed species_aln_human/${species}/unmapped

mkdir -p species_aln_human/human_specific_results

awk '$7=="deletion"' species_aln_human/${species}/human_numts.${species}.deletion_results.tsv > species_aln_human/human_specific_results/${species}_aln_human.deletion_output.tsv

echo "${species}_aln_human complete!"


### HUMAN aln SPECIES #####
mkdir -p human_aln_species
mkdir -p human_aln_species/${species}

# liftOver
liftOver ${species_numts_bed} ${human_aln_species_chain} human_aln_species/${species}/${species}_numts.human.liftOver.bed human_aln_species/${species}/${species}_numts.human.unmapped.bed

## classification
./post_liftover.py ${species} ref human_aln_species/${species}/${species}_numts.human.unmapped.bed ${human_aln_species_bam}

# duplication
liftOver -multiple human_aln_species/${species}/${species}_numts.human.duplication.bed ${human_aln_species_chain} human_aln_species/${species}/${species}_numts.human.duplication.reliftOver.bed human_aln_species/${species}/unmapped

rm human_aln_species/${species}/unmapped

mkdir -p human_aln_species/species_specific_results
awk '$7=="deletion"' human_aln_species/${species}/${species}_numts.human.deletion_results.tsv > human_aln_species/species_specific_results/${species}_numts.human.deletion_output.tsv

echo "human_aln_${species} complete!"


mkdir -p all_human_numts
mkdir -p all_human_numts/${species}

## liftOver
liftOver ${all_human_numts_bed} ${species_aln_human_chain} all_human_numts/${species}/human_numts.${species}.liftOver.bed all_human_numts/${species}/human_numts.${species}.unmapped.bed

## classification
./post_liftover.py ${species} all all_human_numts/${species}/human_numts.${species}.unmapped.bed ${species_aln_human_bam}

## duplication
liftOver -multiple all_human_numts/${species}/human_numts.${species}.duplication.bed ${species_aln_human_chain} all_human_numts/${species}/human_numts.${species}.duplication.reliftOver.bed all_human_numts/${species}/unmapped

mkdir -p all_human_numts/human_specific_results

awk '$7=="deletion"' all_human_numts/${species}/human_numts.${species}.deletion_results.tsv > all_human_numts/human_specific_results/${species}_aln_human.deletion_output.tsv

echo "${species}_aln_human all numts complete!"
