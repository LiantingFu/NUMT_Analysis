#!/usr/bin/env bash

set -e -u -o pipefail

### fixed ###
mkdir -p alignments_results trimmed_results

awk 'NR > 1' ${fixed_numts_info} | while read -r ID nuchrom nustart nuend nulen mtchrom mtstart mtend mtlen strand
do
  nu_region="${nuchrom}\t${nustart}\t${nuend}\t${ID}.nu_seq\t.\t${strand}"
  mt_region="human_chrM_1_chrM_combined\t${mtstart}\t${mtend}\t${ID}.mt_seq\t.\t${strand}"
  ( \
    bedtools getfasta -fi ${nu_fasta} -bed <(echo -e "${nu_region}") -name -s; \
    bedtools getfasta -fi ${mt_fasta} -bed <(echo -e "${mt_region}") -name \
  ) | mafft --localpair --maxiterate 1000 --quiet - > alignments_results/${ID}.mafft_aln.txt

  trimal -in alignments_results/${ID}.mafft_aln.txt -out trimmed_results/${ID}.mafft_aln.trimmed.txt -gt 0.8

  ./get_sequences_similarity.py ${ID} trimmed_results/${ID}.mafft_aln.trimmed.txt
done > fixed_numts.excl_chrY.alignment_results.tsv

awk 'BEGIN{OFS=FS="\t"}{print $0,$2/$3}' fixed_numts.alignment_results.tsv > fixed_numts.excl_chrY.seqs_similarity_results.tsv


### polymorphic ###
cat polymorphic_pannumts.seq_id.list | while read -r sample
do
  trimal -in alignments_results/${sample}.mafft_aln.txt -out trimmed_results/${sample}.mafft_aln.trimmed.txt -gt 0.8
done

cat polymorphic_pannumts.seq_id.list | while read -r sample
do
  ./get_sequences_similarity.py ${sample} trimmed_results/${sample}.mafft_aln.trimmed.txt
done > polymorphic_pannumts.alignment_results.tsv

awk 'BEGIN{OFS=FS="\t"}{split($1,a,"_"); print a[1],$2,$3}' polymorphic_pannumts.alignment_results.tsv | sort -k1,1 | bedtools groupby -g 1 -c 2,3 -o sum,sum | awk 'BEGIN{OFS=FS="\t"}{print $0,$2/$3}' > polymorphic_pannumts.seqs_similarity_results.tsv

## refnumts
cat ${refnumts_info} | while read -r ID nuchrom nustart nuend nulen mtchrom mtstart mtend mtlen strand
do
  nu_region="${nuchrom}\t${nustart}\t${nuend}\t${ID}.nu_seq\t.\t${strand}"
  mt_region="human_chrM_1_chrM_combined\t${mtstart}\t${mtend}\t${ID}.mt_seq\t.\t${strand}"
  ( \
    bedtools getfasta -fi ${nu_fasta} -bed <(echo -e "${nu_region}") -name -s; \
    bedtools getfasta -fi ${mt_fasta} -bed <(echo -e "${mt_region}") -name \
  ) | mafft --localpair --maxiterate 1000 --quiet - > alignments_results/${ID}.mafft_aln.txt

  trimal -in alignments_results/${ID}.mafft_aln.txt -out trimmed_results/${ID}.mafft_aln.trimmed.txt -gt 0.8

  ./get_sequences_similarity.py ${ID} trimmed_results/${ID}.mafft_aln.trimmed.txt
done > polymorphic_refnumts.alignment_results.tsv

awk 'BEGIN{OFS=FS="\t"}{print $0,$2/$3}' polymorphic_refnumts.alignment_results.tsv > polymorphic_refnumts.seqs_similarity_results.tsv
