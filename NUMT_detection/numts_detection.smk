#!/usr/bin/env python

import pandas as pd

configfile: "config.json"

chrM_df = pd.read_csv(config.get("chrM_info", "chrM.tsv"), header=0, index_col="species", sep="\t")
SCRIPT_DIR = config.get("script_dir", "config")
VCF = config.get("vcf", "vcf")
PROJECT_ID = config.get("project_id", "NUMT")


wildcard_constraints:
    chrM_species = "|".join([re.escape(x) for x in chrM_df.index])

def find_mt_seq(wildcards):
    return chrM_df.at[wildcards.species, "fasta"]
def find_mt_repeat2(wildcards):
    return chrM_df.at[wildcards.species, "repeat2_fasta"]
def find_chain(wildcards):
    return chrM_df.at[wildcards.species, "chain"]


rule all:
    input:
        expand("a5-collapse/{PROJECT_ID}_alt_numts_regions.bed", PROJECT_ID = PROJECT_ID)

# filter snv
rule filter_snv:
    input:
        vcf = VCF
    output:
        vcf_filter = "a1-scan/{PROJECT_ID}.biallelic.standard_id.excl_snv.vcf.gz"
    shell:"""
bcftools view {input.vcf} --exclude 'INFO/SVTYPE="SNV"' | bcftools +fill-tags - | bcftools view - --trim-alt-alleles --min-ac 1 --exclude-uncalled -O z -o {output.vcf_filter} && tabix -p vcf {output.vcf_filter}
"""

# get fasta
rule get_alt_seq:
    input:
        vcf_filter = rules.filter_snv.output.vcf_filter
    params:
        alt_blastn_db = "a1-scan/alt_blastn_db/{PROJECT_ID}.alt.blastn_db"
    output:
        alt_fasta = "a1-scan/alt_blastn_db/{PROJECT_ID}.alt.fasta"
    shell:"""
bcftools query {input.vcf_filter} --format '%CHROM\t%POS0\t%ID\t%ALT\n' | awk 'BEGIN{{FS=OFS="\t"}} {{print ">" $3 "\\n" $4}}' > {output.alt_fasta} && samtools faidx {output.alt_fasta}

makeblastdb -in {output.alt_fasta} -dbtype nucl -parse_seqids -out {params.alt_blastn_db}
"""

# blast to scan alt sequences containing NUMTs
rule scan_alt_seq:
    input:
        alt_fasta = rules.get_alt_seq.output.alt_fasta,
        mt_repeat2_seq = find_mt_repeat2
    params:
        alt_blastn_db = "a1-scan/alt_blastn_db/{PROJECT_ID}.alt.blastn_db"
    threads:
        4
    output:
        alt_blastn_result= "a1-scan/alt_blastn_results/{PROJECT_ID}_alt_aln_{species}.blastn.tsv"
    shell:"""
blastn -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -query {input.mt_repeat2_seq} -task blastn -db {params.alt_blastn_db} -out {output.alt_blastn_result} -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qseq sseq" -num_threads {threads} -max_target_seqs 10000
"""

# select variants whose alt seq contain NUMTs
rule select_numts_variants:
    input:
        alt_blastn_results = expand("a1-scan/alt_blastn_results/{PROJECT_ID}_alt_aln_{species}.blastn.tsv", species=chrM_df.index),
        vcf_filter = "a1-scan/{PROJECT_ID}.biallelic.standard_id.excl_snv.vcf.gz"
    output:
        alt_list = "a1-scan/{PROJECT_ID}.alt_numts.list",
        alt_vcf_excl_tr = "a1-scan/{PROJECT_ID}.alt_numts.excl_tr.vcf.gz"
    shell:"""
# for alt sequences
cat {input.alt_blastn_results} | awk 'BEGIN{{FS=OFS="\t"}} !/^#/ && strtonum($3) >= 60 && strtonum($11) <= 1e-3' | cut -f2 | sort -u > {output.alt_list}

bcftools view --include "ID=@{output.alt_list}" -T ^{SCRIPT_DIR}/TR_region.bed {input.vcf_filter} -O z -o {output.alt_vcf_excl_tr} && tabix -p vcf {output.alt_vcf_excl_tr}
"""

# vcf decompose
rule vcf_decompose:
    input:
        alt_vcf_excl_tr = rules.select_numts_variants.output.alt_vcf_excl_tr
    threads:
        24
    params:
        alt_reblastn_db = "a2-redetect/alt_reblastn_db/{PROJECT_ID}.alt.redetect.blastn_db"
    output:
        vcf_decomposed = "a1-scan/{PROJECT_ID}.alt_numts.excl_tr.decomposed.excl_snv.vcf.gz",
        alt_seq_decomposed = "a2-redetect/alt_reblastn_db/{PROJECT_ID}.alt.redetect.fasta"
    shell:"""
vcfwave -t {threads} -I 1000 {input.alt_vcf_excl_tr} | standardize_vcf_id - | bcftools view --exclude 'INFO/SVTYPE="SNV"' - | bcftools +fill-tags - | bcftools view - --trim-alt-alleles --min-ac 1 --exclude-uncalled -O z -o {output.vcf_decomposed}

# get fasta and blastn database
bcftools query {output.vcf_decomposed} --format '%CHROM\t%POS0\t%ID\t%ALT\n' | awk 'BEGIN{{FS=OFS="\t"}} {{print ">" $3 "\\n" $4}}' > {output.alt_seq_decomposed} && samtools faidx {output.alt_seq_decomposed}

makeblastdb -in {output.alt_seq_decomposed} -dbtype nucl -parse_seqids -out {params.alt_reblastn_db}
"""

# redetect numts in variants
rule redetect:
    input:
        alt_seq_decomposed = rules.vcf_decompose.output.alt_seq_decomposed,
        mt_repeat2_seq = find_mt_repeat2
    params:
        alt_reblastn_db = "a2-redetect/alt_reblastn_db/{PROJECT_ID}.alt.redetect.blastn_db"
    threads:
        4
    output:
        alt_reblastn_result = "a2-redetect/alt_reblastn_results/{PROJECT_ID}.alt_aln_{species}.reblastn.tsv"
    shell:"""
blastn -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -query {input.mt_repeat2_seq} -task blastn -db {params.alt_reblastn_db} -out {output.alt_reblastn_result} -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qseq sseq" -num_threads {threads} -max_target_seqs 10000
"""

# reselect
rule reselect_numts_variants:
    input:
        alt_blastn_results = expand("a2-redetect/alt_reblastn_results/{PROJECT_ID}.alt_aln_{species}.reblastn.tsv", species=chrM_df.index),
        vcf_decomposed = rules.vcf_decompose.output.vcf_decomposed
    output:
        alt_list = "a2-redetect/{PROJECT_ID}.alt_numts.redetect.list",
        alt_vcf_excl_tr = "a2-redetect/{PROJECT_ID}.alt_numts.redetect.vcf.gz"
    shell:"""
# for alt sequences
cat {input.alt_blastn_results} | awk 'BEGIN{{FS=OFS="\t"}} !/^#/ && strtonum($3) >= 60 && strtonum($11) <= 1e-3' | cut -f2 | sort -u > {output.alt_list}

bcftools view --include "ID=@{output.alt_list}" -T ^{SCRIPT_DIR}/TR_region.bed {input.vcf_decomposed} -O z -o {output.alt_vcf_excl_tr}
"""

# organize results
rule organize_blast_results:
    input:
        blastn_result = "a2-redetect/alt_reblastn_results/{PROJECT_ID}.alt_aln_{species}.reblastn.tsv",
        mt_seq = find_mt_seq,
        chain = find_chain,
        selected_list = rules.reselect_numts_variants.output.alt_list
    output:
        blastn_result_filtered = temp("a3-liftover/{PROJECT_ID}_alt_aln_{species}.blastn.filtered.tsv"),
        numt_result = temp("a3-liftover/{PROJECT_ID}_alt_aln_{species}.numts_result.tsv"),
        mt_bed = temp("a3-liftover/{PROJECT_ID}_alt_aln_{species}.numts_result.bed"),
        liftover_mt_bed = temp("a3-liftover/{PROJECT_ID}_alt_aln_{species}.mt_liftover.bed"),
        unmapped_bed = temp("a3-liftover/{PROJECT_ID}_alt_aln_{species}.unmapped.bed")
    shell:"""
list_filter -f 2 {input.blastn_result} {input.selected_list} > {output.blastn_result_filtered}

mt_seq_length=$(seqkit fx2tab -n -l {input.mt_seq} | cut -f2)
{SCRIPT_DIR}/numts_dloop_merge.py {output.blastn_result_filtered} $mt_seq_length {wildcards.species} > {output.numt_result}

# convert to 0-based coordinate
awk 'BEGIN{{OFS=FS="\t"}} NR!=1 {{print $4,$5-1,$6,$9,".",$8}}' {output.numt_result} > {output.mt_bed}
liftOver {output.mt_bed} {input.chain} {output.liftover_mt_bed} {output.unmapped_bed}
"""

rule postliftover:
    input:
        numt_result = rules.organize_blast_results.output.numt_result,
        liftover_mt_bed = rules.organize_blast_results.output.liftover_mt_bed
    output:
        liftover_result = "a3-liftover/{PROJECT_ID}_alt_aln_{species}.numts_result.final.tsv"
    run:
        dtypes = {"chrom_adjust": "str", "mt_start_adjust": "int32", "mt_end_adjust": "int32", "id": "str", "score_adjust": "str", "strand_adjust": "str"}
        df_result = pd.read_csv(input.numt_result,sep="\t",header=0)
        df_bed = pd.read_csv(input.liftover_mt_bed,sep="\t",header=None,names=["chrom_adjust", "mt_start_adjust", "mt_end_adjust", "id", "score_adjust", "strand_adjust"],dtype=dtypes)

        # convert human chrM coordinate (repeats3); keep 0-based coordinate
        index1 = (16569 <= df_bed["mt_start_adjust"]) & (df_bed["mt_start_adjust"] < 33138)
        index2 = df_bed["mt_start_adjust"] >= 33138
        df_bed.loc[index1, ["mt_start_adjust", "mt_end_adjust"]] -= 16569
        df_bed.loc[index2, ["mt_start_adjust", "mt_end_adjust"]] -= 33138

        df_combined = pd.merge(df_result,df_bed,on="id",how="left")
        df_combined.to_csv(output.liftover_result, sep="\t", index=False, na_rep='NaN')

rule multi_numt_results_merge:
    input:
        liftover_result = expand("a3-liftover/{PROJECT_ID}_alt_aln_{species}.numts_result.final.tsv", species=chrM_df.index)
    output:
        numt_unmapped = "a4-merge/{PROJECT_ID}_alt_numts_unmapped.tsv",
        numt_bed_igv = "a4-merge/{PROJECT_ID}_alt_numts_result.igv.bed",
        merged_result = temp("a4-merge/{PROJECT_ID}_alt_numts_detection_result.merged.tsv"),
        merged_result_sort = "a4-merge/{PROJECT_ID}_alt_numts_detection_result.merged.sorted.tsv",
        merged_result_sort_bed = "a4-merge/{PROJECT_ID}_alt_numts_detection_result.merged.sorted.bed"
    shell: """
{SCRIPT_DIR}/multi_numt_results_merge.py -s {PROJECT_ID}_alt -w a4-merge {input.liftover_result} > {output.merged_result}
LC_COLLATE=C sort -k1,1 -k2,2n {output.merged_result} > {output.merged_result_sort}
awk 'BEGIN{{OFS=FS="\\t"}}{{print $1,$2,$3,"numt_"NR"."$4":"$5"-"$6,$7,$8}}' {output.merged_result_sort} > {output.merged_result_sort_bed}
"""

# collapse numt associated variants
rule variants_collapse:
    input:
        alt_fasta = rules.vcf_decompose.output.alt_seq_decomposed,
        merged_result_sort = rules.multi_numt_results_merge.output.merged_result_sort
    output:
        alt_numt_bed = "a5-collapse/{PROJECT_ID}_alt_numts_regions.bed",
        alt_numt_seq = "a5-collapse/{PROJECT_ID}_alt_numts_regions.sequence.txt"
    shell:"""
cut -f1-3 {input.merged_result_sort} | bedtools sort | bedtools groupby -g 1 -c 2,3 -o min,max  > {output.alt_numt_bed}

bedtools getfasta -fi {input.alt_fasta} -bed {output.alt_numt_bed} | seqkit fx2tab - | awk 'BEGIN{{OFS=FS="\t"}}{{split($1,a,":"); split(a[1],b,"-"); split(a[2],c,"-"); print b[1],b[2],a[1],c[1],c[2],$2}}' > {output.alt_numt_seq}
"""
