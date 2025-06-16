#!/usr/bin/env python

import pandas as pd
import sys


manifest_df = pd.read_csv("polymorphic_numts.over1kb.info.tsv", header=0, index_col="row_id", sep="\t")
SOFTWARE_DIR = "software"


wildcard_constraints:
    row_id = "|".join([re.escape(x) for x in manifest_df.index])


rule all:
    input:
        [f"results/{row['numt_id']}/{row['hap']}/{row.name}.methyl_recalc_q10.combined.bed.gz" for _, row in manifest_df.iterrows()],
        [f"results/{row['numt_id']}/{row['hap']}/{row.name}.200bp.hifi_methyl_freq.q10.bed" for _, row in manifest_df.iterrows()]



rule get_ref:
    input:
        nu_asm = lambda wildcards: manifest_df.at[wildcards.row_id, "nu_asm"],
        mt_asm = lambda wildcards: manifest_df.at[wildcards.row_id, "mt_asm"]
    output:
        nu_fasta = temp("results/{numt_id}/{hap}/{row_id}.nu.fasta"),
        mt_fasta = temp("results/{numt_id}/{hap}/{row_id}.mt.fasta"),
        mt_repeat2_fasta = temp("results/{numt_id}/{hap}/{row_id}.mt_repeat2.fasta"),
        ref_fasta = "results/{numt_id}/{hap}/{row_id}.ref.fasta",
        ref_mmi = "results/{numt_id}/{hap}/{row_id}.ref.mmi"
    params:
        nu_chrom = lambda wildcards: manifest_df.at[wildcards.row_id, "nu_chrom"],
        mt_chrom = lambda wildcards: manifest_df.at[wildcards.row_id, "mt_chrom"]
    shell:"""
    # chromsome + chrM
    samtools faidx {input.nu_asm} {params.nu_chrom} > {output.nu_fasta}
    samtools faidx {input.mt_asm} {params.mt_chrom} > {output.mt_fasta}
    
    (echo ">chrM_repeat2" && grep -v ">" {output.mt_fasta} | tr -d '\\n' && grep -v ">" {output.mt_fasta} | tr -d '\\n') > {output.mt_repeat2_fasta}
    
    cat {output.nu_fasta} {output.mt_repeat2_fasta} > {output.ref_fasta}

    # pbmm2 index
    {SOFTWARE_DIR}/pbmm2 index {output.ref_fasta} {output.ref_mmi} --preset CCS
    """


rule align_bam:
    input:
        ref_mmi = rules.get_ref.output.ref_mmi,
        numt_bam = lambda wildcards: manifest_df.at[wildcards.row_id, "bam"]
    output:
        numt_bam_align = "results/{numt_id}/{hap}/{row_id}.pbmm2_align.bam"
    shell:"""
    {SOFTWARE_DIR}/pbmm2 align --preset CCS --sort {input.ref_mmi} {input.numt_bam} {output.numt_bam_align}
    """


rule calc_methyl:
    input:
        numt_bam_align = rules.align_bam.output.numt_bam_align
    threads:
        8
    output:
        methyl_bed = "results/{numt_id}/{hap}/{row_id}.methyl_recalc_q10.combined.bed.gz"
    shell:"""
    {SOFTWARE_DIR}/pb-CpG-tools/bin/aligned_bam_to_cpg_scores \
    --bam {input.numt_bam_align} \
    --output-prefix results/{wildcards.numt_id}/{wildcards.hap}/{wildcards.row_id}.methyl_recalc_q10 \
    --min-mapq 10 \
    --threads {threads}
    """


rule calc_window_methyl:
    input:
        methyl_bed = rules.calc_methyl.output.methyl_bed
    output:
        window_bed = "results/{numt_id}/{hap}/{row_id}.200bp.hifi_methyl_freq.q10.bed"
    params:
        region= lambda wildcards: "\t".join([str(manifest_df.at[wildcards.row_id,"nu_chrom"]), str(manifest_df.at[wildcards.row_id, "start"]), str(manifest_df.at[wildcards.row_id, "end"]), str(manifest_df.at[wildcards.row_id, "numt_id"])]),
        hap = lambda wildcards: manifest_df.at[wildcards.row_id, "hap"]
    shell:"""
    echo -e "{params.region}" | bedtools makewindows -b /dev/stdin -w 200 -i srcwinnum | sort -k1,1V -k2,2n | bedtools map -a /dev/stdin -b {input.methyl_bed} -c 4 -o mean | awk -v id="{params.hap}" 'BEGIN{{OFS=FS="\\t"}} {{print $0, id}}' > {output.window_bed}
    """
