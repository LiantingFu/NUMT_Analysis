#!/usr/bin/env python

import re
import pandas as pd


configfile: "config/PanGenie_genotyping.json"

manifest_df = pd.read_csv(config.get("samples", "config/samples.tsv"), header=0, index_col="sample", sep="\t")

multi_vcf = config.get("multi_vcf")
bi_vcf = config.get("bi_vcf")
reference = config.get("reference")
cram_reference = config.get("cram_reference")

PANGENIE_INDEX = config.get("PanGenie-index", "PanGenie-index")
PANGENIE = config.get("PanGenie", "PanGenie")
TABIX = config.get("tabix", "tabix")


wildcard_constraints:
    sample = "|".join([re.escape(x) for x in manifest_df.index])


rule all:
    input:
        multi_vcf = expand("results/{sample}/{sample}_genotyping.sorted.vcf.gz", sample=manifest_df.index),
        bi_vcf = expand("results/{sample}/{sample}_genotyping.biallelic.sorted.vcf.gz", sample=manifest_df.index)


rule PanGenie_index:
    input:
        reference = reference,
        multiallelic_panel_vcf = multi_vcf
    output:
        index_done = "PanGenie_index/PanGenie_index_work.done"
    threads:
        16
    shell:"""
{PANGENIE_INDEX} \
    -r {input.reference} \
    -v {input.multiallelic_panel_vcf} \
    -t {threads} \
    -o PanGenie_index/PanGenie
touch {output.index_done}
"""


rule uncompress_reads:
    input:
        cram = lambda wildcards: manifest_df.at[wildcards.sample, "cram"],
        cram_reference = cram_reference
    output:
        uncompressed_reads = temp("temp/{sample}.fastq")
    group: "genotyping"
    threads: 4
    shell:"""
samtools fastq -@ {threads} --reference {input.cram_reference} {input.cram} -o {output.uncompressed_reads}
"""


rule PanGenie_genotyping:
    input:
        index_done = rules.PanGenie_index.output.index_done,
        uncompressed_reads = rules.uncompress_reads.output.uncompressed_reads
    output:
        vcf = temp("results/{sample}/{sample}_genotyping.vcf"),
        histo = "results/{sample}/{sample}_histogram.histo"
    group: "genotyping"
    threads:
        16
    shell:"""
{PANGENIE} \
    -f PanGenie_index/PanGenie \
    -i {input.uncompressed_reads} \
    -s {wildcards.sample} \
    -o results/{wildcards.sample}/{wildcards.sample} \
    -j {threads} \
    -t {threads}
"""


rule bgzip_vcf:
    input:
        vcf = rules.PanGenie_genotyping.output.vcf,
        reference = reference
    output:
        vcf_compress = temp("results/{sample}/{sample}_genotyping.vcf.gz"),
        vcf_sort = "results/{sample}/{sample}_genotyping.sorted.vcf.gz",
        vcf_sort_index = "results/{sample}/{sample}_genotyping.sorted.vcf.gz.tbi"
    shell:"""
bgzip -c {input.vcf} > {output.vcf_compress}
bcftools reheader {output.vcf_compress} -f {input.reference}.fai | bcftools sort - -O z -o {output.vcf_sort}
{TABIX} -p vcf {output.vcf_sort}
"""


rule to_biallelic:
    input:
        reference = reference,
        vcf = rules.PanGenie_genotyping.output.vcf,
        biallelic_panel_vcf = bi_vcf
    output:
        biallelic_vcf = temp("results/{sample}/{sample}_genotyping.biallelic.vcf.gz"),
        biallelic_vcf_sort = "results/{sample}/{sample}_genotyping.biallelic.sorted.vcf.gz",
        biallelic_vcf_index = "results/{sample}/{sample}_genotyping.biallelic.sorted.vcf.gz.tbi"
    threads:
        12
    shell:"""
cat {input.vcf} \
| python scripts/convert-to-biallelic.py {input.biallelic_panel_vcf} \
| bgzip > {output.biallelic_vcf}

bcftools reheader {output.biallelic_vcf} -f {input.reference}.fai | bcftools sort - -O z -o {output.biallelic_vcf_sort}
{TABIX} -p vcf {output.biallelic_vcf_sort}
"""
