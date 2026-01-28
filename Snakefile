#!usr/bin/env python3

import os
import pandas as pd

configfile: "config/config.json"

mt_manifest_df = pd.read_csv(config.get("chrM_info", "config/chrM_info.tsv"), header=0, index_col="species", sep="\t")

ref_fasta = config.get("ref", "ref")
VCF = config.get("pangenome_vcf", "pangenome_vcf")
PREFIX = config.get("prefix", "pan-NUMT")
PIPELINE_DIR = os.path.dirname(workflow.snakefile)

wildcard_constraints:
    mt_species = "|".join([re.escape(x) for x in mt_manifest_df.index])

def find_mt_seq(wildcards):
    return mt_manifest_df.at[wildcards.species, "fasta"]
def find_mt_repeat2(wildcards):
    return mt_manifest_df.at[wildcards.species, "repeat2_fasta"]
def find_chain(wildcards):
    return mt_manifest_df.at[wildcards.species, "chain"]


include: "workflow/reference-NUMT_detection.smk"
include: "workflow/pangenome-NUMT_detection.smk"
include: "workflow/polymorphic_reference-NUMT_detection.smk"
include: "workflow/merge_NUMT_results.smk"

rule all:
    input:
        expand("results/final-output/{prefix}.NUMT_detection_results.vcf", prefix=PREFIX)
