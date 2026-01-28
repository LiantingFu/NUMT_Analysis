#!/usr/bin/env python

import pandas as pd

# Remove SNVs and non-informative sites (uncalled/monomorphic)
rule pan_clean_variants:
    input:
        vcf = VCF
    output:
        vcf_filtered = "results/PANNUMT-output/1-screen/{prefix}.noSNV.vcf.gz"
    shell:
        """
        bcftools view {input.vcf} --exclude 'INFO/SVTYPE="SNV"' \
            | bcftools +fill-tags - -- -t AC,AN \
            | bcftools view - \
                --trim-alt-alleles \
                --min-ac 1 \
                --exclude-uncalled \
                -O z -o {output.vcf_filtered}

        tabix -p vcf {output.vcf_filtered}
        """


# Extract ALT alleles as FASTA and build BLAST database
rule pan_build_alt_blastdb:
    input:
        vcf = rules.pan_clean_variants.output.vcf_filtered
    output:
        fasta = "results/PANNUMT-output/1-screen/blastdb/{prefix}-ALT.fasta"
    params:
        db_prefix = "results/PANNUMT-output/1-screen/blastdb/{prefix}-ALT.blastdb"
    shell:
        """
        bcftools query {input.vcf} --format '%CHROM\t%POS0\t%ID\t%ALT\n' \
            | awk 'BEGIN{{FS=OFS="\t"}} {{print ">" $3 "\\n" $4}}' \
            > {output.fasta} && samtools faidx {output.fasta}

        mkdir -p results/PANNUMT-output/1-screen/blastdb

        makeblastdb \
            -in {output.fasta} \
            -dbtype nucl \
            -parse_seqids \
            -out {params.db_prefix} \
        """


# Screen for NUMT candidates
rule pan_screen_numt_candidates:
    input:
        fasta = rules.pan_build_alt_blastdb.output.fasta,
        query_seq = find_mt_repeat2
    params:
        db_prefix = "results/PANNUMT-output/1-screen/blastdb/{prefix}-ALT.blastdb"
    threads:
        4
    output:
        blastn_hits = "results/PANNUMT-output/1-screen/blastn_candidates/{prefix}-ALT.aln.{species}.mt_hits.tsv"
    shell:
        """
        blastn -task blastn -reward 2 -penalty -3 -gapopen 5 -gapextend 2 \
            -query {input.query_seq} \
            -db {params.db_prefix} \
            -out {output.blastn_hits} \
            -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qseq sseq" \
            -num_threads {threads} -max_target_seqs 10000
        """


# Filter BLAST hits
rule pan_filter_numt_vcf:
    input:
        blastn_hits = expand("results/PANNUMT-output/1-screen/blastn_candidates/{{prefix}}-ALT.aln.{species}.mt_hits.tsv", species=mt_manifest_df.index),
        vcf = rules.pan_clean_variants.output.vcf_filtered
    output:
        id_list = "results/PANNUMT-output/1-screen/blastn_candidates/{prefix}-ALT.NUMT_candidates_ids.list",
        vcf_filtered = "results/PANNUMT-output/1-screen/blastn_candidates/{prefix}-ALT.NUMT_candidates.vcf.gz"
    shell:
        """
        cat {input.blastn_hits} \
            | awk 'BEGIN{{FS=OFS="\\t"}} !/^#/ && $3 >= 60 && $11 <= 1e-3' \
            | cut -f2 | sort -u > {output.id_list}

        bcftools view \
            --include "ID=@{output.id_list}" \
            {input.vcf} \
            -O z -o {output.vcf_filtered}
        
        tabix -p vcf {output.vcf_filtered}
        """


# Decompose NUMT-assoicated complex variants and update BLAST database
rule pan_vcf_decompose:
    input:
        vcf = rules.pan_filter_numt_vcf.output.vcf_filtered
    output:
        vcf_decomposed = "results/PANNUMT-output/2-refine/{prefix}-ALT.NUMT_candidates.decomposed.noSNV.vcf.gz",
        fasta = "results/PANNUMT-output/2-refine/blastdb/{prefix}-ALT.NUMT_candidates.decomposed.fasta",
    params:
        db_prefix = "results/PANNUMT-output/2-refine/blastdb/{prefix}-ALT.NUMT_candidates.decomposed.blastdb"
    threads:
        24
    shell:
        """
        vcfwave -t {threads} -I 1000 {input.vcf} \
            | {PIPELINE_DIR}/scripts/standardize_vcf_id - \
            | bcftools view --exclude 'INFO/SVTYPE="SNV"' - \
            | bcftools +fill-tags - \
            | bcftools view - \
                --trim-alt-alleles \
                --min-ac 1 \
                --exclude-uncalled \
                -O z -o {output.vcf_decomposed}

        bcftools query {output.vcf_decomposed} --format '%CHROM\t%POS0\t%ID\t%ALT\n' \
            | awk 'BEGIN{{FS=OFS="\\t"}} {{print ">" $3 "\\n" $4}}' \
            > {output.fasta}

        mkdir -p results/PANNUMT-output/2-refine/blastdb

        makeblastdb \
            -in {output.fasta} \
            -dbtype nucl \
            -parse_seqids \
            -out {params.db_prefix}
        """


# Re-scan decomposed NUMT-assoicated variants to refine NUMT boundaries
rule pan_refine_numt_candidates:
    input:
        fasta = rules.pan_vcf_decompose.output.fasta,
        query_seq = find_mt_repeat2
    params:
        db_prefix = "results/PANNUMT-output/2-refine/blastdb/{prefix}-ALT.NUMT_candidates.decomposed.blastdb"
    threads:
        4
    output:
        blastn_hits = "results/PANNUMT-output/2-refine/blastn_candidates_refine/{prefix}-ALT.aln.{species}.mt_hits.refined.tsv"
    shell:
        """
        blastn -task blastn -reward 2 -penalty -3 -gapopen 5 -gapextend 2 \
            -query {input.query_seq} \
            -db {params.db_prefix} \
            -out {output.blastn_hits} \
            -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qseq sseq" \
            -num_threads {threads} -max_target_seqs 10000
        """


# Filter decomposed NUMT-assoicated variants based on refined BLAST hits
rule pan_filter_refined_numts:
    input:
        blastn_hits = expand("results/PANNUMT-output/2-refine/blastn_candidates_refine/{{prefix}}-ALT.aln.{species}.mt_hits.refined.tsv", species=mt_manifest_df.index),
        vcf_decomposed = rules.pan_vcf_decompose.output.vcf_decomposed
    output:
        id_list = "results/PANNUMT-output/2-refine/blastn_candidates_refine/{prefix}-ALT.NUMT_refined_ids.list",
        vcf_final = "results/PANNUMT-output/2-refine/blastn_candidates_refine/{prefix}-ALT.NUMT_refined.vcf.gz"
    shell:
        """
        cat {input.blastn_hits} \
            | awk 'BEGIN{{FS=OFS="\\t"}} !/^#/ && $3 >= 60 && $11 <= 1e-3' \
            | cut -f2 | sort -u > {output.id_list}

        bcftools view \
            --include "ID=@{output.id_list}" \
            {input.vcf_decomposed} \
            -O z -o {output.vcf_final}

        tabix -p vcf {output.vcf_final}
        """


# Filter, Merge D-loop regions, and LiftOver coordinates
rule pan_liftover_blast_results:
    input:
        blastn_hits = rules.pan_refine_numt_candidates.output.blastn_hits,
        valid_ids = rules.pan_filter_refined_numts.output.id_list,
        mt_seq = find_mt_seq,
        chain = find_chain
    output:
        hits_filtered = "results/PANNUMT-output/3-liftover/{prefix}-ALT.{species}.mt_hits.filtered.tsv",
        merged_results = temp("results/PANNUMT-output/3-liftover/{prefix}-ALT.{species}.Dloop-merged.tsv"),
        mt_bed = temp("results/PANNUMT-output/3-liftover/{prefix}-ALT.{species}.mt_coords.bed"),
        lifted_bed = temp("results/PANNUMT-output/3-liftover/{prefix}-ALT.{species}.lifted.bed"),
        unmapped_bed = "results/PANNUMT-output/3-liftover/{prefix}-ALT.{species}.unMapped.bed"
    shell:
        """
        {PIPELINE_DIR}/scripts/list_filter -f 2 {input.blastn_hits} {input.valid_ids} > {output.hits_filtered}

        mt_len=$(seqkit fx2tab -n -l {input.mt_seq} | cut -f2)

        awk 'BEGIN{{OFS=FS="\t"}} $1!~/^#/ {{
        if($9 > $10) {{
            print $2, $10, $9, $1, $7, $8, "-", $3
        }} else {{
            print $2, $9, $10, $1, $7, $8, "+", $3
        }}
    }}' {output.hits_filtered} \
        | {PIPELINE_DIR}/scripts/NUMT_Dloop_merge.py - \
            -l $mt_len \
            -s {wildcards.species} \
            > {output.merged_results}

        awk 'BEGIN{{OFS=FS="\\t"}} NR!=1 {{print $5, $6-1, $7, $1, ".", $8}}' {output.merged_results} > {output.mt_bed}

        liftOver \
            {output.mt_bed} \
            {input.chain} \
            {output.lifted_bed} \
            {output.unmapped_bed}
        """


# Normalize coordinates to reference mtDNA
rule pan_merge_liftover_results:
    input:
        merged_results = rules.pan_liftover_blast_results.output.merged_results,
        lifted_bed = rules.pan_liftover_blast_results.output.lifted_bed
    output:
        lifted_results  = "results/PANNUMT-output/3-liftover/{prefix}-ALT.{species}.NUMT-all.tsv"
    shell:
        """
        {PIPELINE_DIR}/scripts/NUMT_post-liftOver.py -r {input.merged_results} -b {input.lifted_bed} > {output.lifted_results}
        """


# Merge multiple query callsets
rule pan_merge_query_results:
    input:
        liftover_results = expand("results/PANNUMT-output/3-liftover/{{prefix}}-ALT.{species}.NUMT-all.tsv", species=mt_manifest_df.index)
    output:
        unmapped = "results/PANNUMT-output/4-merge/{prefix}-ALT.NUMT.unMapped.tsv",
        igv_bed  = "results/PANNUMT-output/4-merge/{prefix}-ALT.NUMT.IGV.bed",
        merged_tsv = temp("results/PANNUMT-output/4-merge/{prefix}-ALT.NUMT.merged.tsv"),
        numt_tsv = "results/PANNUMT-output/4-merge/{prefix}-ALT.NUMT-all.sorted.tsv",
        numt_bed  = "results/PANNUMT-output/4-merge/{prefix}-ALT.NUMT-all.bed"
    shell:
        """
        {PIPELINE_DIR}/scripts/NUMT_query-callsets_merge.py \
            -p {wildcards.prefix}-ALT.NUMT \
            -w results/PANNUMT-output/4-merge \
            {input.liftover_results} > {output.merged_tsv}

        LC_COLLATE=C sort -k1,1 -k2,2n {output.merged_tsv} > {output.numt_tsv}

        awk 'BEGIN{{OFS=FS="\\t"}} {{print $1, $2, $3, "NUMT_"NR"."$4":"$5"-"$6, $7, $8}}' \
            {output.numt_tsv} > {output.numt_bed}
        """


# collapse NUMT records
rule pan_collapse_records:
    input:
        alt_fasta = rules.pan_vcf_decompose.output.fasta,
        sorted_hits = rules.pan_merge_query_results.output.numt_tsv
    output:
        regions_bed = "results/PANNUMT-output/5-collapse/{prefix}-ALT.NUMT.collapsed.bed",
        regions_seq = "results/PANNUMT-output/5-collapse/{prefix}-ALT.NUMT.collapsed.txt"
    shell:
        """
        cut -f1-3 {input.sorted_hits} \
            | bedtools sort -i - \
            | bedtools groupby -g 1 -c 2,3 -o min,max \
            > {output.regions_bed}

        bedtools getfasta -fi {input.alt_fasta} -bed {output.regions_bed} \
            | seqkit fx2tab - \
            | awk 'BEGIN{{OFS=FS="\\t"}} {{
                split($1, a, ":");
                split(a[1], b, "-");
                split(a[2], c, "-");
                print b[1], b[2], a[1], c[1], c[2], $2
            }}' > {output.regions_seq}
        """


# output final NUMT results
rule pan_output_final_results:
    input:
        vcf_decomposed = rules.pan_vcf_decompose.output.vcf_decomposed,
        regions_seq = rules.pan_collapse_records.output.regions_seq,
        numt_tsv = rules.pan_merge_query_results.output.numt_tsv
    output:
        results_vcf = "results/PANNUMT-output/6-output/{prefix}.polymorphic_pangenome-NUMT_results.vcf",
        id_list = temp("results/PANNUMT-output/6-output/{prefix}.polymorphic_pangenome-NUMT_results.id.list"),
        results_tsv = "results/PANNUMT-output/6-output/{prefix}.polymorphic_pangenome-NUMT_results.tsv"
    shell:
        """
        zcat {input.vcf_decomposed} \
            | {PIPELINE_DIR}/scripts/NUMT_collapse.py \
                -v - --sequence {input.regions_seq} | bcftools sort - > {output.results_vcf}

        bcftools query -f "%ID\n" {output.results_vcf} > {output.id_list}

        list_filter -f1 {input.numt_tsv} {output.id_list} | awk 'BEGIN{{OFS=FS="\t"}}{{split($1,a,"-"); print "PANNUMT_"NR, $1, a[1], a[2], $2, $3, $4, $5, $6, $8}}' > {output.results_tsv}
        """
