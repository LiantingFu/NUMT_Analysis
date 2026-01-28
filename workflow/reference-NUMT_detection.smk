#!/usr/bin/env python

import pandas as pd

# Make BLAST database
rule ref_makedb:
    input:
        fasta = ref_fasta
    params:
        db_prefix = "results/REFNUMT-output/1-mt_hits/blastdb/{prefix}-REF.blastdb"
    output:
        flag = temp("results/REFNUMT-output/1-mt_hits/{prefix}-REF.blastdb.flag")
    shell:
        """
        mkdir -p results/REFNUMT-output/1-mt_hits/blastdb

        makeblastdb \
            -in {input.fasta} \
            -dbtype nucl \
            -parse_seqids \
            -out {params.db_prefix} \
        
        touch {output.flag}
        """


# Detect NUMT candidates
rule ref_detect_numt_candidates:
    input:
        flag = rules.ref_makedb.output.flag,
        query_seq = find_mt_repeat2
    threads:
        4
    params:
        db_prefix = "results/REFNUMT-output/1-mt_hits/blastdb/{prefix}-REF.blastdb"
    output:
        blastn_hits = "results/REFNUMT-output/1-mt_hits/blastn_candidates/{prefix}-REF.aln.{species}.mt_hits.tsv"
    shell:
        """
        blastn -task blastn -reward 2 -penalty -3 -gapopen 5 -gapextend 2 \
            -query {input.query_seq} \
            -db {params.db_prefix} \
            -out {output.blastn_hits} \
            -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qseq sseq" \
            -num_threads {threads} -max_target_seqs 10000
        """


# Filter, Merge D-loop regions, and LiftOver coordinates
rule ref_liftover_blast_results:
    input:
        blastn_hits = rules.ref_detect_numt_candidates.output.blastn_hits,
        mt_seq = find_mt_seq,
        chain = find_chain
    output:
        hits_filtered = "results/REFNUMT-output/2-liftover/{prefix}-REF.{species}.mt_hits.filtered.tsv",
        merged_results = temp("results/REFNUMT-output/2-liftover/{prefix}-REF.{species}.Dloop-merged.tsv"),
        mt_bed = temp("results/REFNUMT-output/2-liftover/{prefix}-REF.{species}.mt_coords.bed"),
        lifted_bed = temp("results/REFNUMT-output/2-liftover/{prefix}.{species}.lifted.bed"),
        unmapped_bed = "results/REFNUMT-output/2-liftover/{prefix}.{species}.unMapped.bed"
    shell:
        """
        awk 'BEGIN{{FS=OFS="\\t"}} !/^#/ && $3 >= 60 && $11 <= 1e-3' {input.blastn_hits} > {output.hits_filtered}

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
rule ref_merge_liftover_results:
    input:
        merged_results = rules.ref_liftover_blast_results.output.merged_results,
        lifted_bed = rules.ref_liftover_blast_results.output.lifted_bed
    output:
        lifted_results  = "results/REFNUMT-output/2-liftover/{prefix}-REF.{species}.NUMT-all.tsv"
    shell:
        """
        {PIPELINE_DIR}/scripts/NUMT_post-liftOver.py -r {input.merged_results} -b {input.lifted_bed} > {output.lifted_results}
        """


# Merge multiple query callsets
rule ref_merge_query_results:
    input:
        liftover_results = expand("results/REFNUMT-output/2-liftover/{{prefix}}-REF.{species}.NUMT-all.tsv", species=mt_manifest_df.index)
    output:
        unmapped = "results/REFNUMT-output/3-merge/{prefix}-REF.NUMT.unMapped.tsv",
        igv_bed = "results/REFNUMT-output/3-merge/{prefix}-REF.NUMT.IGV.bed",
        merged_tsv = temp("results/REFNUMT-output/3-merge/{prefix}-REF.NUMT.merged.tsv"),
        numt_tsv = "results/REFNUMT-output/3-merge/{prefix}-REF.NUMT-all.sorted.tsv",
        numt_bed  = "results/REFNUMT-output/3-merge/{prefix}-REF.NUMT-all.bed",
        results_tsv = "results/REFNUMT-output/4-output/{prefix}.reference-NUMT_results.tsv"
    shell:
        """
        {PIPELINE_DIR}/scripts/NUMT_query-callsets_merge.py \
            -p {wildcards.prefix}-REF.NUMT \
            -w results/REFNUMT-output/3-merge \
            {input.liftover_results} > {output.merged_tsv}

        LC_COLLATE=C sort -k1,1 -k2,2n {output.merged_tsv} > {output.numt_tsv}

        awk 'BEGIN{{OFS=FS="\\t"}} {{print $1, $2, $3, "NUMT_"NR"."$4":"$5"-"$6, $7, $8}}' \
            {output.numt_tsv} > {output.numt_bed}

        awk 'BEGIN{{OFS=FS="\t"}}{{print "REFNUMT_"NR, $1"-"$2"-"$3, $1, int(($2+$3)/2), $2, $3, $4, $5, $6, $8}}' {output.numt_tsv} > {output.results_tsv}
        """
