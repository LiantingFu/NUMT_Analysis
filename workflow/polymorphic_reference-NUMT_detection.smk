#!/usr/bin/env python

# Screen for reference-NUMT_associated variants
rule polyref_screen_variants:
    input:
        vcf = VCF,
        results_tsv = rules.ref_merge_query_results.output.results_tsv
    output:
        id_lst = "results/POLY_REFNUMT-output/1-screen/{prefix}.reference-NUMT.intersect.id.list",
        vcf_filtered = "results/POLY_REFNUMT-output/1-screen/{prefix}.reference-NUMT.intersect.vcf.gz",
    shell:
        """
        awk 'BEGIN{{OFS=FS="\t"}}{{print $3,$5,$6,$1}}' {input.results_tsv} | bedtools intersect -a {input.vcf} -b /dev/stdin -wo | cut -f3 | sort -u > {output.id_lst}

        bcftools view -i "ID=@{output.id_lst}" {input.vcf} \
            -O z -o {output.vcf_filtered}

        tabix -p vcf {output.vcf_filtered}
        """


# Decompose ref-NUMT_assoicated complex variants
rule polyref_vcf_decompose:
    input:
        vcf = rules.polyref_screen_variants.output.vcf_filtered
    output:
        vcf_decomposed = "results/POLY_REFNUMT-output/1-screen/{prefix}.reference-NUMT.intersect.decomposed.vcf.gz"
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
        """


# Define polymorphic reference NUMTs: length left < 30bp || deletion rate >= 50%: completely_deletion
rule polyref_refine_variants:
    input:
        vcf = rules.polyref_screen_variants.output.vcf_filtered,
        results_tsv = rules.ref_merge_query_results.output.results_tsv,
        vcf_decomposed = rules.polyref_vcf_decompose.output.vcf_decomposed
    output:
        intersect_tsv = "results/POLY_REFNUMT-output/1-screen/{prefix}.reference-NUMT.intersect.decomposed.intersect.tsv",
        poly_refnumt_variant_pairs = "results/POLY_REFNUMT-output/2-output/{prefix}.reference-NUMT_associated.tsv",
        poly_refnumt_variant_vcf = "results/POLY_REFNUMT-output/2-output/{prefix}.polymorphic_reference-NUMT_results.vcf",
        poly_refnumt_variant_id = "results/POLY_REFNUMT-output/2-output/{prefix}.polymorphic_reference-NUMT_results.id.list"
    shell:
        """
        awk 'BEGIN{{OFS=FS="\t"}}{{print $3,$5,$6,$1}}' {input.results_tsv} \
            | bedtools intersect -a /dev/stdin -b {input.vcf} -wo \
            | awk 'BEGIN{{OFS=FS="\t"}}{{print $1,$2,$3,$4,$7,($3-$2-$NF),($3-$2-$NF)/($3-$2)}}' \
            > {output.intersect_tsv}
        
        awk 'BEGIN{{OFS=FS="\t"}} $6<30 || $7<0.5 {{print $4,$5}}' {output.intersect_tsv} > {output.poly_refnumt_variant_pairs}

        cut -f2 {output.poly_refnumt_variant_pairs} > {output.poly_refnumt_variant_id}

        bcftools view -i "ID=@{output.poly_refnumt_variant_id}" {input.vcf_decomposed} -O z -o {output.poly_refnumt_variant_vcf}
        """
