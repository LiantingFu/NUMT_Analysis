#!/usr/bin/env python

rule merge_results:
    input:
        input_refnumt = "results/REFNUMT-output/4-output/{prefix}.reference-NUMT_results.tsv",
        input_pannumt = "results/PANNUMT-output/6-output/{prefix}.polymorphic_pangenome-NUMT_results.tsv",
        polyref_variant_pair = "results/POLY_REFNUMT-output/2-output/{prefix}.reference-NUMT_associated.tsv",
        pannumt_variant_vcf = "results/PANNUMT-output/6-output/{prefix}.polymorphic_pangenome-NUMT_results.vcf",
        poly_refnumt_variant_vcf = "results/POLY_REFNUMT-output/2-output/{prefix}.polymorphic_reference-NUMT_results.vcf"
    output:
        poly_refnumt = "results/final-output/{prefix}.polymorphic_reference-NUMT.info.tsv",
        fixed_refnumt = "results/final-output/{prefix}.fixed_reference-NUMT.info.tsv",
        poly_pannumt = "results/final-output/{prefix}.polymorphic_pangenome-NUMT.info.tsv",
        all_numt = "results/final-output/{prefix}.NUMT_detection_results.tsv",
        numt_associated_vcf = "results/final-output/{prefix}.NUMT_detection_results.vcf"
    shell:
        """
        csvtk join -H -t -f "1;1" {input.input_refnumt} {input.polyref_variant_pair} | awk 'BEGIN{{OFS=FS="\t"}}{{print $1,"reference","polymoprhic",$11,$3,$4,$5,$6,$7,$8,$9,$10}}' > {output.poly_refnumt}
        
        cut -f1 {input.polyref_variant_pair} | {PIPELINE_DIR}/scripts/list_excluding -f1 {input.input_refnumt} - | awk 'BEGIN{{OFS=FS="\t"}}{{print $1,"reference","fixed",".",$3,$4,$5,$6,$7,$8,$9,$10}}' > {output.fixed_refnumt}

        awk 'BEGIN{{OFS=FS="\t"}}{{print $1,"pangenome","polymorphic",$2,$3,$4,$5,$6,$7,$8,$9,$10}}' {input.input_pannumt} > {output.poly_pannumt}

        cat {output.fixed_refnumt} {output.poly_refnumt} {output.poly_pannumt} > {output.all_numt}

        bcftools concat {input.pannumt_variant_vcf} {input.poly_refnumt_variant_vcf} > {output.numt_associated_vcf}
        """

