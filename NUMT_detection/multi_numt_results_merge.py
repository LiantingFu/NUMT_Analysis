#!/usr/bin/env python

import sys
import argparse
import pandas as pd
from collections import defaultdict
from collections import Counter
import pybedtools

parser = argparse.ArgumentParser(description="merge numt calling results from different chrM reference sequence")
parser.add_argument("input_numt_results", nargs="+", help="numt results")
parser.add_argument("--sample", "-s", required=True, help="sample name")
parser.add_argument("--work_dir", "-w", required=True, help="work directory")
args = parser.parse_args()

sample = args.sample
input_numt_info = args.input_numt_results
work_dir = args.work_dir

dtypes = {
    "nu_chrom": "str", "nu_start": "Int32", "nu_end": "Int32",
    "mt_chrom": "str", "mt_start": "Int32", "mt_end": "Int32",
    "identity": "float32", "strand": "str", "id": "str",
    "chrom_adjust": "str", "mt_start_adjust": "Int32", "mt_end_adjust": "Int32", "score_adjust": "str", "strand_adjust": "str"
}

input_dfs = [pd.read_csv(file, sep="\t", header=0, dtype=dtypes) for file in input_numt_info]
df_info_merge = pd.concat(input_dfs)

# discard and output unmapped numt records
df_info_merge_unmapped = df_info_merge[df_info_merge.isna().any(axis=1)]
df_info_merge_unmapped.to_csv(f"{work_dir}/{sample}_numts_unmapped.tsv", sep="\t", index=False)
df_info_merge = df_info_merge.dropna()

# output bed file for igv
df_info_merge_igv_bed_columns = df_info_merge[["nu_chrom", "nu_start", "nu_end", "id", "mt_start_adjust", "mt_end_adjust", "score_adjust", "strand_adjust"]]
df_info_merge_igv_bed = df_info_merge_igv_bed_columns.apply(
    lambda row: f"{row['nu_chrom']}\t{row['nu_start']}\t{row['nu_end']}\t{row['id']}.chrM:{row['mt_start_adjust']}-{row['mt_end_adjust']}\t{row['score_adjust']}\t{row['strand_adjust']}", axis=1
)
df_info_merge_igv_bed.to_csv(f"{work_dir}/{sample}_numts_result.igv.bed", index=False, header=False)

# multiple numt results intersect: r>=0.8 or R>=0.8
df_info_merge_bed_format = df_info_merge[["nu_chrom", "nu_start", "nu_end", "id", "score_adjust", "strand_adjust"]].to_csv(sep="\t", index=False, header=False)

df_info_merge_bed = pybedtools.BedTool(df_info_merge_bed_format, from_string=True)

intersect_result = df_info_merge_bed.intersect(df_info_merge_bed, wo=True)

# cluster numts according coordinte on nuclear genome
id_to_group = defaultdict(set)

df_info_merge.set_index("id", inplace=True)
for index, row in df_info_merge.iterrows():
    id_to_group[index].add(index)

for record in intersect_result:
    (chrom1, start1, end1, id1, score1, strand1,
     chrom2, start2, end2, id2, score2, strand2, overlap_length) = record[:13]
    start1, end1, start2, end2, overlap_length = int(start1), int(end1), int(start2), int(end2), int(overlap_length)
    rate1 = overlap_length / (end1 - start1)
    rate2 = overlap_length / (end2 - start2)
    if id1 == id2 or (rate1 < 0.8 and rate2 < 0.8):
        continue

    id1_group = id_to_group.get(id1)
    id2_group = id_to_group.get(id2)

    id1_group.update(id2_group)
    for id in id1_group:
        id_to_group[id] = id1_group

groups = set(map(frozenset, id_to_group.values()))

# get final coordinates: min start and max end for each cluster
for group in groups:
    block_chrom, block_strand = [], []
    final_nu_start, final_nu_end = None, None
    final_mt_start, final_mt_end = None, None

    for numt_id in group:
        nu_start = df_info_merge.loc[numt_id]["nu_start"]
        nu_end = df_info_merge.loc[numt_id]["nu_end"]
        strand = df_info_merge.loc[numt_id]["strand_adjust"]

        if final_nu_start is None or nu_start < final_nu_start:
            final_nu_start = nu_start
            if strand == "+":
                final_mt_start = df_info_merge.loc[numt_id]["mt_start_adjust"]
            else:
                final_mt_end = df_info_merge.loc[numt_id]["mt_end_adjust"]

        if final_nu_end is None or nu_end > final_nu_end:
            final_nu_end = nu_end
            if strand == "+":
                final_mt_end = df_info_merge.loc[numt_id]["mt_end_adjust"]
            else:
                final_mt_start = df_info_merge.loc[numt_id]["mt_start_adjust"]

        block_chrom.append(df_info_merge.loc[numt_id]["nu_chrom"])
        block_strand.append(df_info_merge.loc[numt_id]["strand_adjust"])

    strand_main = Counter(block_strand).most_common(1)[0][0]

    if final_mt_end < final_mt_start:
        final_mt_end = final_mt_end + 16569

    print(block_chrom[0], final_nu_start, final_nu_end, "chrM", final_mt_start, final_mt_end, ".", strand_main, len(group), ",".join(group), sep="\t")
