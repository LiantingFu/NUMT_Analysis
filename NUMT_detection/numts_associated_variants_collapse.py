#!/usr/bin/env python

import sys
import argparse
import collections
import Bio.Align
import pandas as pd
import networkx as nx
from itertools import combinations

aligner = Bio.Align.PairwiseAligner()
aligner.mode = "global"

parser = argparse.ArgumentParser("collapse numts associated variants")
parser.add_argument("-v", "--vcf", default="stdin", type=str, help="stdin or a file [stdin]")
parser.add_argument("--sequence", help="variant sequences to compare")
args = parser.parse_args()


def mapping(record1, record2):
    alignment = aligner.align(record1, record2)[0]
    counts = alignment.counts()
    if alignment.score < (counts.gaps + counts.identities + counts.mismatches) * 0.90:
        return False
    return True


def classify_group(group):
    sequences = group["seq"].tolist()
    G = nx.Graph()

    for seq in sequences:
        G.add_node(seq)

    for seq1, seq2 in combinations(sequences, 2):
        if mapping(seq1, seq2):
            G.add_edge(seq1, seq2)

    clusters = list(nx.connected_components(G))

    seq_to_cluster = {}
    for cluster_id, cluster in enumerate(clusters):
        for seq in cluster:
            seq_to_cluster[seq] = cluster_id

    group["cluster_id"] = group.apply(lambda row: f"{row['chrom']}_{row['pos']}_{seq_to_cluster[row['seq']]}", axis=1)

    return group


def merge_sv_list(sv_list):
    if len(sv_list) == 1:
        sys.stderr.write(f"Directly output SV {', '.join([sv_record[2] for sv_record in sv_list])} on {sv_list[0][0]} {sv_list[0][1]}.\n")
        return sv_list[0][:8] + ["GT"] + [sample_detail for sample_detail in sv_list[0][9:]]

    positive_count = collections.Counter()  # count number of positive genotypes for each record
    merged_sample = ["." for _ in range(len(sv_list[0]) - 9)]  # contains genotypes, one element for each sample

    for record_index, sv_record in enumerate(sv_list):
        for sample_index, genotype in enumerate(sv_record[9:]):
            if genotype == "1":
                positive_count[record_index] += 1
            merged_sample[sample_index] = max(merged_sample[sample_index], genotype)

    if positive_count:
        best_record_index = positive_count.most_common(1)[0][0]
    else:
        best_record_index = 0

    sys.stderr.write(f"Merged SVs {', '.join([sv_record[2] for sv_record in sv_list])} on {sv_list[0][0]} {sv_list[0][1]}.\n")

    return sv_list[best_record_index][:8] + ["GT"] + [genotype for genotype in merged_sample]


dtypes = {
    "chrom": "str", "pos": "Int32", "id": "str",
    "start": "Int32", "end": "Int32", "seq": "str"
}

df = pd.read_csv(args.sequence, sep="\t", dtype=dtypes, names=["chrom", "pos", "id", "start", "end", "seq"])

df_group = df.groupby(["chrom", "pos"], as_index=False)

df_result = df_group.apply(classify_group).reset_index(drop=True)

clustered_ids = df_result.groupby("cluster_id")["id"].apply(list).tolist()

vcf_info = {}

for line in sys.stdin if args.vcf in ("stdin", "-") else open(args.vcf, "r"):
    line = line.strip()
    if not line:
        continue

    if line.startswith("#"):
        print(line)
        continue

    record = line.split("\t")
    flag = 0

    vcf_info[record[2]] = record

keys = list(vcf_info.keys())

for cluster in clustered_ids:
    sv_list = []
    for variant_id in cluster:
        sv_list.append(vcf_info[variant_id])
    merged_sv = merge_sv_list(sv_list)
    print("\t".join(merged_sv))
