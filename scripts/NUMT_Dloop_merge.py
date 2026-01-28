#!/usr/bin/env python

import sys
import pandas as pd
import portion as P
import argparse


parser = argparse.ArgumentParser(description="merge NUMTs split in D-loop region")
parser.add_argument("input_hits", type=argparse.FileType("r"), default=sys.stdin, help="hits <nu_chrom, nu_start, nu_end, mt_chrom, mt_start, mt_end, strand, identity>")
parser.add_argument("--mtDNA_len", "-l", type=int, help="Length of mtDNA sequence")
parser.add_argument("--species", "-s", default="default", help="species of mtDNA")

args = parser.parse_args()

mtDNA_len = args.mtDNA_len

# record1: spanning; record2: front
def merge_dloop(record1, record2):
    if record1["nu_chrom"] != record2["nu_chrom"]:
        return False
    nu_region1 = P.closed(record1["nu_start"], record1["nu_end"])
    nu_region2 = P.closed(record2["nu_start"], record2["nu_end"])
    nu_intersect_region = nu_region1 & nu_region2
    if nu_intersect_region.empty:
        flag = False
        return flag

    nu_intersect_rate = (nu_intersect_region.upper - nu_intersect_region.lower) / (nu_region2.upper - nu_region2.lower)

    mt_region1 = P.closed(1, record1["mt_end"] - mtDNA_len)
    mt_region2 = P.closed(record2["mt_start"], record2["mt_end"])
    mt_intersect_region = mt_region1 & mt_region2
    if mt_intersect_region.empty:
        flag = False
        return flag

    mt_intersect_rate = (mt_intersect_region.upper - mt_intersect_region.lower) / (mt_region2.upper - mt_region2.lower)

    if nu_intersect_rate >= 0.95 and mt_intersect_rate >= 0.95 and record1["strand"] == record2["strand"]:
        flag = True
    else:
        flag = False
    return flag


dtypes = {
    "nu_chrom": "str", "nu_start": "int32", "nu_end": "int32",
    "mt_chrom": "str", "mt_start": "int32", "mt_end": "int32",
    "strand": "str", "identity": "float32"
}

df_input = pd.read_csv(args.input_hits, sep="\t", names=["nu_chrom", "nu_start", "nu_end", "mt_chrom", "mt_start", "mt_end", "strand", "identity"], dtype=dtypes, comment="#")

df_input.sort_values(by=["nu_chrom", "nu_start"], inplace=True)

df_front = df_input[df_input["mt_end"] < mtDNA_len].copy()
df_spanning = df_input[(df_input["mt_start"] <= mtDNA_len) & (df_input["mt_end"] >= mtDNA_len)].copy()

df_front["discard"] = False

for index2, row2 in df_front.iterrows():
    for index, row in df_spanning.iterrows():
        df_front.at[index2, "discard"] = merge_dloop(row, row2)
        if df_front.at[index2, "discard"]:
            break

df_front_not_discarded = df_front[~df_front["discard"]].drop(columns=["discard"])
df_combined = pd.concat([df_spanning, df_front_not_discarded], ignore_index=True)
df_combined.sort_values(by=["nu_chrom", "nu_start"], inplace=True)
df_combined.reset_index(drop=True, inplace=True)
df_combined["numt_id"] = "numt_" + (df_combined.index + 1).astype(str) + "." + args.species


cols_order = ["numt_id", "nu_chrom", "nu_start", "nu_end", "mt_chrom", "mt_start", "mt_end", "strand", "identity"]
df_combined.to_csv(sys.stdout, sep="\t", header=True, index=False, columns=cols_order)
