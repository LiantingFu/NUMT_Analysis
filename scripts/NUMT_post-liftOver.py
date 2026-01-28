#!/usr/bin/env python

import sys
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="post-liftover")
parser.add_argument("--input_records", "-r", help="original NUMT result before liftover")
parser.add_argument("--input_bed", "-b", help="liftOvered NUMT bed")

args = parser.parse_args()

dtypes = {"chrom_adjust": "str", "mt_start_adjust": "int32", "mt_end_adjust": "int32", "numt_id": "str", "score_adjust": "str", "strand_adjust": "str"}
df_result = pd.read_csv(args.input_records,sep="\t",header=0)
df_bed = pd.read_csv(args.input_bed,sep="\t",header=None,names=["chrom_adjust", "mt_start_adjust", "mt_end_adjust", "numt_id", "score_adjust", "strand_adjust"],dtype=dtypes)

# convert human chrM coordinate (repeats3); keep 0-based coordinate
index1 = (16569 <= df_bed["mt_start_adjust"]) & (df_bed["mt_start_adjust"] < 33138)
index2 = df_bed["mt_start_adjust"] >= 33138
df_bed.loc[index1, ["mt_start_adjust", "mt_end_adjust"]] -= 16569
df_bed.loc[index2, ["mt_start_adjust", "mt_end_adjust"]] -= 33138

df_combined = pd.merge(df_result,df_bed,on="numt_id",how="left")

df_combined.to_csv(sys.stdout, sep="\t", index=False, na_rep="NaN")
