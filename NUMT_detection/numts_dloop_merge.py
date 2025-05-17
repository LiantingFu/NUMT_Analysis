#!/usr/bin/env python

import sys
import pandas as pd
import portion as P
import argparse


parser = argparse.ArgumentParser(description="merge numts divided in D-loop region")
parser.add_argument('input_result', type=argparse.FileType('r'), default=sys.stdin, help="Blastn result, file or standard input")
parser.add_argument('mtDNA_len', type=int, help="Length of mtDNA sequence")
parser.add_argument('species', help="species of mtDNA")
args = parser.parse_args()

input_result = args.input_result
mtDNA_len = args.mtDNA_len
species = args.species


def adjust_strand(row):
    if row['nu_start'] > row['nu_end']:
        row['nu_start'], row['nu_end'] = row['nu_end'], row['nu_start']
        row['strand'] = '-'
    else:
        row['strand'] = '+'
    return row


# record1: spanning; record2: front
def merge_dloop(record1, record2):
    if record1['nu_chrom'] != record2['nu_chrom']:
        return False
    nu_region1 = P.closed(record1['nu_start'], record1['nu_end'])
    nu_region2 = P.closed(record2['nu_start'], record2['nu_end'])
    nu_intersect_region = nu_region1 & nu_region2
    if nu_intersect_region.empty:
        flag = False
        return flag

    nu_intersect_rate = (nu_intersect_region.upper - nu_intersect_region.lower) / (nu_region2.upper - nu_region2.lower)

    mt_region1 = P.closed(1, record1['mt_end'] - mtDNA_len)
    mt_region2 = P.closed(record2['mt_start'], record2['mt_end'])
    mt_intersect_region = mt_region1 & mt_region2
    if mt_intersect_region.empty:
        flag = False
        return flag

    mt_intersect_rate = (mt_intersect_region.upper - mt_intersect_region.lower) / (mt_region2.upper - mt_region2.lower)

    if nu_intersect_rate >= 0.95 and mt_intersect_rate >= 0.95 and record1['strand'] == record2['strand']:
        flag = True
    else:
        flag = False
    return flag


dtypes = {
    'mt_chrom': 'str', 'nu_chrom': 'str', 'identity': 'float32',
    'length': 'int32', 'mismatch': 'int32', 'gapopen': 'int32',
    'mt_start': 'int32', 'mt_end': 'int32', 'nu_start': 'int32',
    'nu_end': 'int32', 'evalue': 'float32', 'bitscore': 'str',
    'qcovs': 'str', 'mt_seq': 'str', 'nu_seq': 'str'
}

df_raw = pd.read_csv(input_result, sep='\t', names=['mt_chrom', 'nu_chrom', 'identity', 'length', 'mismatch', 'gapopen', 'mt_start', 'mt_end', 'nu_start', 'nu_end', 'evalue', 'bitscore', 'qcovs', 'mt_seq', 'nu_seq'], dtype=dtypes, comment='#')

df_input = df_raw[(df_raw['identity'] >= 60) & (df_raw['evalue'] <= 1e-3) & (~df_raw['nu_chrom'].str.contains('chrM'))][['nu_chrom', 'nu_start', 'nu_end', 'mt_chrom', 'mt_start', 'mt_end', 'identity']].copy()

df_input = df_input.apply(adjust_strand, axis=1)
df_input.drop_duplicates(inplace=True)
df_input.sort_values(by=['nu_chrom', 'nu_start'], inplace=True)

df_front = df_input[df_input['mt_end'] < mtDNA_len].copy()

df_spanning = df_input[(df_input['mt_start'] <= mtDNA_len) & (df_input['mt_end'] >= mtDNA_len)].copy()

df_front['discard'] = False

for index2, row2 in df_front.iterrows():
    for index, row in df_spanning.iterrows():
        df_front.at[index2, 'discard'] = merge_dloop(row, row2)
        if df_front.at[index2, 'discard']:
            break

df_front_not_discarded = df_front[~df_front['discard']].drop(columns=['discard'])
df_combined = pd.concat([df_spanning, df_front_not_discarded], ignore_index=True)
df_combined.sort_values(by=["nu_chrom", "nu_start"], inplace=True)
df_combined.reset_index(drop=True, inplace=True)
df_combined['id'] = 'numt_' + (df_combined.index + 1).astype(str) + '_' + species

output_result = df_combined.to_csv(sep='\t', header=True, index=False, path_or_buf=None)
sys.stdout.write(output_result)
