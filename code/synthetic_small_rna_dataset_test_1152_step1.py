from io import StringIO
import glob
import os
from os import listdir, path, makedirs
import re
import pandas as pd
import subprocess
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import Levenshtein
import matplotlib.pyplot as plt


def main():
    seq_dict = SeqIO.to_dict(SeqIO.parse("data/final_seq_12_SRR5234277_4N_NEXTflex.fasta", "fasta"))
    sequence_list = [str(record.seq) for record in seq_dict.values() if 12 <= len(record.seq) <= 90]
    df = pd.DataFrame(sequence_list, columns=["sequence"])
    #dfb = df["sequence"]
    dfa = df["sequence"].value_counts().reset_index()
    dfa.columns = ["sequence", "count"]
    print("dfa",dfa.head(5))

    df_csv = pd.read_csv("data/1152_paper_sequences.csv", sep="\t")
    df_csv.columns =["sequence","unkown"]
    df_csv = df_csv.drop(columns=["unkown"])
    #Step 1: Check presence of each sequence from df1 in df2
    dfa['found_in_df2'] = dfa['sequence'].isin(df_csv['sequence'])
    missing_seqs_fasta = dfa.loc[~dfa['found_in_df2'], 'sequence']
    print("Missing sequences:")
    print(len(missing_seqs_fasta))
    match_seqs_fasta = dfa.loc[dfa['found_in_df2'], 'sequence']
    print("match sequences:")
    print(len(match_seqs_fasta))
    print("total sequences:")
    print(len(dfa))

    #Step 1: Check presence of each sequence from df1 in df2
    df_csv['found_in_df2'] = df_csv['sequence'].isin(dfa['sequence'])
    missing_seqs = df_csv.loc[~df_csv['found_in_df2'], 'sequence']
    print("Missing sequences:")
    print(len(missing_seqs))
    match_seqs = df_csv.loc[df_csv['found_in_df2'], 'sequence']
    print("match sequences:")
    print(len(match_seqs))
    print("total sequences:")
    print(len(df_csv))

    umi_dict = SeqIO.to_dict(SeqIO.parse("data/final_umi_SRR5234277_4N_Next.fasta", "fasta"))

    #create list of rows: (seq_id, main_sequence, umi_sequence)
    data = []
    for seq_id in seq_dict:
        if seq_id in umi_dict:
            seq = str(seq_dict[seq_id].seq)
            umi = str(umi_dict[seq_id].seq)
            data.append((seq, umi))

    #Convert to DataFrame
    df = pd.DataFrame(data, columns=["sequence", "umi"])
    #df = df[df['sequence'].str.len().between(16, 25)]
    df_dedup = df.groupby(["sequence", "umi"]).size().reset_index(name="count")
    merged = pd.merge(df_csv, df_dedup, how='left', on='sequence')
    count_df = merged.groupby('sequence')['umi'].count().reset_index(name='umi_count')
    count_df_sorted = count_df.sort_values(by='umi_count')
    #print(count_df_sorted)
    count_df_sorted.to_csv("seq_and_umi_count_4N_Nextflex_1152.csv", index=False)

    #plot bar chart
    plt.figure(figsize=(10, 5))
    x = range(1, len(count_df_sorted) + 1)
    y = count_df_sorted["umi_count"]
    plt.figure(figsize=(10, 5))
    #plt.bar(x, y, color='green', edgecolor='green', label='Bar')
    plt.plot(x, y, marker='o', linestyle='None', color='red', markersize=1)
    plt.xlabel('Sequence Index')
    plt.ylabel('UMI Count(sorted)')
    plt.title('Sorted UMI Counts with Curve among 1152 sequences')
    plt.grid(True)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xticks(ticks=range(0, len(count_df_sorted)+100, 100))

    plt.tight_layout()
    plt.savefig("umi_count_histogram_4N_Nextflex_1152.png", dpi=300)  #Save as PNG


    

if __name__ == "__main__":
    main()  


