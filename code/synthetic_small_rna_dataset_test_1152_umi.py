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


def find_matching_seqs(seq_1, df_csv):
    for seq_id, row in df_csv.iterrows():
        seq_2 = row['sequence']
        # seq_id = row['id']  # assuming df_csv has a column named 'id'
        
        if len(seq_2) >= len(seq_1):
            if seq_1 in seq_2:
                return seq_id
            for i in range(len(seq_2) - len(seq_1) + 1):
                sub_seq = seq_2[i:i+len(seq_1)]
                if Levenshtein.distance(seq_1, sub_seq) <= 1:
                    return seq_id
        # Optional: handle len(seq_2) < len(seq_1) here if needed
    return None


def main():
    seq_dict = SeqIO.to_dict(SeqIO.parse("data/final_seq_12_SRR5234257_4N_A.fasta", "fasta"))
    umi_dict = SeqIO.to_dict(SeqIO.parse("data/final_umi_SRR5234257_4N_A.fasta", "fasta"))

    # Create list of rows: (seq_id, main_sequence, umi_sequence)
    data = []
    for seq_id in seq_dict:
        if seq_id in umi_dict:
            seq = str(seq_dict[seq_id].seq)
            umi = str(umi_dict[seq_id].seq)
            data.append((seq, umi))

    # Convert to DataFrame
    df = pd.DataFrame(data, columns=["sequence", "umi"])
    df_dedup = df.groupby(["sequence", "umi"]).size().reset_index(name="count")
    df.to_csv("seq_with_umi_4N_A.csv", index=False)

    # sequence_list = [str(record.seq) for record in record_dict.values() if 12 <= len(record.seq) <= 90]
    # df = pd.DataFrame(sequence_list, columns=["sequence"])
    # # dfb = df["sequence"]
    # dfa = df.drop_duplicates().reset_index(drop=True)
    # print("dfa",dfa.head(5))
    # num_lines = len(dfa)
    # print(f"The number of sequences is: {num_lines}")
    # print("original mirdeep2 sequences:\t", dfa)
    # num_lines = len(dfa)
    # print(f"The number of sequences is: {num_lines}")

    df_csv = pd.read_csv("data/1152_paper_sequences.csv", sep="\t")
    df_csv.columns =["sequence","unkown"]
   

    df_dedup['seq_id'] = df_dedup['sequence'].apply(lambda seq: find_matching_seqs(seq, df_csv))
    df_grouped = df_dedup.groupby('seq_id')['umi'].apply(lambda x: set(x)).reset_index()
    # Step 2: count how many unique UMIs per set
    df_grouped['umi_count'] = df_grouped['umi'].apply(len)

    # Step 3: sort by seq_id
    df_grouped = df_grouped.sort_values(by='seq_id').reset_index(drop=True)
    df_grouped.to_csv("final_umi_count_4N_A.csv", index=False)
    plt.figure(figsize=(14, 5))  # wide to fit many seq_ids
    df_grouped['x_axis'] = df_grouped['seq_id'] + 1
    plt.bar(df_grouped['x_axis'], df_grouped['umi_count'], color='green', edgecolor='green')
    plt.xlabel('seq_id')
    plt.ylabel('UMI count')
    plt.title('UMI count per seq_id among 1152 sequences')
    # plt.xticks(rotation=90, fontsize=6)  # rotate labels if too dense
    max_id = int(df_grouped['x_axis'].max())
    max_y = int(df_grouped['umi_count'].max())
    plt.xticks(rotation=0, ha='right')  # horizontal ticks
    plt.xlim(left=0)
    plt.xticks(ticks=range(0, max_id+100, 100))
    plt.yticks(ticks=range(0, max_y+1, 1))
    plt.tight_layout()
    plt.savefig("umi_count_histogram_4N_A.png", dpi=300)  # Save as PNG
  
if __name__ == "__main__":
    main()


