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

def map_to_seq_b(seq_a,df_csv):
    for seq_b in df_csv['sequence']:
        if seq_a in seq_b:
            return seq_b  # Return the matching seq_b
    return ''  # Return an empty string if no match is found

def main():
    record_dict = SeqIO.to_dict(SeqIO.parse("data/final_seq_12_SRR6380716_NEBNext.fasta", "fasta"))
    sequence_list = [str(record.seq) for record in record_dict.values() if 16 <= len(record.seq) <= 25]
    df = pd.DataFrame(sequence_list, columns=["sequence"])
    dfa = df["sequence"]
    # dfa = dfb.drop_duplicates().reset_index(drop=True)
    num_lines = len(dfa)
    print(f"The number of sequences is: {num_lines}")
    print("original mirdeep2 sequences:\t", dfa)
    num_lines = len(dfa)
    print(f"The number of sequences is: {num_lines}")

    df_csv = pd.read_csv("data/1152_paper_sequences.csv", sep="\t")
    df_csv.columns =["sequence","unkown"]
    # print(df_csv)
    # df_subset_csv = df_csv['sequence']
    # print(df_subset_csv)

    # df.loc[:, 'is_substring'] = dfa.apply(
    #     lambda seq_a: any(seq_a in seq_b for seq_b in df_subset_csv)
    # )
   


    df['mapped_seq_b'] = dfa.apply(map_to_seq_b, args=(df_csv,))

    mapped_seq_b_counts = df['mapped_seq_b'].value_counts()

    mapped_seq_b_counts_df = mapped_seq_b_counts.reset_index()
    mapped_seq_b_counts_df.columns = ['original_sequence', 'count']

    # df['is_substring'] = df['mapped_seq_b'] != ''

    # result = df[df['is_substring']]
    
    
    mapped_seq_b_counts_df.loc[mapped_seq_b_counts_df['original_sequence'] == '', 'count'] = 0
    real_count = mapped_seq_b_counts_df['count'].sum()
    print("real count",real_count)
   
    count = len(dfa)
    print("df_csv length", count)
    print("mapping percentage",round((real_count /count) * 100,2))
    mapped_seq_b_counts_df.to_csv('data/final_seq_12_SRR6380716_NEBNext_match.csv', sep='\t', index=False)
    # Display the result
    # print(result)

    # not_substring_result = df[~df['is_substring']]
    # print(not_substring_result)
    # # not_substring_result = not_substring_result.copy()
    # # not_substring_result.drop_duplicates(subset=["sequence"], inplace=True)
    # # print(not_substring_result)
    # not_substring_result.to_csv('data/final_seq_12_SRR6380716_NEBNext_unmatch.csv', sep='\t', index=False)

    # csv_file = "eq_not_map_test5_frequency.csv"  # Replace with your CSV file path
    # df = pd.read_csv(csv_file,sep="\t")
    # print(df['sequence'])

if __name__ == "__main__":
    main()


