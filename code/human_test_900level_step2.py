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

#fastq file sequence should only substring of 1152 sequences
def find_matching_seqs(seq_1,df_csv):
    for seq_2 in df_csv['sequence']:
        if len(seq_2) >= len(seq_1):
            if seq_1 in seq_2:
                return True

        
    return False  
# # Function to check if seq_1 is a near-substring (tolerating one error) of any seq_2
def find_matching_seqs_2(seq_1,dfa):
    for seq_2 in dfa['sequence']:
        if len(seq_2) <= len(seq_1):
            if seq_2 in seq_1:
                return True

        
    return False  # No match found


def main():
    record_dict = SeqIO.to_dict(SeqIO.parse("data/final_seq_12_SRR5234257_4N_A.fasta", "fasta"))
    sequence_list = [str(record.seq) for record in record_dict.values() if 16 <= len(record.seq) <= 25]
    df = pd.DataFrame(sequence_list, columns=["sequence"])
    # dfb = df["sequence"]
    dfa = df["sequence"].value_counts().reset_index()
    dfa.columns = ["sequence", "count"]
    print("dfa",dfa.head(5))
    # num_lines = len(dfa)
    # print(f"The number of sequences is: {num_lines}")
    # print("original mirdeep2 sequences:\t", dfa)
    # num_lines = len(dfa)
    # print(f"The number of sequences is: {num_lines}")

    df_csv = pd.read_csv("data/shorter_paper_sequences.csv", sep="\t")
    df_csv.columns =["sequence"]
   

    dfa['matched_seqs'] = dfa['sequence'].apply(lambda seq: find_matching_seqs(seq, df_csv))
    df_match_count = dfa['matched_seqs'].sum()
    dfa['weighted_match'] = dfa['matched_seqs'].astype(int) * dfa['count']
    df_match_count_all = dfa['weighted_match'].sum()
    
    # df_matched = dfa[dfa['matched_seqs'].apply(len) > 0].copy()
    # df_unmatched = dfa[dfa['matched_seqs'].apply(len) == 0].copy()
    
    count = len(dfa)
    total_count = dfa['count'].sum()
    print("total sample count",total_count)
    print("mapping count",df_match_count_all)
    print("mapping percentage",round((df_match_count_all /total_count) * 100,2))
    print("total sample count remove duplicates", count)
    print("mapping count remove duplicates",df_match_count)
    print("not mapped count",count-df_match_count)
    print("mapped percentage",round(100-(count-df_match_count)/count * 100,2))
    print("not mapped percentage",round((count-df_match_count)/count * 100,2))
    df_csv['matched_seqs'] = df_csv['sequence'].apply(lambda seq: find_matching_seqs_2(seq, dfa))
    print("cover paper sequence count",df_csv['matched_seqs'].sum())
    print("not cover paper sequence count",977 - df_csv['matched_seqs'].sum())
    print("mapping percentage",round((df_csv['matched_seqs'].sum() /977) * 100,2))
    print("not mapping percentage",round(100 - (df_csv['matched_seqs'].sum() /977) * 100,2))
    # print("mapping percentage",round((len(df_matched) /1152) * 100,2))
    # df_matched.to_csv('data/final_seq_12_SRR6380716_NEBNext_match.csv', sep='\t', index=False)
    # Display the result
    # print(result)

    # not_substring_result = df[~df['is_substring']]
    # print(not_substring_result)
    # # not_substring_result = not_substring_result.copy()
    # # not_substring_result.drop_duplicates(subset=["sequence"], inplace=True)
    # # print(not_substring_result)
    # df_unmatched.to_csv('data/final_seq_12_SRR6380716_NEBNext_unmatch.csv', sep='\t', index=False)
    # # Apply the function to df_1['sequence']
    # df_csv['total_match'] = df_csv['sequence'].apply(lambda seq: is_near_substring(seq, df_unmatched))
    # print(df_csv)
    # total_match_count = df_csv['total_match'].sum()
    # print("second mapping count",total_match_count)
    # print("second mapping percentage",round(((total_match_count) /1152) * 100,2))
    # Print the result
    # print(df_csv)
    # csv_file = "eq_not_map_test5_frequency.csv"  # Replace with your CSV file path
    # df = pd.read_csv(csv_file,sep="\t")
    # print(df['sequence'])

if __name__ == "__main__":
    main()


