from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
import os
import glob
import pandas as pd
def filter_fasta_by_txt(txt_file, fa_file):
    """
    1. Read a comma-delimited txt file, take unique values from the first column (numbers).
    2. Read a fasta file, find sequences whose header number > 1 among those unique numbers.
    3. Plot bar chart and line plot of sequence counts.
    Returns filtered sequences as a dict {name: sequence}.
    """
    #step 1: Read txt, get unique first-column values
    df = pd.read_csv(txt_file, header=None, sep=",")
    filtered_nums = set(df.iloc[:, 0].dropna().astype(int).unique())
    print("filtered_nums", filtered_nums)

    sequences = {}
    for record in SeqIO.parse(fa_file, "fasta"):
        r_id = int(record.id)
        if r_id in filtered_nums:
            sequences[r_id] = str(record.seq)

    print(f"Filtered FASTA sequences: {len(sequences)}")

    #Step 3: Count how many sequences have each length
    
    length_counts = Counter(len(seq) for seq in sequences.values())
    print(length_counts)
    max_len = max(length_counts.keys())
    x = list(range(1, max_len + 1))
    y = [length_counts.get(l, 0) for l in x]
    print(y)
    # Save x, y to csv
    df_out = pd.DataFrame({"length": x, "count": y})
    df_out.to_csv("length_distribution.csv", index=False)

filter_fasta_by_txt("16th_new_paper_29211829_results/original_100_100_mapping_results/map_genome/blast_score_filter.txt", "16th_new_paper_29211829_results/original_100_100_mapping_results/final_seq_12.fasta")