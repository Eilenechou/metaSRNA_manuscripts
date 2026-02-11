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
    fasta_files_names =["data/final_seq_12_SRR5234257_4N_A.fasta","data/final_seq_12_SRR5234277_4N_NEXTflex.fasta",
    "data/final_seq_12_SRR5234287_CleanTag.fasta","data/final_seq_12_SRR5234307_TrueSeq.fasta","data/final_seq_12_SRR6380716_NEBNext.fasta"]  
    for name in fasta_files_names:
        seq_dict = SeqIO.to_dict(SeqIO.parse(name, "fasta"))
        sequence_list = [str(record.seq) for record in seq_dict.values()]
        df = pd.DataFrame(sequence_list, columns=["sequence"])
        print("sample name", name)
        print("sequence length",len(df))
        # dfb = df["sequence"]
        dfa = df["sequence"].value_counts().reset_index()
        dfa.columns = ["sequence", "count"]
        print("remove duplicates length",len(dfa))
      


    

if __name__ == "__main__":
    main()  


