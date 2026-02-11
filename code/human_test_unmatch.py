import pandas as pd
import Levenshtein

def is_near_substring(seq_1):
    for seq_2 in df_2['original_sequence']:
        for i in range(len(seq_2) - len(seq_1) + 1):  # Slide seq_1 over seq_2
            sub_seq = seq_2[i:i+len(seq_1)]  # Extract substring of same length
            if Levenshtein.distance(seq_1, sub_seq) <= 2:  # Allow one error
                return True  # Match found
    return False  # No match found

def main():
    df_1 = pd.read_csv("data/whole_paper_sequences.csv", sep="\t")
    df_2 = pd.read_csv("csv_data/eq_not_map_test1_frequency.csv", sep="\t")
 
    df_1['original_sequence'] = df_1['sequence'].apply(is_near_substring)
    print(df_1)

   
    print(df_1)

if __name__ == "__main__":
    main()