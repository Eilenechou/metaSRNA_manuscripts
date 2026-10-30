import pandas as pd
import Levenshtein
#function to check if seq_1 is a near-substring (tolerating one error) of any seq_2
def is_near_substring(seq_1):
    for seq_2 in df_2['original_sequence']:
        for i in range(len(seq_2) - len(seq_1) + 1):  #slide seq_1 over seq_2
            sub_seq = seq_2[i:i+len(seq_1)]  #extract substring of same length
            if Levenshtein.distance(seq_1, sub_seq) <= 2:  #allow one error
                return True  #match found
    return False  #no match found

def main():
    df_1 = pd.read_csv("data/whole_paper_sequences.csv", sep="\t")
    df_2 = pd.read_csv("csv_data/eq_not_map_test1_frequency.csv", sep="\t")
    #apply the function to df_1['sequence']
    df_1['original_sequence'] = df_1['sequence'].apply(is_near_substring)
    print(df_1)

    #print the result
    print(df_1)

if __name__ == "__main__":
    main()