import pandas as pd
import Levenshtein
# Function to check if seq_1 is a near-substring (tolerating one error) of any seq_2
def is_near_substring(seq_1,df_2):
    for seq_2 in df_2['original_sequence']:
        for i in range(len(seq_2) - len(seq_1) + 1):  # Slide seq_1 over seq_2
            sub_seq = seq_2[i:i+len(seq_1)]  # Extract substring of same length
            if Levenshtein.distance(seq_1, sub_seq) <= 2:  # Allow one error
                return seq_2  # Match found
    return ''  # No match found

def main():
    df_1 = pd.read_csv("csv_data/eq_not_map_test5_frequency.csv", sep="\t")
    df_2 = pd.read_csv("data/whole_paper_sequences.csv", sep="\t")
    
    df_1['original_sequence'] = df_1['sequence'].apply(is_near_substring, args=(df_2,))
    real_count= df_1.loc[df_1['original_sequence'] != '', 'count'].sum()
    count = df_1['count'].sum()

    # count = len(df_1)
    print("mapping percentage",round((real_count /count) * 100,2))
    print(df_1)
    df_1_not = df_1[df_1['original_sequence'] == '']
    num_rows_not_count_1 = count - (df_1_not['count'] == 1).sum()
    print("final unmapping percentage",round((df_1_not['count'].sum() /count) * 100,2))
    print("final unmapping percentage not count 1",round((num_rows_not_count_1 /count) * 100,2))
    # df_1_not = df_1[['sequence','count']]
    df_1_not.to_csv('eq_notmap_test5_tolerance2_frequency.csv', sep='\t', index=False)
    df_1=df_1[df_1['original_sequence'] != '']
    print(df_1)
    df_1 = df_1[['original_sequence','count']]
    df_1new = df_1.groupby('original_sequence', as_index=False)['count'].sum()
    print(df_1new)
    df_3 = pd.read_csv("csv_data/eq_map_test5_new.csv", sep="\t")
    df_3 = df_3[df_3['original_sequence'].notna()]
    # print(df_3)
    # df_1new = df_1[['original_sequence','count']]
    #print(df_1new)
    
    merged_df = df_1new.merge(df_3, on='original_sequence', how='outer', suffixes=('_df1', '_df2'))
    print(merged_df)
   
    merged_df['count'] = merged_df['count_df1'].fillna(0) + merged_df['count_df2'].fillna(0)
    print(merged_df)

    merged_df = merged_df[['original_sequence', 'count']]
    print(merged_df)
    final_df = df_2.merge(merged_df, on='original_sequence', how='left')
    print(final_df)
    final_df.to_csv('eq_map_test5_tolerance2_frequency.csv', sep='\t', index=False)
if __name__ == "__main__":
    main()