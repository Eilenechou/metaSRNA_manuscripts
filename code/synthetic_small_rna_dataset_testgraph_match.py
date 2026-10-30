import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
from Bio import SeqIO

def main():
    csv_file = "csv_data/eq_map_test5_tolerance2_frequency.csv"
    df = pd.read_csv(csv_file,sep="\t")
    csv_whole_file = "data/whole_paper_sequences.csv"
    df_w = pd.read_csv(csv_whole_file,sep="\t")
    #print(df['sequence'])
    merged_df = df_w.merge(df, on='original_sequence', how='left')
    merged_df.fillna(0, inplace=True)
    print(merged_df)
    #number of subplots (only 1 in this case)
    num_columns = 1
    fig, axes = plt.subplots(num_columns, 1, figsize=(10, 6 * num_columns))

    #since there's only one subplot, `axes` is not an array but a single Axes object
    if num_columns == 1:
        axes = [axes]

    #plot histogram on the first subplot
    axes[0].bar(merged_df.index, merged_df['count'].values, edgecolor='black', alpha=0.7)

    #set y-axis tick spacing to 250,000
    #axes[0].yaxis.set_major_locator(ticker.MultipleLocator(250000))
    #axes[0].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x:,.0f}'))

    #set titles and labels
    axes[0].set_title('Frequency Plot for Sequences')
    axes[0].set_xlabel('Sequence id from 1 to 1152')
    axes[0].set_ylabel('Frequency')

    #adjust layout to prevent overlap
    plt.tight_layout()

    #save the plot
    plt.savefig("csv_data/frequency_plot_test5_map_tolerance2.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
