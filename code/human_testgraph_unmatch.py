import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
from Bio import SeqIO

def main():
    csv_file = "csv_data/eq_notmap_test5_tolerance2_frequency.csv"  
    df = pd.read_csv(csv_file,sep="\t")
    # print(df['sequence'])
    # Number of subplots (only 1 in this case)
    num_columns = 1
    fig, axes = plt.subplots(num_columns, 1, figsize=(10, 6 * num_columns))  

    
    if num_columns == 1:
        axes = [axes] 

   
    axes[0].bar(df.index, df['count'].values, edgecolor='black', alpha=0.7)

    # Set y-axis tick spacing to 250,000
    # axes[0].yaxis.set_major_locator(ticker.MultipleLocator(250000))
    # axes[0].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x:,.0f}'))

    
    axes[0].set_title('Frequency Plot for Sequences')
    axes[0].set_xlabel('Sequence id ordered by frequency')
    axes[0].set_ylabel('Frequency')

    
    plt.tight_layout()

    
    plt.savefig("csv_data/frequency_plot_test5_notmap_torlerance2.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    main()
