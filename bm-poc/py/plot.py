import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('miacl-eval.csv', comment='#')

df.columns = ['Signers', 'Attrs', 'Issuance', 'Token Verification', 'Credential Presentation', 'Credential Verification']

def plot_metric(df, metric_name):
    plt.figure(figsize=(8, 6))
    
    # For each unique number of attributes
    for attr in sorted(df['Attrs'].unique()):
        subset = df[df['Attrs'] == attr]
        subset = subset.sort_values(by='Signers')  # Ensure x-axis order
        plt.plot(subset['Signers'], subset[metric_name], marker='o',markersize=2, linewidth=1, label=f'{attr} Attrs', linestyle="--")

    plt.title(metric_name)
    plt.xlabel('Number of Signers')
    plt.ylabel('Time (ms)')
    
    plt.legend(
        title='Attributes'
    )
    
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{metric_name.replace(" ", "_").lower()}.png')
    plt.show()

# Generate the four plots
for metric in ['Issuance', 'Token Verification', 'Credential Presentation', 'Credential Verification']:
    plot_metric(df, metric)