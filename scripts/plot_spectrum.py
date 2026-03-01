import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


def get_kvmer_error_type_spectrum(df):
    # insertion, deletion, substitution
    spectrum = {"Insertion": 0, "Deletion": 0, "Substitution": 0}
    for col in df.columns:
        if ">-" in col:
            spectrum["Deletion"] += df[col].item()
        elif "->" in col:
            spectrum["Insertion"] += df[col].item()
        elif ">" in col:
            spectrum["Substitution"] += df[col].item()
    
    return spectrum


def plot_spectrum(skiver_report_file, output_file, normalize=True):
    color = 'slategray'

    # Range of the heatmap values for consistent coloring across plots
    cmap = sns.light_palette(color, as_cmap=True)

    cols = ["A", "C", "G", "T", "-"]
    bases = ["A", "C", "G", "T"]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    ## First subplot: spectrum as a matrix
    spectrum_df = pd.read_csv(skiver_report_file)
    
    # Estimate per-base error rate
    est_lambda = spectrum_df["lambda"].item()
    est_lambda_range = spectrum_df["lambda_5-95th_percentile"].item().split("~")
    est_lambda_range = [float(x) for x in est_lambda_range]
    per_base_error_rate = 1 - np.exp(-est_lambda)
    per_base_error_rate_range = [1 - np.exp(-est_lambda_range[0]), 1 - np.exp(-est_lambda_range[1])]
    print(f"Estimated per-base error rate: {per_base_error_rate:.4f} (95% CI: {per_base_error_rate_range[0]:.4f} ~ {per_base_error_rate_range[1]:.4f})")


    spectrum_matrix = np.zeros((5, 5))
    for from_base in cols:
        for to_base in cols:
            for prev_base in bases:
                for next_base in bases:
                    if from_base == to_base:
                        continue
                    key = f"{prev_base}[{from_base}>{to_base}]{next_base}"
                    spectrum_matrix[cols.index(from_base), cols.index(to_base)] += spectrum_df[key].item()
        
    # Normalize
    if normalize:
        spectrum_matrix /= spectrum_matrix.sum()
        print("The spectrum is normalized so that the sum of all entries is 1.")
    else:
        # Calculate the per-base error frequencies
        spectrum_matrix /= spectrum_matrix.sum()
        spectrum_matrix *= per_base_error_rate
    
    spectrum_matrix *= 100  # Convert to percentage for better readability


    # Plot the headmap
    ax1 = axes[0]
    mask_matrix = np.eye(5, dtype=bool)
    sns.heatmap(spectrum_matrix, annot=True, fmt=".3f", 
                xticklabels=cols, yticklabels=cols,
                cbar=True, mask=mask_matrix,
                linewidths=3, linecolor='white',
                cmap=cmap, ax=ax1)
    
    if normalize:
        ax1.set_title("Normalized error spectrum (%)")
    else:
        ax1.set_title("Error rates (%)")
    ax1.set_ylabel("Original base")
    ax1.set_xlabel("Observed base")
    plt.tight_layout()

    ## Second subplot: substitution/indel rates
    ax2 = axes[1]
    error_rate_df = pd.DataFrame(get_kvmer_error_type_spectrum(spectrum_df), index=[0])
    #print(error_rate_df)

    total_sum = error_rate_df.sum(axis=1).item()
    
    if normalize:
        error_rate_df /= total_sum
    else:
        error_rate_df /= total_sum
        error_rate_df *= per_base_error_rate
    
    error_rate_df *= 100  # Convert to percentage for better readability

    print(error_rate_df)

    #error_rate_df = error_rate_df.melt(var_name="Error Type", value_name="Frequency")
    #bars = ax.barplot(data=error_rate_df, x="Error Type", y="Frequency", ax=ax)
    x = range(len(error_rate_df.columns))
    x_labels = ["Insertion", "Deletion", "Substitution"]
    error_rates = error_rate_df.iloc[0].values
    ax2.set_xticks(x)
    ax2.set_xticklabels(x_labels)
    bars = ax2.bar(x, error_rates, color=color, width=0.6)
    ax2.bar_label(bars, fmt='%.3f', padding=3)

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # set y-axis limit to maximum error rate + 10% margin
    max_error_rate = max(error_rates)
    ax2.set_ylim(0, max_error_rate * 1.1)

    if normalize:
        ax2.set_title("Normalized error type distribution (%)")
    else:
        ax2.set_title("Error type distribution (%)")
    ax2.set_ylabel("Frequency (%)")
    ax2.set_xlabel("Error Type")

    plt.subplots_adjust(wspace=0.2, hspace=0.3)
    plt.savefig(output_file)
    #plt.show()
    plt.close()


if __name__ == "__main__":
    #skiver_report_file = "./skiver_report.csv"
    #output_file = "./output.png"
    #plot_spectrum(skiver_report_file, output_file, normalize=False)

    import argparse
    parser = argparse.ArgumentParser(description="Plot the error spectrum from Skiver report.")
    parser.add_argument("skiver_report_file", help="Path to the Skiver report CSV file.")
    parser.add_argument("output_file", help="Path to save the output plot image.")
    parser.add_argument("--normalize", action="store_true", help="Whether to normalize the spectrum (default: False). If not normalized, the spectrum will be scaled to reflect the estimated per-base error rate.")
    args = parser.parse_args()

    plot_spectrum(args.skiver_report_file, args.output_file, normalize=args.normalize)


