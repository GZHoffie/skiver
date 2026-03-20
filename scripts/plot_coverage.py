import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


def plot_coverage_histogram(verbose_output_file, skiver_report_file, output_file):
    # Load the Skiver report data
    report_df = pd.read_csv(skiver_report_file)
    coverage_df = pd.read_csv(verbose_output_file)

    # Extract the estimated lambda and beta
    est_lambda = report_df["lambda"].item()
    est_beta = report_df["beta"].item()

    # Collect k
    k = len(coverage_df.iloc[0]["key"])

    # Estimate S(k)
    est_S_k = np.exp(-est_lambda * (k ** est_beta))

    print("k =", k)
    print(f"Estimated S(k) = {est_S_k:.4f}")

    all_coverages = np.array(coverage_df["total_count"].values) / est_S_k
    passing_mask = coverage_df["passes_filter"].values.astype(bool)
    passing_coverages = all_coverages[passing_mask]

    # Exclude the top 0.01% of all coverages for better visualization
    coverage_threshold = np.percentile(all_coverages, 99.99)
    all_coverages = all_coverages[all_coverages <= coverage_threshold]
    passing_coverages = passing_coverages[passing_coverages <= coverage_threshold]

    # Print the estimated true coverage (median, and 5-95th percentile)
    median_coverage = np.median(passing_coverages)
    coverage_5th_percentile = np.percentile(passing_coverages, 5)
    coverage_95th_percentile = np.percentile(passing_coverages, 95)
    print(f"Estimated true coverage (median): {median_coverage:.2f}")
    print(f"Estimated true coverage (5-95th percentile): {coverage_5th_percentile:.2f} ~ {coverage_95th_percentile:.2f}")

    # Shared bins so both histograms align on the same ticks
    bin_min = all_coverages.min()
    bin_max = all_coverages.max()
    bins = np.linspace(bin_min, bin_max, 101)  # 100 bins of equal width

    # Plot the histogram of key coverages
    plt.figure(figsize=(10, 4))
    plt.rc('axes.spines', **{'bottom':True, 'left':True, 'right':False, 'top':False})
    sns.histplot(data=all_coverages, bins=bins, color='slategray', edgecolor='none', stat='count', alpha=0.5, label='All keys')
    sns.histplot(data=passing_coverages, bins=bins, color='steelblue', edgecolor='none', stat='count', alpha=0.7, label='Passing filter')

    plt.xlabel("Coverage")
    plt.ylabel("Count")
    plt.title("Estimated true coverage")
    plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file)

    #plt.show()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Estimate the true coverage histogram from Skiver report.")
    parser.add_argument("kvmer_csv", help="Path to the verbose kvmer summary CSV file from Skiver.")
    parser.add_argument("summary_error_rate_csv", help="Path to the Skiver error rate report CSV file.")
    parser.add_argument("output_file", help="Path to save the output plot image.")
    args = parser.parse_args()
    plot_coverage_histogram(args.kvmer_csv, args.summary_error_rate_csv, args.output_file)