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

    key_coverages = np.array(coverage_df["total_count"].values) / est_S_k

    # Exclude the top 0.01% of coverages for better visualization
    coverage_threshold = np.percentile(key_coverages, 99.99)
    key_coverages = key_coverages[key_coverages <= coverage_threshold]

    # Print the estimated true coverage (median, and 5-95th percentile)
    median_coverage = np.median(key_coverages)
    coverage_5th_percentile = np.percentile(key_coverages, 5)
    coverage_95th_percentile = np.percentile(key_coverages, 95)
    print(f"Estimated true coverage (median): {median_coverage:.2f}")
    print(f"Estimated true coverage (5-95th percentile): {coverage_5th_percentile:.2f} ~ {coverage_95th_percentile:.2f}")




    # Plot the histogram of key coverages
    plt.figure(figsize=(10, 4))
    plt.rc('axes.spines', **{'bottom':True, 'left':True, 'right':False, 'top':False})
    sns.histplot(data=key_coverages, bins=100, color='slategray', edgecolor='white', stat='probability', log_scale=(False, True))

    plt.xlabel("Coverage")
    plt.ylabel("Probability")
    plt.title("Estimated true coverage")
    plt.yscale('log')  # Use logarithmic scale for better visibility of low-frequency bins
    plt.tight_layout()
    plt.savefig(output_file)

    #plt.show()

if __name__ == "__main__":
    #verbose_output_file = "./verbose_output.csv"
    #skiver_report_file = "./skiver_report.csv"
    #output_file = "./coverage_histogram.png"
    #plot_coverage_histogram(verbose_output_file, skiver_report_file, output_file)

    import argparse
    parser = argparse.ArgumentParser(description="Estimate the true coverage histogram from Skiver report.")
    parser.add_argument("verbose_output_file", help="Path to the verbose output CSV file from Skiver.")
    parser.add_argument("skiver_report_file", help="Path to the Skiver report CSV file.")
    parser.add_argument("output_file", help="Path to save the output plot image.")
    args = parser.parse_args()
    plot_coverage_histogram(args.verbose_output_file, args.skiver_report_file, args.output_file)