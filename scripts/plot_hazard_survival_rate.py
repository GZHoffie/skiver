import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


def plot_hazard_survival_rate(hazard_rate_file, skiver_report_file, output_file, t_min=1, t_max=100, log_scale=False):
    color1 = 'slategray'
    color2 = 'indianred'

    # Load the Skiver report data
    report_df = pd.read_csv(skiver_report_file)
    hr_df = pd.read_csv(hazard_rate_file)

    # Extract the estimated lambda and beta
    est_lambda = report_df["lambda"].item()
    est_beta = report_df["beta"].item()

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)

    # First subplot: hazard rate curve
    plt.sca(axes[0])
    
    plt.plot(hr_df["t"], hr_df["hazard_ratio"], color=color1, linewidth=3, label="Estimated hazard rate")
    # plt the range between hr_df["5th_percentile"] and hr_df["95th_percentile"] as a shaded area
    plt.fill_between(hr_df["t"], hr_df["5th_percentile"], hr_df["95th_percentile"], color=color1, alpha=0.3, label="5%-95%  percentile")

    # Plot the estimated hazard rate curve based on the estimated lambda and beta
    t_values = hr_df["t"].values
    est_hazard_rate = 1 - np.exp(- est_lambda * (t_values ** est_beta - (t_values - 1) ** est_beta))
    plt.plot(t_values, est_hazard_rate, color=color2, linestyle='--', linewidth=3, label="Fitted hazard rate")
    plt.legend()
    if log_scale:
        plt.yscale('log') 


    plt.title("Estimated hazard rate")
    plt.xlabel("$t$")
    plt.ylabel("$h(t)$")

    # Second subplot: survival rate curve
    t_values_all = np.arange(t_min, t_max + 1)
    est_survival_rate = np.exp(- est_lambda * (t_values_all ** est_beta))
    # print the estimated survival rate as a table
    print("t,estimated_survival_rate")
    for t, s in zip(t_values_all, est_survival_rate):
        print(f"{t},{s:.6f}")

    plt.sca(axes[1])
    plt.plot(t_values_all, est_survival_rate, color=color2, linestyle='--', linewidth=3, label="Estimated survival rate \n from fitted hazard rate")
    plt.title("Estimated survival rate")
    plt.xlabel("$t$")
    plt.ylabel("$S(t)$")
    if log_scale:
        plt.yscale('log') 
    plt.legend()

 
    plt.tight_layout()
    plt.savefig(output_file)

    #plt.show()

if __name__ == "__main__":
    #hazard_rate_file = "./hazard_rate.csv"
    #skiver_report_file = "./skiver_report.csv"
    #output_file = "./coverage_histogram.png"
    #plot_hazard_survival_rate(hazard_rate_file, skiver_report_file, output_file)

    import argparse
    parser = argparse.ArgumentParser(description="Estimate the true coverage histogram from Skiver report.")
    parser.add_argument("hazard_rate_file", help="Path to the hazard rate CSV file.")
    parser.add_argument("skiver_report_file", help="Path to the Skiver report CSV file.")
    parser.add_argument("output_file", help="Path to save the output plot image.")
    parser.add_argument("--log_scale", action="store_true", help="Use logarithmic scale for y-axis.")
    parser.add_argument("-t", type=int, default=1, help="Minimum t value for survival rate curve (default: 1).")
    parser.add_argument("-T", type=int, default=100, help="Maximum t value for survival rate curve (default: 100).")
    args = parser.parse_args()
    plot_hazard_survival_rate(args.hazard_rate_file, args.skiver_report_file, args.output_file, t_min=args.t, t_max=args.T, log_scale=args.log_scale)