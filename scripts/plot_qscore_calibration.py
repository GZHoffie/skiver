import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def plot_qscore_calibration(calibration_file, output_file, log_scale=False, min_coverage=100):
    color_empirical  = 'slategray'
    color_theoretical = 'indianred'

    df = pd.read_csv(calibration_file)

    # Drop rows with zero observations
    df = df[(df["num_correct"] + df["num_error"]) > 0].copy()
    df["num_total"] = df["num_correct"] + df["num_error"]

    # Filter by minimum coverage
    df = df[df["num_total"] >= min_coverage]

    q_vals = df["qscore"].values
    emp_rate = df["error_rate"].values
    #ci_lower = df["5th_percentile"].values
    #ci_upper = df["95th_percentile"].values
    counts   = df["num_total"].values

    # Theoretical Phred error rate: P(error) = 10^(-Q/10)
    q_theory = np.linspace(1, max(q_vals), 300)
    theory_rate = 10 ** (-q_theory / 10.0)

    # --- layout: tall main panel + shorter histogram panel ---
    fig, (ax_main, ax_hist) = plt.subplots(
        2, 1,
        figsize=(8, 7),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=False,
    )

    # ── Main panel ────────────────────────────────────────────
    # Theoretical line
    ax_main.plot(q_theory, theory_rate, color="black",
                 linestyle='--', linewidth=3, label="Theoretical ($10^{-Q/10}$)", zorder=3)

    # Confidence band
    #ax_main.fill_between(q_vals, ci_lower, ci_upper,
    #                     color=color_empirical, alpha=0.25, label="5%–95% CI")

    # Empirical line
    ax_main.plot(q_vals, emp_rate, color=color_empirical,
                 linewidth=3, marker='o', 
                 label="Empirical error rate", zorder=4)

    if log_scale:
        ax_main.set_yscale('log')
    #ax_main.set_xlim(q_min, q_max)
    if log_scale:
        ax_main.set_ylabel("Error rate (log scale)")
    else:
        ax_main.set_ylabel("Error rate")
    ax_main.set_title("Quality-score calibration")
    ax_main.legend(frameon=False)
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)

    # ── Histogram panel ───────────────────────────────────────
    ax_hist.bar(q_vals, counts, width=0.8, color=color_empirical, alpha=0.7)
    #ax_hist.set_xlim(q_min, q_max)
    ax_hist.set_xlabel("Phred quality score ($Q$)")
    ax_hist.set_ylabel("# bases used for estimation")
    ax_hist.spines['top'].set_visible(False)
    ax_hist.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Plot empirical vs. theoretical error rate by Phred quality score."
    )
    parser.add_argument("summary_phred_csv",
                        help="Path to the summary Phred CSV file.")
    parser.add_argument("output_file", help="Path to save the output plot image.")
    parser.add_argument("--log", action="store_true",
                        help="Use logarithmic scale for the y-axis (default: False).")
    parser.add_argument("--min-coverage", type=int, default=100,
                        help="Minimum coverage for plotting data points (default: 100).")
    args = parser.parse_args()

    plot_qscore_calibration(args.summary_phred_csv, args.output_file,
                            log_scale=args.log, min_coverage=args.min_coverage)
