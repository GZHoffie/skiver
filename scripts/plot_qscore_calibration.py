import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def plot_qscore_calibration(calibration_file, output_file, log_scale=False, min_coverage=100):
    color_sub = 'steelblue'
    color_ins = 'darkorange'
    color_del = 'forestgreen'
    color_empirical = 'slategray'
    color_theoretical = 'black'

    df = pd.read_csv(calibration_file)
    df = df[df["num_observed"] >= min_coverage].copy()

    q_vals = df["qscore"].values
    sub_rate = df["normalized_substitution_rate"].values
    ins_rate = df["normalized_insertion_rate"].values
    del_rate = df["normalized_deletion_rate"].values
    err_rate = df["normalized_error_rate"].values
    counts   = df["num_observed"].values

    q_theory = np.linspace(q_vals.min(), q_vals.max(), 300)
    theory_rate = 10 ** (-q_theory / 10.0)

    fig, (ax_sub, ax_err, ax_hist) = plt.subplots(
        3, 1,
        figsize=(8, 9),
        gridspec_kw={"height_ratios": [2, 2, 1]},
        sharex=True,
    )

    # ── Subplot 1: substitution / insertion / deletion rates ──
    ax_sub.plot(q_vals, sub_rate, color=color_sub,
                linewidth=2, marker='o', markersize=3, label="Substitution")
    ax_sub.plot(q_vals, ins_rate, color=color_ins,
                linewidth=2, marker='s', markersize=3, label="Insertion")
    ax_sub.plot(q_vals, del_rate, color=color_del,
                linewidth=2, marker='^', markersize=3, label="Deletion")
    if log_scale:
        ax_sub.set_yscale('log')
    ax_sub.set_ylabel("Normalized rate (log)" if log_scale else "Normalized rate")
    ax_sub.set_title("Error component rates by quality score")
    ax_sub.legend(frameon=False)
    ax_sub.spines['top'].set_visible(False)
    ax_sub.spines['right'].set_visible(False)

    # ── Subplot 2: theoretical vs. empirical error rate ───────
    ax_err.plot(q_theory, theory_rate, color=color_theoretical,
                linestyle='--', linewidth=2.5, label="Theoretical ($10^{-Q/10}$)", zorder=3)
    ax_err.plot(q_vals, err_rate, color=color_empirical,
                linewidth=2.5, marker='o', markersize=3,
                label="Normalized error rate", zorder=4)
    if log_scale:
        ax_err.set_yscale('log')
    ax_err.set_ylabel("Error rate (log)" if log_scale else "Error rate")
    ax_err.set_title("Quality-score calibration")
    ax_err.legend(frameon=False)
    ax_err.spines['top'].set_visible(False)
    ax_err.spines['right'].set_visible(False)

    # ── Subplot 3: histogram of observed bases ────────────────
    ax_hist.bar(q_vals, counts, width=0.8, color=color_empirical, alpha=0.7)
    ax_hist.set_xlabel("Phred quality score ($Q$)")
    ax_hist.set_ylabel("# observed bases")
    ax_hist.spines['top'].set_visible(False)
    ax_hist.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Plot error component rates and calibration by Phred quality score."
    )
    parser.add_argument("summary_phred_csv",
                        help="Path to the summary Phred CSV file.")
    parser.add_argument("output_file", help="Path to save the output plot image.")
    parser.add_argument("--log", action="store_true",
                        help="Use logarithmic scale for the y-axes (default: False).")
    parser.add_argument("--min-coverage", type=int, default=100,
                        help="Minimum number of observed bases for plotting (default: 100).")
    args = parser.parse_args()

    plot_qscore_calibration(args.summary_phred_csv, args.output_file,
                            log_scale=args.log, min_coverage=args.min_coverage)
