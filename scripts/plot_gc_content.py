import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def plot_gc_content(gc_content_file, output_file, log_scale=False, min_bases=500):
    color_line = 'slategray'

    df = pd.read_csv(gc_content_file)

    df = df[(df["num_correct"] + df["num_error"]) > 0].copy()
    df["num_total"] = df["num_correct"] + df["num_error"]
    df["gc_content_mid"] = (df["gc_content_min"] + df["gc_content_max_exclusive"] - 1) / 2.0

    df_hist = df.copy()
    df_line = df[df["num_total"] >= min_bases].copy()

    fig, (ax_main, ax_hist) = plt.subplots(
        2, 1,
        figsize=(8, 7),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True,
    )

    # ── Main panel ────────────────────────────────────────────
    ax_main.plot(
        df_line["gc_content_mid"], df_line["per_base_error_rate"],
        color=color_line, linewidth=3, marker='o',
        label="Empirical error rate",
    )

    if log_scale:
        ax_main.set_yscale('log')
        ax_main.set_ylabel("Error rate (log scale)")
    else:
        ax_main.set_ylabel("Error rate")
    ax_main.set_title("Error rate vs. GC content")
    ax_main.legend(frameon=False)
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)

    if len(df_line) > 0:
        half_bin = (df_line["gc_content_max_exclusive"] - df_line["gc_content_min"]).iloc[0] / 2.0
        ax_main.set_xlim(
            df_line["gc_content_mid"].min() - half_bin,
            df_line["gc_content_mid"].max() + half_bin,
        )

    # ── Histogram panel ───────────────────────────────────────
    bin_width = (df_hist["gc_content_max_exclusive"] - df_hist["gc_content_min"]).iloc[0] \
        if len(df_hist) > 0 else 1.0
    ax_hist.bar(
        df_hist["gc_content_mid"], df_hist["num_total"],
        width=bin_width * 0.9, color=color_line, alpha=0.7,
    )
    ax_hist.axhline(min_bases, color='indianred', linestyle='--', linewidth=1,
                    label=f"min_bases = {min_bases}")
    ax_hist.set_xlabel("GC content (%)")
    ax_hist.set_ylabel("# bases")
    ax_hist.legend(frameon=False, fontsize=8)
    ax_hist.spines['top'].set_visible(False)
    ax_hist.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Plot per-base error rate by GC content."
    )
    parser.add_argument("summary_gc_content_csv",
                        help="Path to the summary GC content CSV file.")
    parser.add_argument("output_file", help="Path to save the output plot image.")
    parser.add_argument("--log", action="store_true",
                        help="Use logarithmic scale for the y-axis (default: False).")
    parser.add_argument("--min-bases", type=int, default=500,
                        help="Minimum number of bases for a GC bin to appear in the line plot (default: 500).")
    args = parser.parse_args()

    plot_gc_content(args.summary_gc_content_csv, args.output_file,
                    log_scale=args.log, min_bases=args.min_bases)
