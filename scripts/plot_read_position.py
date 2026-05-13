import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.ndimage import uniform_filter1d


def _wilson_interval(num_error, num_total, z=1.6448536269514722):
    """Return a two-sided 90% Wilson interval, matching 5th-95th bands."""
    n = np.asarray(num_total, dtype=float)
    valid = n > 0
    n_safe = np.where(valid, n, 1.0)
    p = np.divide(
        np.asarray(num_error, dtype=float),
        n,
        out=np.zeros_like(n, dtype=float),
        where=valid,
    )
    z2 = z ** 2
    denom = 1 + z2 / n_safe
    center = (p + z2 / (2 * n_safe)) / denom
    half_width = (
        z
        * np.sqrt((p * (1 - p) / n_safe) + (z2 / (4 * n_safe ** 2)))
        / denom
    )
    lower = np.clip(center - half_width, 0, 1)
    upper = np.clip(center + half_width, 0, 1)
    lower[~valid] = 0
    upper[~valid] = 0
    return lower, upper


def plot_read_position(read_position_file, output_file, num_bases=100):
    color = 'slategray'
    color_smooth = 'indianred'

    df = pd.read_csv(read_position_file)
    df["num_total"] = df["num_correct"] + df["num_error"]

    df_start = df[df["from_start"] == True].sort_values("index")
    df_end = df[df["from_start"] == False].sort_values("index")

    df_start = df_start[df_start["index"] <= num_bases]
    df_end = df_end[df_end["index"] <= num_bases]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8),
                             gridspec_kw={"height_ratios": [3, 1]})

    ci_start = _wilson_interval(df_start["num_error"], df_start["num_total"])
    ci_end = _wilson_interval(df_end["num_error"], df_end["num_total"])
    max_y = max(ci_start[1].max(), ci_end[1].max()) * 1.1

    for col, (df_pos, label, ci) in enumerate([
        (df_start, "from start", ci_start),
        (df_end, "from end", ci_end),
    ]):
        ax_line = axes[0][col]
        ax_bar = axes[1][col]

        x = df_pos["index"].values
        y = df_pos["error_rate"].values
        counts = df_pos["num_total"].values
        ci_lower, ci_upper = ci

        ax_line.fill_between(
            x, ci_lower, ci_upper,
            color=color, alpha=0.25, label="5%-95% Wilson interval",
        )

        # Raw error rate line
        ax_line.plot(x, y, color=color, linewidth=1.5, alpha=0.6, label="Error rate")

        # Smoothed line (uniform filter, window = 5% of num_bases, min 3)
        window = max(3, num_bases // 10)
        if len(y) >= window:
            y_smooth = uniform_filter1d(y, size=window)
            ax_line.plot(x, y_smooth, color=color_smooth, linewidth=2.5,
                         label=f"Smoothed (window={window})")

        ax_line.set_title(f"Error rate ({label})")
        ax_line.set_ylabel("Error rate")
        ax_line.set_xlim(x.min(), x.max())
        ax_line.set_ylim(0, max_y)
        ax_line.legend(frameon=False)
        ax_line.spines['top'].set_visible(False)
        ax_line.spines['right'].set_visible(False)

        # Number of bases bar
        ax_bar.bar(x, counts, width=1.0, color=color, alpha=0.7)
        ax_bar.set_xlabel("Position in read")
        ax_bar.set_ylabel("# bases used for estimation")
        ax_bar.set_xlim(x.min(), x.max())
        ax_bar.spines['top'].set_visible(False)
        ax_bar.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Plot error rate by read position from summary_read_position.csv."
    )
    parser.add_argument("summary_read_position_csv",
                        help="Path to the summary_read_position CSV file.")
    parser.add_argument("output_file", help="Path to save the output plot image.")
    parser.add_argument("--num-bases", type=int, default=100,
                        help="Number of bases from each end to plot (default: 100).")
    args = parser.parse_args()

    plot_read_position(args.summary_read_position_csv, args.output_file,
                       num_bases=args.num_bases)
