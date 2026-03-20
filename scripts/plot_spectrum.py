import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


COLS = ["A", "C", "G", "T", "-"]
BASES = ["A", "C", "G", "T"]


def _build_spectrum_matrix(df, count_col):
    """Build a 5×5 (from_base × to_base) spectrum matrix from the long-format df."""
    matrix = np.zeros((5, 5))
    for _, row in df.iterrows():
        op = row["operation"]          # e.g. "A>C", "->A", "A>-"
        from_base, to_base = op.split(">")
        fi = COLS.index(from_base) if from_base in COLS else COLS.index("-")
        ti = COLS.index(to_base)   if to_base   in COLS else COLS.index("-")
        matrix[fi, ti] += row[count_col]
    return matrix


def _error_type_spectrum(df, count_col):
    """Sum counts by error type from the long-format df."""
    spec = {"Insertion": 0.0, "Deletion": 0.0, "Substitution": 0.0}
    for _, row in df.iterrows():
        op = row["operation"]
        count = row[count_col]
        if op.startswith("->"):
            spec["Insertion"] += count
        elif op.endswith(">-"):
            spec["Deletion"] += count
        else:
            spec["Substitution"] += count
    return spec


def _plot_row(axes, df, count_col, row_title_suffix, target_sum, cmap, normalize=False):
    """Fill one row of the 2×2 figure (heatmap + bar chart).

    target_sum: the value that the spectrum entries should sum to after scaling.
                1.0 → relative proportions; per_base_error_rate → absolute rates.
    """
    ax_heat, ax_bar = axes

    # ── Spectrum matrix ────────────────────────────────────────────────────
    matrix = _build_spectrum_matrix(df, count_col)
    total = matrix.sum()
    if total > 0:
        matrix = matrix / total * target_sum * 100

    mask = np.eye(5, dtype=bool)
    sns.heatmap(matrix, annot=True, fmt=".3f",
                xticklabels=COLS, yticklabels=COLS,
                cbar=True, mask=mask,
                linewidths=3, linecolor='white',
                cmap=cmap, ax=ax_heat)

    normalized_suffix = " (normalized)" if normalize else ""
    ax_heat.set_title(f"Error spectrum{normalized_suffix} — {row_title_suffix}")
    ax_heat.set_ylabel("Original base")
    ax_heat.set_xlabel("Observed base")

    # ── Error type bar chart ───────────────────────────────────────────────
    spec = _error_type_spectrum(df, count_col)
    values = np.array([spec["Insertion"], spec["Deletion"], spec["Substitution"]], dtype=float)
    spec_total = values.sum()
    if spec_total > 0:
        values = values / spec_total * target_sum * 100

    x = range(3)
    bars = ax_bar.bar(x, values, color='slategray', width=0.6)
    ax_bar.bar_label(bars, fmt='%.3f', padding=3)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(["Insertion", "Deletion", "Substitution"])
    ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)
    ax_bar.set_ylim(0, max(values) * 1.15 if values.max() > 0 else 1)
    ax_bar.set_title(f"Error type distribution{normalized_suffix} — {row_title_suffix}")
    ax_bar.set_ylabel("Proportion (%)" if normalize else "Error rate (%)")
    ax_bar.set_xlabel("Error Type")


def plot_spectrum(skiver_spectrum_file, skiver_error_rate_file, output_file, normalize=False):
    """Plot the error spectrum as a 2×2 figure.

    Top row    – both strands combined (``total`` column).
    Bottom row – forward strand only   (``forward`` column).

    normalize=False (default): each panel is scaled so entries sum to per_base_error_rate from skiver_error_rate_file.
    normalize=True:            each panel is scaled so entries sum to 1.
    """
    color = 'slategray'
    cmap = sns.light_palette(color, as_cmap=True)

    df = pd.read_csv(skiver_spectrum_file)

    if not normalize:
        rate_df = pd.read_csv(skiver_error_rate_file)
        # [FIXME] use effective_error_rate instead? Which one is more appropriate?
        target_sum = rate_df["per_base_error_rate"].item()
    else:
        target_sum = 1.0

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    _plot_row(axes[0], df, "total",   "both strands",  target_sum, cmap, normalize)
    _plot_row(axes[1], df, "forward", "forward strand", target_sum, cmap, normalize)

    plt.subplots_adjust(wspace=0.3, hspace=0.45)
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Plot the error spectrum from Skiver summary_error_spectrum.csv."
    )
    parser.add_argument("summary_error_spectrum_csv",
                        help="Path to the summary_error_spectrum CSV file.")
    parser.add_argument("summary_error_rate_csv",
                        help="Path to the summary_error_rate CSV file.")
    parser.add_argument("output_file", help="Path to save the output plot image.")
    parser.add_argument("--normalize", action="store_true",
                        help="Scale entries so they sum to 1. If not set, the entries are scaled to sum to the per_base_error_rate from the summary_error_rate file.")
    args = parser.parse_args()

    plot_spectrum(args.summary_error_spectrum_csv, args.summary_error_rate_csv,
                  args.output_file, normalize=args.normalize)
