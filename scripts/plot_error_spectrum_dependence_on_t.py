import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import numpy as np


def classify_op(op):
    if ">" in op and not op.startswith("-") and not op.endswith("-"):
        return "Substitution"
    elif op.startswith("->"):
        return "Insertion"
    elif op.endswith(">-"):
        return "Deletion"
    return "Other"


def plot_error_spectrum_dependence_on_t(csv_file, output_file):
    df = pd.read_csv(csv_file)

    pos_cols = [c for c in df.columns if c.startswith("freq_at_t")]
    positions = [int(c.split("t")[-1]) for c in pos_cols]

    df["type"] = df["operation"].apply(classify_op)

    # For each position column, compute column total and proportions
    col_totals = df[pos_cols].sum(axis=0)

    # Filter positions with nonzero total
    valid_mask = col_totals > 0
    valid_pos_cols = [c for c, v in zip(pos_cols, valid_mask) if v]
    valid_positions = [p for p, v in zip(positions, valid_mask) if v]
    valid_t = [p for p in valid_positions]
    valid_totals = col_totals[valid_pos_cols]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # --- Left subplot: proportion by error type (Substitution / Insertion / Deletion) ---
    ax1 = axes[0]
    type_colors = {"Substitution": "steelblue", "Insertion": "tomato", "Deletion": "seagreen"}
    grand_total = df["total"].sum()
    for etype, color in type_colors.items():
        mask = df["type"] == etype
        type_counts = df.loc[mask, valid_pos_cols].sum(axis=0)
        proportions = type_counts.values / valid_totals.values
        ax1.plot(valid_t, proportions, label=etype, color=color, linewidth=2, marker='o')
        if grand_total > 0:
            total_proportion = df.loc[mask, "total"].sum() / grand_total
            ax1.axhline(total_proportion, color=color, linewidth=1.2, linestyle="--")

    ax1.set_xlabel("$t$")
    ax1.set_ylabel("Proportion")
    ax1.set_title("Error type proportion by $t$")
    ax1.legend()
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # --- Right subplot: proportion of each edit operation ---
    ax2 = axes[1]
    unique_ops = df["operation"].unique()

    # Substitution ops get one colormap, insertions another, deletions another
    sub_ops = [o for o in unique_ops if classify_op(o) == "Substitution"]
    ins_ops = [o for o in unique_ops if classify_op(o) == "Insertion"]
    del_ops = [o for o in unique_ops if classify_op(o) == "Deletion"]

    def ops_with_colors(ops, cmap_name, n):
        cmap = plt.get_cmap(cmap_name, max(n, 1))
        return [(op, cmap(i)) for i, op in enumerate(ops)]

    colored_ops = (
        ops_with_colors(sub_ops, "Blues", len(sub_ops)) +
        ops_with_colors(ins_ops, "Reds", len(ins_ops)) +
        ops_with_colors(del_ops, "Greens", len(del_ops))
    )

    for op, color in colored_ops:
        mask = df["operation"] == op
        op_counts = df.loc[mask, valid_pos_cols].sum(axis=0)
        proportions = op_counts.values / valid_totals.values
        ax2.plot(valid_t, proportions, label=op, color=color, linewidth=1.2)

    ax2.set_xlabel("$t$")
    ax2.set_ylabel("Proportion")
    ax2.set_title("Edit operation proportion by $t$")
    ax2.legend(fontsize=6, ncol=2)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)


    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Plot error spectrum dependence on value position."
    )
    parser.add_argument("summary_error_rate_dependence_csv",
                        help="Path to summary_error_spectrum_dependence_on_t.csv")
    parser.add_argument("output_file", help="Path to save the output plot image.")
    args = parser.parse_args()

    plot_error_spectrum_dependence_on_t(args.summary_error_rate_dependence_csv, args.output_file)
