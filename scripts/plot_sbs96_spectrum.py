import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

dna_bases = ['A', 'C', 'G', 'T']

def plot_sbs96_spectrum(data_file, output_file):
    # Load the long-format SBS96 spectrum data
    # Columns: operation, prev_base, next_base, total, forward
    df = pd.read_csv(data_file)

    # Keep only substitution operations (exclude insertions/deletions)
    df = df[~df["operation"].str.startswith("->") & ~df["operation"].str.endswith(">-")].copy()

    ops_present = df["operation"].unique()

    # Canonical SBS96: only pyrimidine-context mutations C>* and T>*
    # Non-canonical: all 12 substitution types
    non_canonical_ops = {'A>C', 'A>G', 'A>T', 'G>A', 'G>C', 'G>T'}
    canonical = not any(op in ops_present for op in non_canonical_ops)

    if canonical:
        mutations = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
    else:
        mutations = ['C>T', 'G>A', 'G>T', 'G>C', 'C>A', 'T>A',
                     'T>C', 'A>G', 'T>G', 'C>G', 'A>C', 'A>T']

    def get_values_for_mut(mut):
        vals = []
        for first_base in dna_bases:
            for last_base in dna_bases:
                mask = (
                    (df["operation"] == mut) &
                    (df["prev_base"] == first_base) &
                    (df["next_base"] == last_base)
                )
                vals.append(float(df.loc[mask, "total"].sum()))
        return vals

    max_freq = 0
    for mut in mutations:
        values = get_values_for_mut(mut)
        max_freq = max(max_freq, max(values))

    bar_width = 0.1
    group_width = bar_width * 16
    fig_size = (6, 3) if canonical else (6, 5)

    if canonical:
        fig, ax = plt.subplots(figsize=fig_size)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        for i, mut in enumerate(mutations):
            values = get_values_for_mut(mut)
            x = [i * group_width + j * bar_width for j in range(16)]
            ax.bar(x, values, width=bar_width, color=plt.cm.tab10(i), label=mut)

        ax.set_xticks([i * group_width + 8 * bar_width for i in range(len(mutations))])
        ax.set_xticklabels(mutations)
        ax.set_ylabel('Frequency')
        ax.legend()

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close()

    else:
        top_muts = mutations[:6]
        bot_muts = mutations[6:]

        fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=fig_size)

        ax_top.spines['top'].set_visible(False)
        ax_top.spines['right'].set_visible(False)
        ax_bot.spines['top'].set_visible(True)
        ax_bot.spines['bottom'].set_visible(False)
        ax_bot.spines['right'].set_visible(False)

        for i, mut in enumerate(top_muts):
            values = get_values_for_mut(mut)
            x = [i * group_width + j * bar_width for j in range(16)]
            ax_top.bar(x, values, width=bar_width, color=plt.cm.tab10(i))

        ax_top.set_ylim(0, max_freq * 1.1)

        for i, mut in enumerate(bot_muts):
            values = get_values_for_mut(mut)
            x = [i * group_width + j * bar_width for j in range(16)]
            ax_bot.bar(x, [-v for v in values], width=bar_width, color=plt.cm.tab10(i))

        ax_bot.set_ylim(-max_freq * 1.1, 0)

        tick_positions = [i * group_width + 8 * bar_width for i in range(6)]
        ax_top.set_xticks(tick_positions)
        ax_top.set_xticklabels(top_muts)
        ax_bot.set_xticks(tick_positions)
        ax_bot.set_xticklabels(bot_muts)

        ax_top.tick_params(axis="x", top=False, labeltop=True, bottom=False, labelbottom=False)
        ax_bot.tick_params(axis="x", top=False, labeltop=False, bottom=False, labelbottom=True)

        from matplotlib.ticker import FuncFormatter
        formatter = FuncFormatter(lambda y, _: f"{abs(y):.2g}")
        ax_top.yaxis.set_major_formatter(formatter)
        ax_bot.yaxis.set_major_formatter(formatter)

        fig.text(0.02, 0.5, "Frequency", va="center", rotation="vertical")
        plt.tight_layout()
        fig.subplots_adjust(left=0.16)
        plt.savefig(output_file)
        plt.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot the SBS spectrum from Skiver summary_error_spectrum.csv.")
    parser.add_argument("summary_error_spectrum_csv", help="Path to the summary_error_spectrum CSV file.")
    parser.add_argument("output_file", help="Path to save the output plot image.")

    args = parser.parse_args()
    plot_sbs96_spectrum(args.summary_error_spectrum_csv, args.output_file)
