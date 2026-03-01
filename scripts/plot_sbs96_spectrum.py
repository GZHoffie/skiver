import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

dna_bases = ['A', 'C', 'G', 'T']

def plot_sbs96_spectrum(data_file, output_file):
    # Load the SBS96 spectrum data
    df = pd.read_csv(data_file)

    # Check if we include canonical substitutions.
    canonical = True
    if 'A[A>C]A' in df.columns:
        canonical = False

    if canonical:
        mutations = ['C>A', 'C>G', 'C>T',
                     'T>A', 'T>C', 'T>G']
    else:
        mutations = ['C>T', 'G>A', 'G>T', 'G>C', 'C>A', "T>A",
                     'T>C', 'A>G', 'T>G', 'C>G', 'A>C', 'A>T']
        
    def get_values_for_mut(mut: str):
        vals = []
        for first_base in dna_bases:
            for last_base in dna_bases:
                context = f"{first_base}[{mut}]{last_base}"
                vals.append(float(df.at[0, context]))
        return vals

    # Get the maximum frequency for scaling y-axis
    max_freq = 0
    total_freq = 0
    for mut in mutations:
        values = get_values_for_mut(mut)
        max_freq = max(max_freq, max(values))
        total_freq += sum(values)



    # Plotting parameters (same style as canonical)
    bar_width = 0.1
    group_width = bar_width * 16
    fig_size = (6, 3) if canonical else (6, 5)

    

    if canonical:
        fig, ax = plt.subplots(figsize=fig_size)

        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        for i, mut in enumerate(mutations):
            values = get_values_for_mut(mut)
            values = [v for v in values]
            x = [i * group_width + j * bar_width for j in range(16)]
            ax.bar(x, values, width=bar_width, color=plt.cm.tab10(i), label=mut)

        ax.set_xticks([i * group_width + 8 * bar_width for i in range(len(mutations))])
        ax.set_xticklabels(mutations)
        ax.set_ylabel('Frequency')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.show()
        plt.close()

    else:
        top_muts = mutations[:6]
        bot_muts = mutations[6:]

        fig, (ax_top, ax_bot) = plt.subplots(
            2, 1,
            figsize=fig_size,
            #sharex=True,
            #gridspec_kw={"height_ratios": [1, 1], "hspace": 0.05}
        )

        ax_top.spines['top'].set_visible(False)
        ax_top.spines['right'].set_visible(False)
        ax_bot.spines['bottom'].set_visible(False)
        ax_bot.spines['right'].set_visible(False)

        # Top bars
        for i, mut in enumerate(top_muts):
            values = get_values_for_mut(mut)
            values = [v for v in values]
            x = [i * group_width + j * bar_width for j in range(16)]
            ax_top.bar(x, values, width=bar_width, color=plt.cm.tab10(i))
        #print(values)

        ax_top.set_ylim(0, max_freq * 1.1)

        # Bottom bars
        for i, mut in enumerate(bot_muts):
            values = get_values_for_mut(mut)
            values = [v for v in values]
            x = [i * group_width + j * bar_width for j in range(16)]
            ax_bot.bar(x, [-v for v in values], width=bar_width, color=plt.cm.tab10(i))

        ax_bot.set_ylim(-max_freq * 1.1, 0)
        

        # ticks
        tick_positions = [i * group_width + 8 * bar_width for i in range(6)]
        ax_top.set_xticks(tick_positions)
        ax_top.set_xticklabels(top_muts)
        ax_bot.set_xticks(tick_positions)
        ax_bot.set_xticklabels(bot_muts)

        ax_top.tick_params(
            axis="x",
            top=False,
            labeltop=True,
            bottom=False,
            labelbottom=False
        )

        ax_bot.tick_params(
            axis="x",
            top=False,
            labeltop=False,
            bottom=False,
            labelbottom=True
        )

        from matplotlib.ticker import MaxNLocator, FuncFormatter

        formatter = FuncFormatter(lambda y, _: f"{abs(y):.2g}")
        ax_top.yaxis.set_major_formatter(formatter)
        ax_bot.yaxis.set_major_formatter(formatter)

        # shared label
        fig.text(0.02, 0.5, "Frequency", va="center", rotation="vertical")
        plt.tight_layout()
        fig.subplots_adjust(left=0.16)
        plt.savefig(output_file)
        plt.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot the SBS spectrum from Skiver report.")
    parser.add_argument("skiver_report_file", help="Path to the Skiver report CSV file.")
    parser.add_argument("output_file", help="Path to save the output plot image.")

    args = parser.parse_args()
    plot_sbs96_spectrum(args.skiver_report_file, args.output_file)
    