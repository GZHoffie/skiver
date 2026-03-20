"""
plot_all.py – Generate all skiver plots from a set of output CSV file prefixes.

Usage:
    python plot_all.py <input_prefix> [<input_prefix> ...] -o <output_prefix> [options]

Each <input_prefix> corresponds to a skiver run whose output files follow the
naming convention:
    <prefix>.summary_error_rate.csv
    <prefix>.kvmer.csv
    <prefix>.hazard_rate.csv
    <prefix>.summary_error_spectrum.csv
    <prefix>.summary_error_spectrum_dependence_on_t.csv
    <prefix>.summary_phred.csv
    <prefix>.summary_read_position.csv

When a single input prefix is given the output plots are named:
    <output_prefix>_spectrum.png
    <output_prefix>_coverage.png
    ...

When multiple input prefixes are given each is distinguished by a numeric suffix
derived from its position in the argument list (0-indexed):
    <output_prefix>_0_spectrum.png
    <output_prefix>_1_spectrum.png
    ...
"""

import os
import sys
import argparse

# ---------------------------------------------------------------------------
# Import individual plot functions from sibling scripts.
# ---------------------------------------------------------------------------
_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _script_dir)

from plot_spectrum import plot_spectrum
from plot_coverage import plot_coverage_histogram
from plot_hazard_survival_rate import plot_hazard_survival_rate
from plot_qscore_calibration import plot_qscore_calibration
from plot_read_position import plot_read_position
from plot_sbs96_spectrum import plot_sbs96_spectrum
from plot_error_spectrum_dependence_on_t import plot_error_spectrum_dependence_on_t


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _exists(path):
    if os.path.isfile(path):
        return True
    print(f"  [skip] file not found: {path}", file=sys.stderr)
    return False


def generate_plots(input_prefix, output_prefix,
                   normalize=False, log_scale=False,
                   t_min=1, t_max=100,
                   num_bases=100, min_coverage=100):
    """Generate all available plots for a single *input_prefix*."""

    report   = f"{input_prefix}.summary_error_rate.csv"
    kvmer    = f"{input_prefix}.kvmer.csv"
    hazard   = f"{input_prefix}.hazard_rate.csv"
    spectrum = f"{input_prefix}.summary_error_spectrum.csv"
    dep_t    = f"{input_prefix}.summary_error_spectrum_dependence_on_t.csv"
    phred    = f"{input_prefix}.summary_phred.csv"
    readpos  = f"{input_prefix}.summary_read_position.csv"

    plots = [
        # (description, required_files, thunk)
        (
            "error spectrum",
            [spectrum, report],
            lambda: plot_spectrum(spectrum, report,
                                  f"{output_prefix}_spectrum.png",
                                  normalize=normalize),
        ),
        (
            "coverage histogram",
            [kvmer, report],
            lambda: plot_coverage_histogram(kvmer, report,
                                            f"{output_prefix}_coverage.png"),
        ),
        (
            "hazard/survival rate",
            [hazard, report],
            lambda: plot_hazard_survival_rate(hazard, report,
                                              f"{output_prefix}_hazard_survival.png",
                                              t_min=t_min, t_max=t_max,
                                              log_scale=log_scale),
        ),
        (
            "SBS-96 spectrum",
            [spectrum],
            lambda: plot_sbs96_spectrum(spectrum,
                                        f"{output_prefix}_sbs96_spectrum.png"),
        ),
        (
            "error spectrum dependence on t",
            [dep_t],
            lambda: plot_error_spectrum_dependence_on_t(dep_t,
                                                        f"{output_prefix}_error_spectrum_dep_t.png"),
        ),
        (
            "quality-score calibration",
            [phred],
            lambda: plot_qscore_calibration(phred,
                                            f"{output_prefix}_qscore_calibration.png",
                                            log_scale=log_scale,
                                            min_coverage=min_coverage),
        ),
        (
            "read position error rate",
            [readpos],
            lambda: plot_read_position(readpos,
                                       f"{output_prefix}_read_position.png",
                                       num_bases=num_bases),
        ),
    ]

    for description, required, thunk in plots:
        if all(_exists(f) for f in required):
            print(f"  Plotting {description}...")
            thunk()
        # missing-file warning already printed by _exists()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate all skiver plots from one or more sets of output CSV files. "
            "Each INPUT_PREFIX should match the -o/--output-prefix used when running skiver."
        )
    )
    parser.add_argument(
        "input_prefixes",
        nargs="+",
        metavar="INPUT_PREFIX",
        help="One or more skiver output prefixes (e.g. ./results/sample1).",
    )
    parser.add_argument(
        "-o", "--output-prefix",
        required=True,
        metavar="OUTPUT_PREFIX",
        help="Prefix for the generated plot files.",
    )

    # Forwarded to individual plot scripts
    parser.add_argument("--normalize", action="store_true",
                        help="Normalize the error spectrum (passed to plot_spectrum).")
    parser.add_argument("--log-scale", action="store_true",
                        help="Use log scale where applicable.")
    parser.add_argument("-t", type=int, default=1,
                        help="Minimum t for survival rate curve (default: 1).")
    parser.add_argument("-T", type=int, default=100,
                        help="Maximum t for survival rate curve (default: 100).")
    parser.add_argument("--num-bases", type=int, default=100,
                        help="Number of bases from each end to plot for read position (default: 100).")
    parser.add_argument("--min-coverage", type=int, default=100,
                        help="Minimum coverage for qscore calibration plot (default: 100).")

    args = parser.parse_args()

    multi = len(args.input_prefixes) > 1

    for idx, input_prefix in enumerate(args.input_prefixes):
        if multi:
            out_prefix = f"{args.output_prefix}_{idx}"
            print(f"[{idx}] input prefix: {input_prefix}  →  output prefix: {out_prefix}")
        else:
            out_prefix = args.output_prefix
            print(f"input prefix: {input_prefix}  →  output prefix: {out_prefix}")

        generate_plots(
            input_prefix=input_prefix,
            output_prefix=out_prefix,
            normalize=args.normalize,
            log_scale=args.log_scale,
            t_min=args.t,
            t_max=args.T,
            num_bases=args.num_bases,
            min_coverage=args.min_coverage
        )

    print("Done.")


if __name__ == "__main__":
    main()
