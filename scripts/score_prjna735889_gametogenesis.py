#!/usr/bin/env python3
"""Call spermatogenesis status for the PRJNA735889 juveniles from the unadjusted DESeq2 run.

PC1 of the unadjusted run (about half the variance) is not pH but a male gametogenesis
programme: axonemal dyneins, calaxin, ropporin, synaptonemal-complex genes and cyclin B3.
The module is the expressed genes (mean log2 normalized count > 3) whose expression
correlates with PC1 at r < -0.9; each oyster's score is its mean centred log2 expression
over the module. Scores are bimodal, and oysters above the widest gap are called
``active``. The call is a transcriptional proxy for maturation, not a histological sex call.

Writes ``gametogenesis_scores.tsv`` and a copy of the design sheet with a
``gametogenesis`` column for ``~ gametogenesis + condition``.
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np


MODULE_CORRELATION = -0.9
MIN_MEAN_LOG2 = 3.0


def module_scores(normalized_counts, pca):
    with open(normalized_counts, newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    with open(pca, newline="") as handle:
        pc1 = {row["name"]: float(row["PC1"]) for row in csv.DictReader(handle)}
    samples = [name for name in rows[0] if name in pc1]
    columns = [rows[0].index(name) for name in samples]
    expression = np.log2(np.array([[float(row[i]) for i in columns] for row in rows[1:]]) + 1)
    expression = expression[expression.mean(1) > MIN_MEAN_LOG2]
    centred = expression - expression.mean(1, keepdims=True)
    axis = np.array([pc1[name] for name in samples])
    axis = axis - axis.mean()
    correlation = centred @ axis / np.sqrt((centred ** 2).sum(1) * (axis ** 2).sum())
    module = correlation < MODULE_CORRELATION
    # PC1's sign is arbitrary; orient the score so the module is high in active oysters.
    return samples, centred[module].mean(0), int(module.sum())


def widest_gap_threshold(scores):
    ordered = np.sort(scores)
    gaps = np.diff(ordered)
    return float(ordered[gaps.argmax()] + gaps.max() / 2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unadjusted", type=Path, required=True, help="unadjusted DESeq2 results directory")
    parser.add_argument("--design-sheet", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    samples, scores, size = module_scores(
        args.unadjusted / "acidified_vs_control_normalized_counts.tsv", args.unadjusted / "sample_pca.csv",
    )
    threshold = widest_gap_threshold(scores)
    status = {name: "active" if score > threshold else "inactive" for name, score in zip(samples, scores)}
    with (args.output_dir / "gametogenesis_scores.tsv").open("w") as handle:
        handle.write("sample\tspermatogenesis_module_score\tgametogenesis\n")
        for name, score in sorted(zip(samples, scores)):
            handle.write("{}\t{:.3f}\t{}\n".format(name, score, status[name]))
    with args.design_sheet.open(newline="") as handle:
        design = list(csv.DictReader(handle))
    output = args.output_dir / "deseq2_samplesheet_gametogenesis.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[*design[0], "gametogenesis"], lineterminator="\n")
        writer.writeheader()
        for row in design:
            writer.writerow({**row, "gametogenesis": status[row["sample"]]})
    for condition in ("control", "acidified"):
        called = [status[row["sample"]] for row in design if row["condition"] == condition]
        print("{}: {} of {} active".format(condition, called.count("active"), len(called)))
    print("module genes: {}; threshold {:.2f}; wrote {}".format(size, threshold, output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
