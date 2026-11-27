"""Generate deterministic synthetic quantitative proteomics LFQ matrix and design.csv."""

import argparse
import csv
import random
from pathlib import Path


GENE_ID_BASE = 105500000


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--proteins", type=int, default=150)
    args = parser.parse_args()

    rng = random.Random(20261128)
    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    samples = []
    design = []
    for cond in ("control", "heat"):
        for rep in "ABC":
            sample = "{}_{}".format(cond, rep)
            samples.append((sample, cond))
            design.append({"sample": sample, "condition": cond, "replicate": rep})

    header = ["feature_id", "feature_id_standardized", "gene_symbol", "description"] + [s for s, _ in samples]
    n_up = args.proteins // 10
    n_down = args.proteins // 10

    rows = []
    for idx in range(args.proteins):
        gid = GENE_ID_BASE + idx
        prot_id = "XP_{:09d}.1".format(11000000 + idx)
        base_log2 = rng.uniform(16.0, 24.0)
        row = {
            "feature_id": prot_id,
            "feature_id_standardized": "NCBI:GeneID:{}".format(gid),
            "gene_symbol": "HSP{}".format(idx) if idx < n_up else "LOC{}".format(gid),
            "description": "Heat shock cognate candidate {}".format(idx) if idx < n_up else "Synthetic bivalve protein {}".format(idx),
        }
        for sample, cond in samples:
            shift = 0.0
            if cond == "heat":
                if idx < n_up:
                    shift = 2.0
                elif idx < n_up + n_down:
                    shift = -2.0
            intensity = 2.0 ** max(1.0, base_log2 + shift + rng.gauss(0, 0.25))
            row[sample] = "{:.2f}".format(intensity)
        rows.append(row)

    with (out / "proteomics_abundance.tsv").open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=header, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    with (out / "design.csv").open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=["sample", "condition", "replicate"])
        w.writeheader()
        w.writerows(design)

    print("wrote synthetic proteomics dataset to {}".format(out))


if __name__ == "__main__":
    main()
