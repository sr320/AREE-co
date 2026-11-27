"""Generate deterministic synthetic metabolomics peak matrix and design.csv."""

import argparse
import csv
import random
from pathlib import Path


CHEBI_BASE = 15300


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metabolites", type=int, default=100)
    args = parser.parse_args()

    rng = random.Random(20261129)
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
    n_up = args.metabolites // 10
    n_down = args.metabolites // 10

    rows = []
    for idx in range(args.metabolites):
        cid = "CHEBI:{}".format(CHEBI_BASE + idx)
        base_log2 = rng.uniform(12.0, 19.0)
        row = {
            "feature_id": cid,
            "feature_id_standardized": cid,
            "gene_symbol": "taurine" if idx == 0 else "betaine" if idx == 1 else "MET_{:03d}".format(idx),
            "description": "Synthetic bivalve osmolyte/metabolite {}".format(idx),
        }
        for sample, cond in samples:
            shift = 0.0
            if cond == "heat":
                if idx < n_up:
                    shift = 1.8
                elif idx < n_up + n_down:
                    shift = -1.8
            val = 2.0 ** max(1.0, base_log2 + shift + rng.gauss(0, 0.22))
            row[sample] = "{:.2f}".format(val)
        rows.append(row)

    with (out / "metabolomics_abundance.tsv").open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=header, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    with (out / "design.csv").open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=["sample", "condition", "replicate"])
        w.writeheader()
        w.writerows(design)

    print("wrote synthetic metabolomics dataset to {}".format(out))


if __name__ == "__main__":
    main()
