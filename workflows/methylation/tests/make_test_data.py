"""Generate deterministic synthetic WGBS CpG methylation coverage matrix and design.csv."""

import argparse
import csv
import random
from pathlib import Path


GENE_ID_BASE = 105500000


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sites", type=int, default=200)
    args = parser.parse_args()

    rng = random.Random(20261127)
    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    samples = []
    design = []
    for cond in ("control", "heat"):
        for rep in "ABC":
            sample = "{}_{}".format(cond, rep)
            samples.append((sample, cond))
            design.append({"sample": sample, "condition": cond, "replicate": rep})

    header = ["feature_id", "feature_id_standardized", "gene_symbol", "description"]
    for sample, _ in samples:
        header.extend([f"{sample}_M", f"{sample}_U"])

    rows = []
    n_hyper = args.sites // 10
    n_hypo = args.sites // 10
    for idx in range(args.sites):
        gid = GENE_ID_BASE + idx
        cpg_id = f"NC_TEST:{1000 + idx * 50}"
        std_id = f"NCBI:GeneID:{gid}"
        base_beta = rng.uniform(0.25, 0.70)
        row = {
            "feature_id": cpg_id,
            "feature_id_standardized": std_id,
            "gene_symbol": f"LOC{gid}",
            "description": f"Synthetic promoter CpG site near LOC{gid}",
        }
        for sample, cond in samples:
            depth = rng.randint(25, 45)
            beta = base_beta
            if cond == "heat":
                if idx < n_hyper:
                    beta = min(0.95, base_beta + 0.35)
                elif idx < n_hyper + n_hypo:
                    beta = max(0.05, base_beta - 0.35)
            beta = min(0.98, max(0.02, beta + rng.gauss(0, 0.03)))
            m_reads = int(round(depth * beta))
            u_reads = max(0, depth - m_reads)
            row[f"{sample}_M"] = m_reads
            row[f"{sample}_U"] = u_reads
        rows.append(row)

    with (out / "methylation_coverage.tsv").open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=header, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    with (out / "design.csv").open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=["sample", "condition", "replicate"])
        w.writeheader()
        w.writerows(design)

    print("wrote synthetic methylation dataset to {}".format(out))


if __name__ == "__main__":
    main()
