"""Summarize versioned RefSeq annotations for the PRJNA694496 DESeq2 result."""

import argparse
import csv
import math
from collections import Counter
from pathlib import Path


def log_choose(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeometric_upper_tail(population, successes, draws, observed):
    log_denominator = log_choose(population, draws)
    log_probabilities = [
        log_choose(successes, value) + log_choose(population - successes, draws - value)
        - log_denominator
        for value in range(observed, min(successes, draws) + 1)
        if 0 <= draws - value <= population - successes
    ]
    if not log_probabilities:
        return 0.0
    maximum = max(log_probabilities)
    return min(1.0, math.exp(maximum) * sum(math.exp(value - maximum) for value in log_probabilities))


def benjamini_hochberg(p_values):
    order = sorted(range(len(p_values)), key=p_values.__getitem__)
    adjusted = [1.0] * len(p_values)
    running = 1.0
    for rank_index in range(len(order) - 1, -1, -1):
        index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[index] * len(p_values) / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def summarize(results_path, annotations_path, output_dir):
    results = list(csv.DictReader(results_path.open(newline=""), delimiter="\t"))
    annotations = {
        row["feature_id_standardized"]: row
        for row in csv.DictReader(annotations_path.open(newline=""), delimiter="\t")
    }
    if len(annotations) != len(results):
        raise ValueError("Annotation table must have exactly one row per tested gene")
    groups = {"tested": [], "selected_higher_fdr_lt_0.05": [], "selected_lower_fdr_lt_0.05": []}
    for result in results:
        annotation = annotations[result["feature_id_standardized"]]
        merged = dict(result)
        merged.update(annotation)
        groups["tested"].append(merged)
        if result["padj"] and float(result["padj"]) < 0.05:
            group = "selected_higher_fdr_lt_0.05" if float(result["log2FoldChange"]) > 0 else "selected_lower_fdr_lt_0.05"
            groups[group].append(merged)

    output_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    for group, rows in groups.items():
        counts = Counter(row["gene_biotype"] or "missing" for row in rows)
        for biotype, count in sorted(counts.items()):
            summary_rows.append({
                "gene_set": group, "gene_biotype": biotype, "count": count,
                "fraction": count / len(rows), "gene_set_size": len(rows),
            })
    with (output_dir / "refseq_biotype_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    population = len(groups["tested"])
    background_counts = Counter(row["gene_biotype"] or "missing" for row in groups["tested"])
    enrichment_rows = []
    for group in ("selected_higher_fdr_lt_0.05", "selected_lower_fdr_lt_0.05"):
        group_counts = Counter(row["gene_biotype"] or "missing" for row in groups[group])
        draws = len(groups[group])
        for biotype, background_count in sorted(background_counts.items()):
            observed = group_counts[biotype]
            enrichment_rows.append({
                "gene_set": group,
                "gene_biotype": biotype,
                "observed": observed,
                "expected": draws * background_count / population,
                "fold_enrichment": (observed / draws) / (background_count / population),
                "p_value_upper_tail": hypergeometric_upper_tail(
                    population, background_count, draws, observed
                ),
            })
    adjusted = benjamini_hochberg([row["p_value_upper_tail"] for row in enrichment_rows])
    for row, value in zip(enrichment_rows, adjusted):
        row["adjusted_p_value_bh"] = value
    with (output_dir / "refseq_biotype_enrichment.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=enrichment_rows[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(enrichment_rows)

    lines = [
        "# PRJNA694496 RefSeq functional annotations", "",
        "Descriptions and gene biotypes come directly from GCF_963853765.1 / RS_2024_06. This is annotation and ranking, not ontology or pathway enrichment.", "",
    ]
    for group in ("selected_higher_fdr_lt_0.05", "selected_lower_fdr_lt_0.05"):
        rows = groups[group]
        characterized = [r for r in rows if r["description"] and not r["description"].lower().startswith("uncharacterized")]
        lines.extend([
            f"## {group.replace('_', ' ').title()}", "",
            f"Genes: {len(rows):,}; RefSeq-characterized descriptions: {len(characterized):,} ({100 * len(characterized) / len(rows):.1f}%).", "",
        ])
        ranked = sorted(characterized, key=lambda r: float(r["padj"]))[:15]
        for row in ranked:
            lines.append(
                f"- {row['gene_symbol']} ({row['feature_id_standardized']}): {row['description']} "
                f"(log2FC={float(row['log2FoldChange']):.2f}; adjusted p={float(row['padj']):.3g})"
            )
        lines.append("")
    lines.extend([
        "## Gene-biotype over-representation", "",
        "One-sided hypergeometric tests compare each significant direction with all genes having DESeq2 results; Benjamini-Hochberg adjustment covers both directions and all biotypes.", "",
    ])
    for row in enrichment_rows:
        if row["adjusted_p_value_bh"] < 0.05:
            lines.append(
                f"- {row['gene_set'].replace('_', ' ')}: {row['gene_biotype']} "
                f"({row['observed']} observed, {row['expected']:.1f} expected; "
                f"fold={row['fold_enrichment']:.2f}; adjusted p={row['adjusted_p_value_bh']:.3g})"
            )
    lines.extend(["",
        "Many oyster loci retain automated or uncharacterized RefSeq labels. Functional claims should be based on formal ontology/pathway enrichment or orthology-aware analysis rather than names alone.", "",
    ])
    (output_dir / "refseq_annotation_report.md").write_text("\n".join(lines))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--annotations", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summarize(args.results, args.annotations, args.output_dir)
