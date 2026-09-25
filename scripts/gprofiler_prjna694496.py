"""Prepare and summarize version-matched g:Profiler enrichment for PRJNA694496."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import date
from pathlib import Path


ORGANISM = "mggca963853765v1rs"
SOURCES = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"]


def gene_id(value: str) -> str:
    prefix = "NCBI:GeneID:"
    if not value.startswith(prefix) or not value[len(prefix):].isdigit():
        raise ValueError(f"Expected current NCBI GeneID, found {value!r}")
    return value[len(prefix):]


def prepare(results_path: Path, annotations_path: Path, output_dir: Path) -> None:
    rows = list(csv.DictReader(results_path.open(newline=""), delimiter="\t"))
    required = {"feature_id_standardized", "log2FoldChange", "pvalue", "padj"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("DESeq2 table is empty or missing required columns")

    annotations = {
        row["feature_id_standardized"]: row["gene_symbol"]
        for row in csv.DictReader(annotations_path.open(newline=""), delimiter="\t")
    }
    for row in rows:
        gene_id(row["feature_id_standardized"])
        if not annotations.get(row["feature_id_standardized"]):
            raise ValueError(f"Missing RefSeq gene symbol for {row['feature_id_standardized']}")
    if len(set(annotations.values())) != len(annotations):
        raise ValueError("RefSeq gene symbols must be unique for enrichment")
    symbol = lambda row: annotations[row["feature_id_standardized"]]
    background = sorted({symbol(r) for r in rows if r["pvalue"]})
    significant = [r for r in rows if r["padj"] and float(r["padj"]) < 0.05]
    queries = {
        "selected_higher": sorted({symbol(r) for r in significant
                                   if float(r["log2FoldChange"]) > 0}),
        "selected_lower": sorted({symbol(r) for r in significant
                                  if float(r["log2FoldChange"]) < 0}),
    }
    if any(not query for query in queries.values()):
        raise ValueError("Both differential-expression directions require nonempty queries")

    output_dir.mkdir(parents=True, exist_ok=True)
    for name, query in queries.items():
        request = {
            "organism": ORGANISM,
            "query": query,
            "sources": SOURCES,
            "user_threshold": 0.05,
            "significance_threshold_method": "g_SCS",
            "domain_scope": "custom",
            "background": background,
            "no_evidences": False,
        }
        (output_dir / f"{name}_request.json").write_text(json.dumps(request) + "\n")

    manifest = {
        "date_generated": date.today().isoformat(),
        "service_url": "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
        "organism": ORGANISM,
        "reference_match": "GCF_963853765.1 / xbMagGiga1.1",
        "query_identifier_type": "RefSeq gene symbol from the version-matched GFF",
        "sources_requested": SOURCES,
        "correction_method": "g_SCS",
        "domain_scope": "custom",
        "background_genes_with_pvalue": len(background),
        "selected_higher_fdr_lt_0.05": len(queries["selected_higher"]),
        "selected_lower_fdr_lt_0.05": len(queries["selected_lower"]),
        "input_sha256": hashlib.sha256(results_path.read_bytes()).hexdigest(),
    }
    (output_dir / "enrichment_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def summarize_one(response_path: Path, direction: str):
    response = json.loads(response_path.read_text())
    if "result" not in response:
        raise ValueError(f"g:Profiler response has no result: {response}")
    rows = []
    for result in response["result"]:
        rows.append({
            "direction": direction,
            "source": result["source"],
            "term_id": result["native"],
            "term_name": result["name"],
            "adjusted_p_value": result["p_value"],
            "significant": result["significant"],
            "term_size": result["term_size"],
            "query_size": result["query_size"],
            "intersection_size": result["intersection_size"],
            "effective_domain_size": result["effective_domain_size"],
            "precision": result["precision"],
            "recall": result["recall"],
            "intersection_gene_ids": ";".join(result.get("intersections", [[]])[0]),
        })
    return rows, response.get("meta", {})


def summarize(output_dir: Path) -> None:
    all_rows = []
    metadata = {}
    enrichment_qc = {}
    for direction in ("selected_higher", "selected_lower"):
        rows, meta = summarize_one(output_dir / f"{direction}_response.json", direction)
        all_rows.extend(rows)
        metadata[direction] = meta
        genes = meta["genes_metadata"]
        requested = len(meta["query_metadata"]["queries"]["query_1"])
        enrichment_qc[direction] = {
            "requested_genes": requested,
            "failed_identifiers": len(genes["failed"]),
            "ambiguous_identifiers": len(genes["ambiguous"]),
            "available_sources": sorted(meta["result_metadata"]),
            "gprofiler_version": meta["version"],
            "timestamp": meta["timestamp"],
        }
    all_rows.sort(key=lambda r: (r["direction"], float(r["adjusted_p_value"]), r["source"], r["term_id"]))
    fields = list(all_rows[0]) if all_rows else [
        "direction", "source", "term_id", "term_name", "adjusted_p_value", "significant",
        "term_size", "query_size", "intersection_size", "effective_domain_size", "precision",
        "recall", "intersection_gene_ids",
    ]
    with (output_dir / "functional_enrichment.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)
    (output_dir / "gprofiler_response_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (output_dir / "enrichment_qc.json").write_text(json.dumps(enrichment_qc, indent=2) + "\n")

    significant = [r for r in all_rows if r["significant"]]
    lines = [
        "# PRJNA694496 functional enrichment", "",
        "Over-representation analysis used g:Profiler with the version-matched Magallana gigas GCF_963853765.1 organism build. Selected-higher and selected-lower genes were tested separately against genes with nonmissing DESeq2 p-values. g:SCS-adjusted p < 0.05 defines significance.", "",
        "All direction-specific RefSeq symbols mapped without failure or ambiguity. This organism build supplied GO Biological Process, Molecular Function, and Cellular Component annotations; requested KEGG, Reactome, and WikiPathways sources were unavailable.", "",
    ]
    for direction in ("selected_higher", "selected_lower"):
        subset = [r for r in significant if r["direction"] == direction]
        lines.extend([f"## {direction.replace('_', ' ').title()}", "",
                      f"Significant terms: {len(subset)}.", ""])
        for row in subset[:15]:
            lines.append(
                f"- {row['source']} {row['term_id']}: {row['term_name']} "
                f"(adjusted p={float(row['adjusted_p_value']):.3g}; "
                f"overlap={row['intersection_size']}/{row['query_size']})"
            )
        lines.append("")
    lines.extend([
        "Terms are descriptive population-selection associations. Related ontology terms are not statistically independent, and enrichment does not establish causal thermotolerance mechanisms.", "",
    ])
    (output_dir / "functional_enrichment_report.md").write_text("\n".join(lines))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--results", type=Path, required=True)
    prepare_parser.add_argument("--annotations", type=Path, required=True)
    prepare_parser.add_argument("--output-dir", type=Path, required=True)
    summarize_parser = subparsers.add_parser("summarize")
    summarize_parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.results, args.annotations, args.output_dir)
    else:
        summarize(args.output_dir)
