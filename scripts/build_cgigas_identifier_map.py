"""Build a conservative oyster_v9 CGI-to-current-NCBI-Gene mapping release."""

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_attributes(text):
    attributes = {}
    for item in text.rstrip().split(";"):
        if not item:
            continue
        key, _, value = item.partition("=")
        attributes[key] = unquote(value)
    return attributes


def load_feature_ids(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return {row["feature_id_original"] for row in csv.DictReader(handle, delimiter="\t")}


def load_assembly_aliases(path):
    aliases = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 7 and fields[6] != "na":
                aliases[fields[0]] = fields[6]
    return aliases


def load_legacy_genes(path, selected_ids):
    genes = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = parse_attributes(fields[8])
            gene_id = attributes.get("gene_id")
            if gene_id not in selected_ids:
                continue
            genes[gene_id] = {
                "seqid": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "description": attributes.get("description", ""),
            }
    return genes


def load_ncbi_genes(path):
    by_seqid = defaultdict(list)
    by_gene_id = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = parse_attributes(fields[8])
            dbxref = attributes.get("Dbxref", "")
            gene_id = None
            for item in dbxref.split(","):
                if item.startswith("GeneID:"):
                    gene_id = item.split(":", 1)[1]
                    break
            if gene_id is None:
                continue
            gene = {
                "gene_id": gene_id,
                "seqid": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "symbol": attributes.get("Name", attributes.get("gene", "")),
                "description": attributes.get("description", ""),
            }
            by_seqid[fields[0]].append(gene)
            by_gene_id[gene_id] = gene
    for genes in by_seqid.values():
        genes.sort(key=lambda gene: gene["start"])
    return by_seqid, by_gene_id


def load_gene_replacements(path):
    if path is None:
        return {}
    with Path(path).open(encoding="utf-8") as handle:
        result = json.load(handle)["result"]
    replacements = {}
    for old_gene_id in result.get("uids", []):
        current_gene_id = result[old_gene_id].get("currentid")
        if current_gene_id:
            replacements[old_gene_id] = str(current_gene_id)
    return replacements


def overlap_metrics(source, target):
    overlap = max(0, min(source["end"], target["end"]) - max(source["start"], target["start"]) + 1)
    source_length = source["end"] - source["start"] + 1
    target_length = target["end"] - target["start"] + 1
    return overlap, overlap / source_length, overlap / target_length


def map_feature(
    feature_id,
    legacy,
    aliases,
    legacy_ncbi_by_seq,
    current_by_gene_id,
    gene_replacements,
    release_id,
):
    unresolved = {
        "feature_id_original": feature_id,
        "feature_id_standardized": feature_id,
        "ortholog_reference": "",
        "mapping_confidence": "unresolved",
        "mapping_release": release_id,
    }
    if legacy is None:
        return {
            **unresolved,
            "mapping_evidence": "legacy_novel_gene_no_genomic_model",
        }

    refseq = aliases.get(legacy["seqid"])
    if refseq is None:
        return {
            **unresolved,
            "mapping_evidence": "legacy_scaffold_has_no_refseq_alias",
            "legacy_scaffold": legacy["seqid"],
            "legacy_start": legacy["start"],
            "legacy_end": legacy["end"],
            "legacy_strand": legacy["strand"],
        }

    candidates = []
    for target in legacy_ncbi_by_seq.get(refseq, []):
        if target["start"] > legacy["end"]:
            break
        if target["end"] < legacy["start"] or target["strand"] != legacy["strand"]:
            continue
        overlap, source_coverage, target_coverage = overlap_metrics(legacy, target)
        if overlap:
            candidates.append((min(source_coverage, target_coverage), overlap, source_coverage, target_coverage, target))
    candidates.sort(key=lambda item: (item[0], item[1]), reverse=True)
    qualifying = [item for item in candidates if item[2] >= 0.8 or item[0] >= 0.5]

    common = {
        "legacy_scaffold": legacy["seqid"],
        "legacy_refseq_accession": refseq,
        "legacy_start": legacy["start"],
        "legacy_end": legacy["end"],
        "legacy_strand": legacy["strand"],
        "legacy_description": legacy["description"],
    }
    if not qualifying:
        return {
            **unresolved,
            **common,
            "mapping_evidence": "no_unique_same_strand_legacy_coverage_ge_0.8_or_reciprocal_ge_0.5",
            "candidate_overlap_count": len(candidates),
        }
    if len(qualifying) > 1:
        return {
            **unresolved,
            **common,
            "mapping_evidence": "ambiguous_multiple_same_strand_legacy_coverage_ge_0.8_or_reciprocal_ge_0.5",
            "candidate_overlap_count": len(qualifying),
        }

    _, overlap, source_coverage, target_coverage, target = qualifying[0]
    current_gene_id = target["gene_id"]
    current = current_by_gene_id.get(current_gene_id)
    if current is None:
        replacement_gene_id = gene_replacements.get(target["gene_id"])
        replacement = current_by_gene_id.get(replacement_gene_id)
        if replacement is not None:
            return {
                "feature_id_original": feature_id,
                "feature_id_standardized": "NCBI:GeneID:{}".format(replacement_gene_id),
                "ortholog_reference": "",
                "mapping_confidence": "inferred",
                "mapping_release": release_id,
                "mapping_evidence": "unique_coordinate_overlap_to_discontinued_ncbi_gene_then_ncbi_currentid",
                **common,
                "legacy_ncbi_gene_id": target["gene_id"],
                "ncbi_replacement_gene_id": replacement_gene_id,
                "legacy_ncbi_start": target["start"],
                "legacy_ncbi_end": target["end"],
                "overlap_bp": overlap,
                "legacy_coverage": round(source_coverage, 6),
                "legacy_ncbi_coverage": round(target_coverage, 6),
                "current_symbol": replacement["symbol"],
                "current_description": replacement["description"],
                "current_refseq_accession": replacement["seqid"],
                "current_start": replacement["start"],
                "current_end": replacement["end"],
                "current_strand": replacement["strand"],
                "candidate_overlap_count": 1,
            }
        return {
            **unresolved,
            **common,
            "mapping_evidence": (
                "ncbi_replacement_gene_absent_from_current_annotation"
                if replacement_gene_id
                else "legacy_ncbi_gene_discontinued_without_current_replacement"
            ),
            "legacy_ncbi_gene_id": target["gene_id"],
            "ncbi_replacement_gene_id": replacement_gene_id or "",
            "overlap_bp": overlap,
            "legacy_coverage": round(source_coverage, 6),
            "legacy_ncbi_coverage": round(target_coverage, 6),
        }

    exact_coordinates = (
        legacy["start"] == target["start"]
        and legacy["end"] == target["end"]
        and legacy["strand"] == target["strand"]
    )
    return {
        "feature_id_original": feature_id,
        "feature_id_standardized": "NCBI:GeneID:{}".format(target["gene_id"]),
        "ortholog_reference": "",
        "mapping_confidence": "exact" if exact_coordinates else "inferred",
        "mapping_release": release_id,
        "mapping_evidence": "exact_gene_coordinates" if exact_coordinates else "unique_same_strand_legacy_coverage_ge_0.8_or_reciprocal_ge_0.5",
        **common,
        "legacy_ncbi_gene_id": target["gene_id"],
        "legacy_ncbi_start": target["start"],
        "legacy_ncbi_end": target["end"],
        "overlap_bp": overlap,
        "legacy_coverage": round(source_coverage, 6),
        "legacy_ncbi_coverage": round(target_coverage, 6),
        "current_symbol": current["symbol"],
        "current_description": current["description"],
        "current_refseq_accession": current["seqid"],
        "current_start": current["start"],
        "current_end": current["end"],
        "current_strand": current["strand"],
        "candidate_overlap_count": 1,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--feature-table", required=True)
    parser.add_argument("--legacy-gff", required=True)
    parser.add_argument("--assembly-report", required=True)
    parser.add_argument("--legacy-ncbi-gff", required=True)
    parser.add_argument("--current-ncbi-gff", required=True)
    parser.add_argument("--gene-summary-json")
    parser.add_argument("--release-id", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    feature_ids = load_feature_ids(args.feature_table)
    aliases = load_assembly_aliases(args.assembly_report)
    legacy_genes = load_legacy_genes(args.legacy_gff, feature_ids)
    legacy_ncbi_by_seq, _ = load_ncbi_genes(args.legacy_ncbi_gff)
    _, current_by_gene_id = load_ncbi_genes(args.current_ncbi_gff)
    gene_replacements = load_gene_replacements(args.gene_summary_json)

    records = [
        map_feature(
            feature_id,
            legacy_genes.get(feature_id),
            aliases,
            legacy_ncbi_by_seq,
            current_by_gene_id,
            gene_replacements,
            args.release_id,
        )
        for feature_id in sorted(feature_ids)
    ]

    # Do not allow two selected legacy features to collapse onto one current
    # Gene ID. Such cases can represent annotation merges, but treating them as
    # resolved would double-count a standardized feature downstream.
    target_to_records = defaultdict(list)
    for record in records:
        if record["mapping_confidence"] != "unresolved":
            target_to_records[record["feature_id_standardized"]].append(record)
    for target_id, target_records in target_to_records.items():
        if len(target_records) < 2:
            continue
        for record in target_records:
            record["candidate_current_gene_id"] = target_id
            record["feature_id_standardized"] = record["feature_id_original"]
            record["mapping_confidence"] = "unresolved"
            record["mapping_evidence"] = "multiple_selected_legacy_features_map_to_same_current_gene"

    fieldnames = [
        "feature_id_original",
        "feature_id_standardized",
        "ortholog_reference",
        "mapping_confidence",
        "mapping_release",
        "mapping_evidence",
        "candidate_current_gene_id",
        "legacy_scaffold",
        "legacy_refseq_accession",
        "legacy_start",
        "legacy_end",
        "legacy_strand",
        "legacy_description",
        "legacy_ncbi_gene_id",
        "ncbi_replacement_gene_id",
        "legacy_ncbi_start",
        "legacy_ncbi_end",
        "overlap_bp",
        "legacy_coverage",
        "legacy_ncbi_coverage",
        "current_symbol",
        "current_description",
        "current_refseq_accession",
        "current_start",
        "current_end",
        "current_strand",
        "candidate_overlap_count",
    ]
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)

    counts = defaultdict(int)
    for record in records:
        counts[record["mapping_confidence"]] += 1
    print("wrote {} mappings to {}".format(len(records), output_path))
    print("confidence counts: {}".format(dict(sorted(counts.items()))))


if __name__ == "__main__":
    main()
