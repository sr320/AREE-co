"""Multi-species bivalve OrthoFinder 3.1.5 ortholog mapping builders."""

import csv
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple


DEFAULT_RELEASE = "orthofinder3.1.5_8bivalves_2026-11-26_v1"
TARGET_SPECIES_CODE = "mgigas"
TARGET_ASSEMBLY = "GCF_963853765.1"

SPECIES_METADATA = {
    "mgigas": ("Magallana gigas", "GCF_963853765.1"),
    "cvirg": ("Crassostrea virginica", "GCF_053477285.1"),
    "medul": ("Mytilus edulis", "GCF_963676685.1"),
    "mgall": ("Mytilus galloprovincialis", "GCF_965363235.1"),
    "mcali": ("Mytilus californianus", "GCF_021869535.1"),
    "mtros": ("Mytilus trossulus", "GCF_036588685.1"),
    "mmerc": ("Mercenaria mercenaria", "GCF_021730395.1"),
    "rphil": ("Ruditapes philippinarum", "GCF_026571515.1"),
}

MAPPING_COLUMNS = [
    "feature_id_original",
    "feature_id_standardized",
    "ortholog_reference",
    "mapping_confidence",
    "mapping_release",
    "mapping_evidence",
    "source_species",
    "source_assembly",
    "target_assembly",
    "orthogroup_single_copy_8sp",
    "n_source_genes_in_orthologue",
    "n_mgigas_genes_in_orthologue",
    "candidate_mgigas_gene_ids",
]


def _split_gene_list(cell: str) -> List[str]:
    if not cell:
        return []
    return [token.strip() for token in str(cell).split(",") if token.strip()]


def _format_gene_id(raw_id: str) -> str:
    raw_id = str(raw_id).strip()
    if raw_id.startswith("NCBI:GeneID:"):
        return raw_id
    if raw_id.isdigit():
        return "NCBI:GeneID:{}".format(raw_id)
    return raw_id


def _sort_gene_key(raw_id: str) -> Tuple[int, str]:
    clean = raw_id.replace("NCBI:GeneID:", "")
    return (0, "{:012d}".format(int(clean))) if clean.isdigit() else (1, clean)


def load_single_copy_orthogroups(path: Path) -> Set[str]:
    """Load Orthogroups_SingleCopyOrthologues.txt into a set of OG IDs."""
    return {line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()}


def build_species_to_mgigas_map(
    orthologues_tsv: Path,
    species_code: str,
    single_copy_ogs: Optional[Set[str]] = None,
    output_path: Optional[Path] = None,
    release: str = DEFAULT_RELEASE,
) -> Dict[str, object]:
    """Parse ``mgigas__v__<sp>.tsv`` and emit a 1-row-per-source-gene AREE mapping TSV."""
    if species_code not in SPECIES_METADATA or species_code == TARGET_SPECIES_CODE:
        raise ValueError("Unsupported partner species code: {!r}".format(species_code))
    species_name, source_assembly = SPECIES_METADATA[species_code]
    single_copy_ogs = single_copy_ogs or set()

    # Aggregate pairwise rows per source gene so every feature_id_original is strictly unique.
    per_source = {}  # type: Dict[str, Dict[str, object]]
    with Path(orthologues_tsv).open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Orthogroup", TARGET_SPECIES_CODE, species_code}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise ValueError("Missing columns {} in {}".format(sorted(required), orthologues_tsv))
        for row in reader:
            og = row["Orthogroup"].strip()
            mg_genes = _split_gene_list(row[TARGET_SPECIES_CODE])
            sp_genes = _split_gene_list(row[species_code])
            if not mg_genes or not sp_genes:
                continue
            for sp_gene in sp_genes:
                entry = per_source.setdefault(
                    sp_gene,
                    {"ogs": set(), "mg_targets": set(), "sp_siblings": set()},
                )
                entry["ogs"].add(og)
                entry["mg_targets"].update(mg_genes)
                entry["sp_siblings"].update(sp_genes)

    rows = []
    conf_counts = Counter()
    for sp_gene in sorted(per_source, key=_sort_gene_key):
        entry = per_source[sp_gene]
        mg_sorted = sorted(entry["mg_targets"], key=_sort_gene_key)
        sp_sorted = sorted(entry["sp_siblings"], key=_sort_gene_key)
        og_sorted = sorted(entry["ogs"])
        primary_og = og_sorted[0]
        is_sc = any(og in single_copy_ogs for og in og_sorted)

        n_src = len(sp_sorted)
        n_mg = len(mg_sorted)
        if n_src == 1 and n_mg == 1:
            conf = "one-to-one ortholog"
            ev = "orthofinder3.1.5_1to1_single_copy_8sp" if is_sc else "orthofinder3.1.5_1to1"
        elif n_src > 1 and n_mg == 1:
            conf = "many-to-one ortholog"
            ev = "orthofinder3.1.5_many_to_one"
        elif n_src == 1 and n_mg > 1:
            conf = "one-to-many ortholog"
            ev = "orthofinder3.1.5_one_to_many"
        else:
            conf = "inferred"
            ev = "orthofinder3.1.5_many_to_many_lowest_target"

        conf_counts[conf] += 1
        rows.append({
            "feature_id_original": _format_gene_id(sp_gene),
            "feature_id_standardized": _format_gene_id(mg_sorted[0]),
            "ortholog_reference": primary_og,
            "mapping_confidence": conf,
            "mapping_release": release,
            "mapping_evidence": ev,
            "source_species": species_name,
            "source_assembly": source_assembly,
            "target_assembly": TARGET_ASSEMBLY,
            "orthogroup_single_copy_8sp": "true" if is_sc else "false",
            "n_source_genes_in_orthologue": n_src,
            "n_mgigas_genes_in_orthologue": n_mg,
            "candidate_mgigas_gene_ids": ";".join(_format_gene_id(g) for g in mg_sorted),
        })

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline="", encoding="utf-8") as out:
            writer = csv.DictWriter(out, fieldnames=MAPPING_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    return {
        "species_code": species_code,
        "species": species_name,
        "source_assembly": source_assembly,
        "total_mapped_genes": len(rows),
        "one_to_one": conf_counts["one-to-one ortholog"],
        "many_to_one": conf_counts["many-to-one ortholog"],
        "one_to_many": conf_counts["one-to-many ortholog"],
        "many_to_many_inferred": conf_counts["inferred"],
        "output_path": str(output_path) if output_path else None,
    }


def build_mgigas_orthogroup_index(
    orthogroups_tsv: Path,
    single_copy_ogs: Optional[Set[str]] = None,
    output_path: Optional[Path] = None,
) -> int:
    """Build a *Magallana gigas*-anchored orthogroup index TSV across all 8 cataloged bivalves."""
    single_copy_ogs = single_copy_ogs or set()
    partner_codes = ["cvirg", "medul", "mgall", "mcali", "mtros", "mmerc", "rphil"]
    all_codes = [TARGET_SPECIES_CODE] + partner_codes
    rows = []

    with Path(orthogroups_tsv).open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            og = row["Orthogroup"].strip()
            mg_genes = _split_gene_list(row.get(TARGET_SPECIES_CODE, ""))
            if not mg_genes:
                continue
            species_present = sum(1 for c in all_codes if _split_gene_list(row.get(c, "")))
            partner_cells = {
                "{}_genes".format(c): ";".join(_format_gene_id(g) for g in _split_gene_list(row.get(c, "")))
                for c in partner_codes
            }
            is_sc = "true" if og in single_copy_ogs else "false"
            for mg_gene in sorted(mg_genes, key=_sort_gene_key):
                rec = {
                    "feature_id_standardized": _format_gene_id(mg_gene),
                    "ortholog_reference": og,
                    "is_single_copy_8sp": is_sc,
                    "species_present_count": species_present,
                    "mgigas_paralogs_in_og": len(mg_genes),
                }
                rec.update(partner_cells)
                rows.append(rec)

    rows.sort(key=lambda r: _sort_gene_key(r["feature_id_standardized"]))
    fieldnames = [
        "feature_id_standardized",
        "ortholog_reference",
        "is_single_copy_8sp",
        "species_present_count",
        "mgigas_paralogs_in_og",
    ] + ["{}_genes".format(c) for c in partner_codes]

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with output_path.open("w", newline="", encoding="utf-8") as out:
            writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    return len(rows)


def build_all_multispecies_maps(
    results_dir: Path,
    output_dir: Path,
    release: str = DEFAULT_RELEASE,
    species_codes: Optional[Iterable[str]] = None,
) -> List[Dict[str, object]]:
    """Build pairwise ortholog mapping TSVs, the *M. gigas* orthogroup index, and ``release_summary.tsv``."""
    results_dir = Path(results_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    sc_path = results_dir / "Orthogroups" / "Orthogroups_SingleCopyOrthologues.txt"
    single_copy = load_single_copy_orthogroups(sc_path) if sc_path.exists() else set()

    og_path = results_dir / "Orthogroups" / "Orthogroups.tsv"
    n_index = build_mgigas_orthogroup_index(
        og_path,
        single_copy_ogs=single_copy,
        output_path=output_dir / "mgigas_orthogroups_v1.tsv",
    )

    codes = list(species_codes) if species_codes is not None else [c for c in SPECIES_METADATA if c != TARGET_SPECIES_CODE]
    summaries = []
    for code in codes:
        pair_path = results_dir / "Orthologues" / "Orthologues_{}".format(TARGET_SPECIES_CODE) / "{}__{}__{}.tsv".format(
            TARGET_SPECIES_CODE, "v", code
        )
        out_tsv = output_dir / "{}_to_{}_orthofinder_v1.tsv".format(code, TARGET_SPECIES_CODE)
        summary = build_species_to_mgigas_map(
            pair_path,
            species_code=code,
            single_copy_ogs=single_copy,
            output_path=out_tsv,
            release=release,
        )
        summary["single_copy_8sp_orthogroups"] = len(single_copy)
        summary["mgigas_indexed_genes"] = n_index
        summaries.append(summary)

    summary_cols = [
        "species_code",
        "species",
        "source_assembly",
        "total_mapped_genes",
        "one_to_one",
        "many_to_one",
        "one_to_many",
        "many_to_many_inferred",
        "single_copy_8sp_orthogroups",
        "mgigas_indexed_genes",
    ]
    with (output_dir / "release_summary.tsv").open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=summary_cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summaries)

    return summaries
