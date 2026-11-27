"""Export AREE registry, harmonized evidence, synthesis, and catalogs to SQLite, Parquet, and JSON-LD."""

import json
import os
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import yaml

from aree.paths import root_path
from aree.raw.reference import list_references


REPO_URL = "https://github.com/sr320/AREE-co"

SCHEMA_SQL = """-- AREE-co Relational Export Schema (SQLite 3)

CREATE TABLE IF NOT EXISTS studies (
    study_id TEXT PRIMARY KEY,
    bioproject TEXT,
    species TEXT,
    taxonomy_id TEXT,
    stressor_class TEXT,
    stressor TEXT,
    phenotype TEXT,
    tissue TEXT,
    life_stage TEXT,
    data_status TEXT,
    biological_replication TEXT,
    study_design TEXT,
    experimental_treatment TEXT,
    control_condition TEXT,
    phenotype_direction TEXT,
    limitations_json TEXT
);

CREATE TABLE IF NOT EXISTS refseq_references (
    accession TEXT PRIMARY KEY,
    species TEXT,
    species_code TEXT,
    ncbi_taxid INTEGER,
    assembly_id TEXT,
    annotation_release TEXT,
    tx2gene_genes INTEGER,
    tx2gene_transcripts INTEGER
);

CREATE TABLE IF NOT EXISTS harmonized_evidence (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    study_id TEXT,
    feature_id_original TEXT,
    feature_id_standardized TEXT,
    feature_type TEXT,
    gene_symbol TEXT,
    description TEXT,
    sample_comparison TEXT,
    molecular_direction TEXT,
    effect_size REAL,
    effect_size_type TEXT,
    standard_error REAL,
    p_value REAL,
    adjusted_p_value REAL,
    sample_size INTEGER,
    tissue TEXT,
    life_stage TEXT,
    stressor TEXT,
    phenotype TEXT,
    resilience_classification TEXT,
    inference_status TEXT,
    ortholog_reference TEXT,
    mapping_confidence TEXT,
    quality_flags TEXT,
    workflow_version TEXT,
    genome_assembly TEXT,
    annotation_version TEXT,
    source_file TEXT,
    input_checksum TEXT
);

CREATE INDEX IF NOT EXISTS idx_ev_study ON harmonized_evidence(study_id);
CREATE INDEX IF NOT EXISTS idx_ev_feature ON harmonized_evidence(feature_id_standardized);
CREATE INDEX IF NOT EXISTS idx_ev_fdr ON harmonized_evidence(adjusted_p_value);

CREATE TABLE IF NOT EXISTS synthesis_strata (
    stratum_id TEXT PRIMARY KEY,
    title TEXT,
    phenotype TEXT,
    stressor TEXT,
    n_contrasts INTEGER,
    n_bioprojects INTEGER,
    n_independent_groups INTEGER,
    multi_group_pooled_features INTEGER,
    pooled_fdr_lt_0_05 INTEGER,
    concordant_fdr_lt_0_05 INTEGER,
    median_i2_percent REAL
);

CREATE TABLE IF NOT EXISTS synthesis_meta (
    stratum_id TEXT,
    feature_id_standardized TEXT,
    gene_symbol TEXT,
    description TEXT,
    n_studies INTEGER,
    n_independent_groups INTEGER,
    method TEXT,
    shared_control_rho REAL,
    pooled_effect REAL,
    pooled_standard_error REAL,
    pooled_ci_low REAL,
    pooled_ci_high REAL,
    p_value REAL,
    adjusted_p_value REAL,
    i2_percent REAL,
    tau2 REAL,
    direction_concordant TEXT,
    PRIMARY KEY (stratum_id, feature_id_standardized)
);

CREATE INDEX IF NOT EXISTS idx_synth_fdr ON synthesis_meta(stratum_id, adjusted_p_value);

CREATE TABLE IF NOT EXISTS ortholog_release_summary (
    species_code TEXT PRIMARY KEY,
    species TEXT,
    source_assembly TEXT,
    total_mapped_genes INTEGER,
    one_to_one INTEGER,
    many_to_one INTEGER,
    one_to_many INTEGER,
    many_to_many_inferred INTEGER,
    single_copy_8sp_orthogroups INTEGER,
    mgigas_indexed_genes INTEGER
);
"""


def _resolve_root(root: Optional[Path] = None) -> Path:
    return Path(root) if root else root_path()


def _load_studies_df(root: Path, include_simulated: bool = False) -> pd.DataFrame:
    registry_path = root / "registry" / "study_registry.csv"
    if not registry_path.exists():
        return pd.DataFrame()
    reg = pd.read_csv(registry_path, dtype=str)
    rows = []
    for row in reg.itertuples(index=False):
        status = str(getattr(row, "data_status", "") or "")
        if not include_simulated and status.lower().startswith("simulated"):
            continue
        sid = row.study_id
        yml_path = root / "registry" / "studies" / "{}.yaml".format(sid)
        meta = {}
        if yml_path.exists():
            meta = yaml.safe_load(yml_path.read_text(encoding="utf-8")) or {}
        accessions = meta.get("accessions", {})
        rows.append({
            "study_id": sid,
            "bioproject": accessions.get("bioproject", ""),
            "species": getattr(row, "species", meta.get("species", "")),
            "taxonomy_id": str(meta.get("taxonomy_id", "")),
            "stressor_class": meta.get("stressor_class", getattr(row, "stressor_class", "")),
            "stressor": getattr(row, "stressor", meta.get("stressor", "")),
            "phenotype": getattr(row, "phenotype", meta.get("phenotype", "")),
            "tissue": getattr(row, "tissue", meta.get("tissue", "")),
            "life_stage": meta.get("life_stage", ""),
            "data_status": status,
            "biological_replication": meta.get("biological_replication", ""),
            "study_design": meta.get("study_design", ""),
            "experimental_treatment": meta.get("experimental_treatment", ""),
            "control_condition": meta.get("control_condition", ""),
            "phenotype_direction": meta.get("phenotype_direction", ""),
            "limitations_json": json.dumps(meta.get("limitations", [])),
        })
    return pd.DataFrame(rows)


def _load_references_df(root: Path) -> pd.DataFrame:
    cat_path = root / "config" / "references.yaml"
    if not cat_path.exists():
        return pd.DataFrame()
    rows = []
    for spec in list_references(cat_path):
        rows.append({
            "accession": spec["accession"],
            "species": spec["species"],
            "species_code": spec.get("species_code", ""),
            "ncbi_taxid": int(spec["ncbi_taxid"]) if spec.get("ncbi_taxid") is not None else None,
            "assembly_id": spec.get("assembly_id", ""),
            "annotation_release": spec.get("annotation_release", ""),
            "tx2gene_genes": spec.get("tx2gene_genes"),
            "tx2gene_transcripts": spec.get("tx2gene_transcripts"),
        })
    return pd.DataFrame(rows)


def _load_harmonized_evidence_df(root: Path, study_ids: Optional[List[str]] = None) -> pd.DataFrame:
    harm_dir = root / "data" / "harmonized"
    if not harm_dir.exists():
        return pd.DataFrame()
    frames = []
    allowed = set(study_ids) if study_ids is not None else None
    for path in sorted(harm_dir.glob("*.tsv")):
        sid = path.stem
        if sid.startswith("demo_") or (allowed is not None and sid not in allowed):
            continue
        frames.append(pd.read_csv(path, sep="\t", low_memory=False))
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    for col in ("effect_size", "standard_error", "p_value", "adjusted_p_value"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if "sample_size" in df.columns:
        df["sample_size"] = pd.to_numeric(df["sample_size"], errors="coerce")
    return df


def _load_synthesis_dfs(root: Path):
    strata_path = root / "reports" / "synthesis" / "stratum_summary.tsv"
    strata_df = pd.read_csv(strata_path, sep="\t") if strata_path.exists() else pd.DataFrame()
    synth_dir = root / "data" / "synthesis"
    frames = []
    if synth_dir.exists():
        for path in sorted(synth_dir.glob("*_meta.tsv")):
            sid = path.stem[:-5] if path.stem.endswith("_meta") else path.stem
            df = pd.read_csv(path, sep="\t", low_memory=False)
            if "stratum_id" not in df.columns:
                df.insert(0, "stratum_id", sid)
            frames.append(df)
    meta_df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return strata_df, meta_df


def _load_ortholog_summary_df(root: Path) -> pd.DataFrame:
    path = root / "data" / "mappings" / "orthologs" / "release_summary.tsv"
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame()


def _collect_tables(root: Path, include_simulated: bool = False) -> Dict[str, pd.DataFrame]:
    studies_df = _load_studies_df(root, include_simulated=include_simulated)
    study_ids = studies_df["study_id"].tolist() if not studies_df.empty else []
    ev_df = _load_harmonized_evidence_df(root, study_ids=study_ids)
    strata_df, meta_df = _load_synthesis_dfs(root)
    return {
        "studies": studies_df,
        "refseq_references": _load_references_df(root),
        "harmonized_evidence": ev_df,
        "synthesis_strata": strata_df,
        "synthesis_meta": meta_df,
        "ortholog_release_summary": _load_ortholog_summary_df(root),
    }


def export_schema_sql(output_path: Path) -> Path:
    """Write standalone SQLite DDL schema to ``output_path``."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(SCHEMA_SQL, encoding="utf-8")
    return output_path


def export_sqlite(
    output_path: Path,
    root: Optional[Path] = None,
    include_simulated: bool = False,
) -> dict:
    """Build a normalized SQLite database covering studies, references, evidence, synthesis, and orthologs."""
    base = _resolve_root(root)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()

    tables = _collect_tables(base, include_simulated=include_simulated)
    conn = sqlite3.connect(str(tmp_path))
    counts = {}
    try:
        conn.executescript(SCHEMA_SQL)
        for name, df in tables.items():
            if df.empty:
                counts[name] = 0
                continue
            cur_cols = [info[1] for info in conn.execute("PRAGMA table_info({})".format(name)).fetchall()]
            usable_cols = [c for c in df.columns if c in cur_cols]
            sub = df[usable_cols].copy()
            if name == "synthesis_meta" and "direction_concordant" in sub.columns:
                sub["direction_concordant"] = sub["direction_concordant"].astype(str)
            sub.to_sql(name, conn, if_exists="append", index=False)
            counts[name] = int(len(sub))
        conn.commit()
    finally:
        conn.close()

    os.replace(str(tmp_path), str(output_path))
    return {
        "format": "sqlite",
        "output_path": str(output_path),
        "bytes": output_path.stat().st_size,
        "tables": counts,
    }


def _frictionless_type(dtype) -> str:
    s = str(dtype)
    if "int" in s:
        return "integer"
    if "float" in s:
        return "number"
    if "bool" in s:
        return "boolean"
    return "string"


def export_parquet_bundle(
    output_dir: Path,
    root: Optional[Path] = None,
    include_simulated: bool = False,
) -> dict:
    """Export all tables as ZSTD Parquet files with a Frictionless ``datapackage.json`` manifest."""
    try:
        import pyarrow  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "pyarrow is required for Parquet export; install with `pip install -e '.[site]'`"
        ) from exc

    base = _resolve_root(root)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tables = _collect_tables(base, include_simulated=include_simulated)

    resources = []
    counts = {}
    for name, df in tables.items():
        if df.empty:
            counts[name] = 0
            continue
        pq_name = "{}.parquet".format(name)
        pq_path = output_dir / pq_name
        df.to_parquet(pq_path, index=False, compression="zstd")
        counts[name] = int(len(df))
        resources.append({
            "name": name,
            "path": pq_name,
            "format": "parquet",
            "mediatype": "application/vnd.apache.parquet",
            "row_count": int(len(df)),
            "schema": {
                "fields": [
                    {"name": str(col), "type": _frictionless_type(df[col].dtype)}
                    for col in df.columns
                ]
            },
        })

    datapackage = {
        "profile": "tabular-data-package",
        "name": "aree-co-evidence-bundle",
        "title": "Aquaculture Resilience Evidence Companion (AREC) Database Bundle",
        "license": "CC-BY-4.0",
        "homepage": REPO_URL,
        "resources": resources,
    }
    pkg_path = output_dir / "datapackage.json"
    pkg_path.write_text(json.dumps(datapackage, indent=2) + "\n", encoding="utf-8")
    return {
        "format": "parquet",
        "output_dir": str(output_dir),
        "datapackage_path": str(pkg_path),
        "tables": counts,
    }


def export_jsonld_catalog(
    output_path: Path,
    root: Optional[Path] = None,
) -> dict:
    """Write a Schema.org DataCatalog JSON-LD manifest of real studies, strata, references, and orthologs."""
    base = _resolve_root(root)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    studies_df = _load_studies_df(base, include_simulated=False)
    refs_df = _load_references_df(base)
    strata_df, _ = _load_synthesis_dfs(base)
    ortho_df = _load_ortholog_summary_df(base)

    datasets = []
    if not studies_df.empty:
        for row in studies_df.itertuples(index=False):
            harm_path = base / "data" / "harmonized" / "{}.tsv".format(row.study_id)
            if not harm_path.exists():
                continue
            datasets.append({
                "@type": "Dataset",
                "@id": "{}/blob/main/data/harmonized/{}.tsv".format(REPO_URL, row.study_id),
                "identifier": row.study_id,
                "name": "{} — {} ({})".format(row.study_id, row.species, row.stressor),
                "description": "Harmonized differential molecular association table for BioProject {} ({}; {}).".format(
                    row.bioproject, row.stressor, row.phenotype
                ),
                "distribution": {
                    "@type": "DataDownload",
                    "encodingFormat": "text/tab-separated-values",
                    "contentUrl": "{}/raw/main/data/harmonized/{}.tsv".format(REPO_URL, row.study_id),
                },
                "variableMeasured": [row.stressor, row.phenotype, row.tissue],
            })

    if not strata_df.empty:
        for row in strata_df.itertuples(index=False):
            datasets.append({
                "@type": "Dataset",
                "@id": "{}/blob/main/data/synthesis/{}_meta.tsv".format(REPO_URL, row.stratum_id),
                "identifier": row.stratum_id,
                "name": "Synthesis stratum: {}".format(row.title),
                "description": "REML + shared-control block-GLS cross-study synthesis across {} contrasts and {} BioProjects ({} concordant FDR < 0.05).".format(
                    int(row.n_contrasts), int(row.n_bioprojects), int(row.concordant_fdr_lt_0_05)
                ),
                "distribution": {
                    "@type": "DataDownload",
                    "encodingFormat": "text/tab-separated-values",
                    "contentUrl": "{}/raw/main/data/synthesis/{}_meta.tsv".format(REPO_URL, row.stratum_id),
                },
            })

    if not ortho_df.empty:
        datasets.append({
            "@type": "Dataset",
            "@id": "{}/tree/main/data/mappings/orthologs".format(REPO_URL),
            "identifier": "orthofinder3.1.5_8bivalves_2026-11-26_v1",
            "name": "8-Bivalve OrthoFinder 3.1.5 Ortholog Mapping Release",
            "description": "Pairwise bivalve-to-Magallana gigas gene ortholog maps across {} species.".format(len(ortho_df)),
            "distribution": {
                "@type": "DataDownload",
                "encodingFormat": "text/tab-separated-values",
                "contentUrl": "{}/raw/main/data/mappings/orthologs/release_summary.tsv".format(REPO_URL),
            },
        })

    catalog = {
        "@context": "https://schema.org",
        "@type": "DataCatalog",
        "@id": REPO_URL,
        "name": "Aquaculture Resilience Evidence Companion (AREC)",
        "description": "Curated, provenance-tracked molecular aquaculture resilience evidence catalog covering {} real study entries, {} synthesis strata, and {} RefSeq bivalve assemblies.".format(
            len(studies_df), len(strata_df), len(refs_df)
        ),
        "url": "https://sr320.github.io/AREE-co/",
        "license": "https://creativecommons.org/licenses/by/4.0/",
        "hasPart": datasets,
    }
    output_path.write_text(json.dumps(catalog, indent=2) + "\n", encoding="utf-8")
    return {
        "format": "jsonld",
        "output_path": str(output_path),
        "n_datasets": len(datasets),
    }
