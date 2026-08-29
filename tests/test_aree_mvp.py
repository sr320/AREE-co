from pathlib import Path

import pytest
import yaml

from aree.harmonize.identifiers import map_identifier
from aree.harmonize.processed import harmonize_demo, harmonize_processed
from aree.intake.registry import register_study
from aree.meta_analysis.random_effects import random_effects, run_meta_analysis
from aree.prioritize.scoring import score_candidates
from aree.reporting.demo_report import build_demo_report
from aree.reporting.evidence_cards import build_evidence_cards
from aree.validation.schemas import validate_study_file


ROOT = Path(__file__).resolve().parents[1]


def test_schema_validation_accepts_demo_study():
    study = validate_study_file(ROOT / "registry/studies/CGIG_HEAT_RNASEQ_001.yaml")
    assert study["study_id"] == "CGIG_HEAT_RNASEQ_001"


def test_schema_validation_rejects_malformed_metadata(tmp_path):
    bad = {"study_id": "bad id"}
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(bad))
    with pytest.raises(Exception):
        validate_study_file(path)


def test_duplicate_study_ids_are_rejected(tmp_path):
    registry = tmp_path / "registry.csv"
    study = ROOT / "registry/studies/CGIG_HEAT_RNASEQ_001.yaml"
    register_study(study, registry)
    with pytest.raises(ValueError):
        register_study(study, registry)


def test_required_provenance_fields_present():
    study = validate_study_file(ROOT / "registry/studies/CGIG_HEAT_RNASEQ_001.yaml")
    assert study["provenance"]["source_links"]
    assert study["provenance"]["curator"]
    assert study["provenance"]["curation_date"]


def test_second_real_study_manifest_matches_published_run_totals():
    import pandas as pd

    study = validate_study_file(
        ROOT / "registry/studies/CGIG_THERMOTOL_RNASEQ_PRJNA694496.yaml"
    )
    manifest = pd.read_csv(
        ROOT / "data/manifests/CGIG_THERMOTOL_RNASEQ_PRJNA694496_runs.tsv",
        sep="\t",
    )
    assert study["analysis_status"] == "raw_reanalysis_manifest_ready"
    assert len(manifest) == 6
    assert manifest["run_accession"].is_unique
    assert manifest["condition"].value_counts().to_dict() == {"selected": 3, "control": 3}
    assert manifest["read_count"].sum() == 132413989
    assert manifest[["fastq_1_bytes", "fastq_2_bytes"]].to_numpy().sum() == 16252181463
    assert manifest["assignment_evidence"].str.contains("exactly match").all()


def test_fastq_preflight_uses_nearest_existing_ancestor(tmp_path):
    from scripts.prepare_prjna694496_fastqs import nearest_existing_ancestor

    nested_destination = tmp_path / "not-yet-created" / "fastq"
    assert nearest_existing_ancestor(nested_destination) == tmp_path.resolve()


def test_prjna694496_workflow_sheets_and_tx2gene_are_schema_safe(tmp_path):
    import csv
    import gzip

    from scripts.prepare_gcf963853765_reference import build_tx2gene
    from scripts.prepare_prjna694496_fastqs import (
        load_manifest,
        write_design_samplesheet,
        write_nfcore_samplesheet,
    )

    rows = load_manifest(
        ROOT / "data/manifests/CGIG_THERMOTOL_RNASEQ_PRJNA694496_runs.tsv"
    )
    nfcore = tmp_path / "nfcore.csv"
    design = tmp_path / "design.csv"
    write_nfcore_samplesheet(rows, tmp_path / "fastq", nfcore)
    write_design_samplesheet(rows, design)
    with nfcore.open(newline="") as handle:
        nfcore_rows = list(csv.DictReader(handle))
    with design.open(newline="") as handle:
        design_rows = list(csv.DictReader(handle))
    assert list(nfcore_rows[0]) == ["sample", "fastq_1", "fastq_2", "strandedness"]
    assert list(design_rows[0]) == ["sample", "condition", "replicate", "run_accession"]

    gff = tmp_path / "mini.gff.gz"
    with gzip.open(gff, "wt") as handle:
        handle.write(
            "NC_1\tGnomon\tmRNA\t1\t10\t.\t+\t.\t"
            "ID=rna-XM_1.1;Dbxref=GeneID:123,GenBank:XM_1.1;"
            "gene=example;transcript_id=XM_1.1\n"
        )
    tx2gene = tmp_path / "tx2gene.tsv"
    build_tx2gene(gff, tx2gene)
    assert tx2gene.read_text().splitlines() == [
        "transcript_id\tgene_id\tgene_symbol",
        "XM_1.1\tNCBI:GeneID:123\texample",
    ]


def test_identifier_mapping_confidence_assignment():
    exact = map_identifier("CGI_10001")
    unresolved = map_identifier("NOT_IN_MAP")
    assert exact["mapping_confidence"] == "exact"
    assert unresolved["mapping_confidence"] == "unresolved"


def test_explicit_mapping_release_is_preserved(tmp_path):
    import pandas as pd

    mapping_path = tmp_path / "mapping.tsv"
    pd.DataFrame(
        [
            {
                "feature_id_original": "CGI_10001",
                "feature_id_standardized": "NCBI:GeneID:123",
                "ortholog_reference": "",
                "mapping_confidence": "inferred",
                "mapping_release": "test_release_v1",
                "mapping_evidence": "test_evidence",
            }
        ]
    ).to_csv(mapping_path, sep="\t", index=False)
    output = tmp_path / "evidence.tsv"
    harmonize_processed(
        "CGIG_HEAT_RNASEQ_001",
        ROOT / "data/demo/processed/CGIG_HEAT_RNASEQ_001_rnaseq.tsv",
        output,
        mapping_path,
    )
    evidence = pd.read_csv(output, sep="\t")
    mapped = evidence[evidence["feature_id_original"] == "CGI_10001"].iloc[0]
    assert mapped["feature_id_standardized"] == "NCBI:GeneID:123"
    assert mapped["mapping_release"] == "test_release_v1"
    assert mapped["mapping_evidence"] == "test_evidence"


def test_real_study_mapping_release_is_unique_and_conservative():
    import pandas as pd

    mapping = pd.read_csv(
        ROOT / "data/mappings/cgigas_cgi_to_ncbi_gene_rs2024_06_v1.tsv",
        sep="\t",
    )
    resolved = mapping[mapping["mapping_confidence"] == "inferred"]
    assert len(mapping) == 150
    assert mapping["feature_id_original"].is_unique
    assert len(resolved) == 60
    assert resolved["feature_id_standardized"].is_unique
    assert resolved["feature_id_standardized"].str.startswith("NCBI:GeneID:").all()


def test_harmonize_processed_outputs_provenance(tmp_path):
    output = tmp_path / "evidence.tsv"
    path = harmonize_processed(
        "CGIG_HEAT_RNASEQ_001",
        ROOT / "data/demo/processed/CGIG_HEAT_RNASEQ_001_rnaseq.tsv",
        output,
    )
    text = path.read_text()
    assert "input_checksum" in text
    assert "workflow_version" in text


def test_effect_size_meta_analysis_calculates_pooled_effect():
    import pandas as pd

    group = pd.DataFrame(
        {
            "effect_size": [1.0, 1.5],
            "standard_error": [0.2, 0.3],
            "study_id": ["A", "B"],
        }
    )
    result = random_effects(group)
    assert result["n_studies"] == 2
    assert result["pooled_effect"] > 1.0


def test_demo_meta_analysis_and_scoring(tmp_path):
    evidence = harmonize_demo()
    meta = run_meta_analysis(evidence_path=evidence, output_path=tmp_path / "meta.tsv")
    scores = score_candidates(evidence_path=evidence, meta_path=meta, output_path=tmp_path / "scores.tsv")
    assert meta.exists()
    assert scores.exists()
    assert "candidate_id" in scores.read_text()


def test_candidate_score_reproducibility(tmp_path):
    evidence = harmonize_demo()
    meta = run_meta_analysis(evidence_path=evidence, output_path=tmp_path / "meta.tsv")
    one = score_candidates(evidence_path=evidence, meta_path=meta, output_path=tmp_path / "scores1.tsv").read_text()
    two = score_candidates(evidence_path=evidence, meta_path=meta, output_path=tmp_path / "scores2.tsv").read_text()
    assert one == two


def test_evidence_card_generation(tmp_path):
    evidence = harmonize_demo()
    meta = run_meta_analysis(evidence_path=evidence, output_path=tmp_path / "meta.tsv")
    scores = score_candidates(evidence_path=evidence, meta_path=meta, output_path=tmp_path / "scores.tsv")
    cards = build_evidence_cards(evidence_path=evidence, scores_path=scores, output_dir=tmp_path / "cards")
    assert cards
    assert "not a validated biomarker" in cards[0].read_text()


def test_demo_report_build(tmp_path):
    evidence = harmonize_demo()
    meta = run_meta_analysis(evidence_path=evidence, output_path=tmp_path / "meta.tsv")
    score_candidates(evidence_path=evidence, meta_path=meta, output_path=ROOT / "data/demo/candidate_scores.tsv")
    report = build_demo_report(output_path=tmp_path / "report.md")
    assert report.exists()
    assert "AREE Demo Report" in report.read_text()
