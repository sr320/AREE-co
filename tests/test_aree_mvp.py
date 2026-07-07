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


def test_identifier_mapping_confidence_assignment():
    exact = map_identifier("CGI_10001")
    unresolved = map_identifier("NOT_IN_MAP")
    assert exact["mapping_confidence"] == "exact"
    assert unresolved["mapping_confidence"] == "unresolved"


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

