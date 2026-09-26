from pathlib import Path
import sys

import jsonschema
import pandas as pd
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
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def test_schema_validation_accepts_demo_study():
    study = validate_study_file(ROOT / "registry/studies/CGIG_HEAT_RNASEQ_001.yaml")
    assert study["study_id"] == "CGIG_HEAT_RNASEQ_001"


def test_schema_validation_rejects_malformed_metadata(tmp_path):
    bad = {"study_id": "bad id"}
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(bad))
    with pytest.raises(jsonschema.ValidationError):
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
    assert study["analysis_status"] == "raw_reanalysis_harmonized"
    assert study["genome_assembly"] == "GCF_963853765.1"
    assert study["annotation_version"] == "RS_2024_06"
    assert study["data_availability"]["processed"] is True
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


def test_salmon_qc_summary_follows_design_order(tmp_path):
    import csv
    import json

    from scripts.summarize_salmon_qc import summarize_salmon_qc

    design = tmp_path / "design.csv"
    design.write_text(
        "sample,condition,replicate,run_accession\n"
        "selected_A,selected,A,SRR1\n"
    )
    metadata_dir = tmp_path / "salmon" / "selected_A" / "aux_info"
    metadata_dir.mkdir(parents=True)
    (metadata_dir / "meta_info.json").write_text(
        json.dumps(
            {
                "salmon_version": "1.10.3",
                "library_types": ["IU"],
                "num_processed": 100,
                "num_mapped": 75,
                "percent_mapped": 75.0,
                "num_decoy_fragments": 5,
                "frag_length_mean": 250.0,
                "frag_length_sd": 40.0,
                "seq_bias_correct": True,
                "gc_bias_correct": True,
            }
        )
    )
    output = summarize_salmon_qc(tmp_path / "salmon", design, tmp_path / "qc.tsv")
    with output.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["sample"] == "selected_A"
    assert rows[0]["library_type"] == "IU"
    assert rows[0]["mapping_rate_percent"] == "75.000000"
    assert rows[0]["decoy_rate_percent"] == "5.000000"


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


def test_meta_analysis_p_value_keeps_precision_for_large_effects():
    import pandas as pd

    group = pd.DataFrame({"effect_size": [10.0, 10.0], "standard_error": [0.5, 0.5], "study_id": ["A", "B"]})
    p_value = random_effects(group)["p_value"]
    assert 0.0 < p_value < 1e-150


def test_demo_meta_analysis_and_scoring(tmp_path):
    evidence = harmonize_demo(tmp_path / "evidence.tsv")
    meta = run_meta_analysis(evidence_path=evidence, output_path=tmp_path / "meta.tsv")
    scores = score_candidates(evidence_path=evidence, output_path=tmp_path / "scores.tsv")
    assert meta.exists()
    assert scores.exists()
    assert "candidate_id" in scores.read_text()


def test_candidate_score_reproducibility(tmp_path):
    evidence = harmonize_demo(tmp_path / "evidence.tsv")
    one = score_candidates(evidence_path=evidence, output_path=tmp_path / "scores1.tsv").read_text()
    two = score_candidates(evidence_path=evidence, output_path=tmp_path / "scores2.tsv").read_text()
    assert one == two


def test_evidence_card_generation(tmp_path):
    evidence = harmonize_demo(tmp_path / "evidence.tsv")
    scores = score_candidates(evidence_path=evidence, output_path=tmp_path / "scores.tsv")
    cards = build_evidence_cards(evidence_path=evidence, scores_path=scores, output_dir=tmp_path / "cards")
    assert cards
    assert "not a validated biomarker" in cards[0].read_text()


def test_demo_report_build(tmp_path):
    evidence = harmonize_demo(tmp_path / "evidence.tsv")
    scores = score_candidates(evidence_path=evidence, output_path=tmp_path / "scores.tsv")
    report = build_demo_report(output_path=tmp_path / "report.md", evidence_path=evidence, scores_path=scores)
    assert report.exists()
    assert "AREE Demo Report" in report.read_text()


def test_harmonize_rerun_is_byte_identical_and_portable(tmp_path):
    evidence = harmonize_demo(tmp_path / "evidence.tsv")
    table = pd.read_csv(evidence, sep="\t")
    evidence.write_text(evidence.read_text().replace(table["date_generated"].iloc[0], "2000-01-01"))
    before = evidence.read_text()
    harmonize_demo(evidence)
    assert evidence.read_text() == before
    assert all(not str(path).startswith("/") for path in table["source_file"])


def test_changed_input_refreshes_generation_date(tmp_path):
    source = ROOT / "data/demo/processed/CGIG_HEAT_RNASEQ_001_rnaseq.tsv"
    evidence = harmonize_processed("CGIG_HEAT_RNASEQ_001", source, tmp_path / "evidence.tsv")
    table = pd.read_csv(evidence, sep="\t")
    evidence.write_text(evidence.read_text().replace(table["date_generated"].iloc[0], "2000-01-01"))
    changed = tmp_path / "changed.tsv"
    processed = pd.read_csv(source, sep="\t")
    processed.loc[0, "effect_size"] += 1.0
    processed.to_csv(changed, sep="\t", index=False)
    harmonize_processed("CGIG_HEAT_RNASEQ_001", changed, evidence)
    assert "2000-01-01" not in set(pd.read_csv(evidence, sep="\t")["date_generated"])


def test_real_study_requires_explicit_output():
    with pytest.raises(ValueError, match="explicit output path"):
        harmonize_processed(
            "CGIG_HEAT_RNASEQ_PRJNA516762", ROOT / "data/processed/CGIG_HEAT_RNASEQ_PRJNA516762_rnaseq.tsv"
        )


def test_cli_downstream_commands_accept_real_evidence_paths(tmp_path):
    from typer.testing import CliRunner

    from aree.cli import app

    evidence = tmp_path / "evidence.tsv"
    source = ROOT / "data/processed/CGIG_HEAT_RNASEQ_PRJNA516762_rnaseq.tsv"
    runner = CliRunner()
    missing_output = runner.invoke(app, ["harmonize", "--study", "CGIG_HEAT_RNASEQ_PRJNA516762", "--input", str(source)])
    assert missing_output.exit_code != 0
    commands = [
        ["harmonize", "--study", "CGIG_HEAT_RNASEQ_PRJNA516762", "--input", str(source),
         "--mapping", str(ROOT / "data/mappings/cgigas_cgi_to_ncbi_gene_rs2024_06_v1.tsv"), "--output", str(evidence)],
        ["meta-analyze", "--evidence", str(evidence), "--output", str(tmp_path / "meta.tsv")],
        ["build-evidence-cards", "--evidence", str(evidence),
         "--scores", str(tmp_path / "scores.tsv"), "--output-dir", str(tmp_path / "cards")],
        ["build-demo-report", "--evidence", str(evidence), "--scores", str(tmp_path / "scores.tsv"),
         "--output", str(tmp_path / "report.md")],
    ]
    for command in commands:
        result = runner.invoke(app, command)
        assert result.exit_code == 0, result.output
    scores = pd.read_csv(tmp_path / "scores.tsv", sep="\t")
    assert len(scores) == pd.read_csv(evidence, sep="\t")["feature_id_standardized"].nunique()
    assert list((tmp_path / "cards").glob("*.md"))


def _evidence_rows(tmp_path, rows):
    base = {
        "feature_type": "gene",
        "effect_size_type": "log2_fold_change",
        "sample_size": 10,
        "resilience_classification": "resilience_associated",
        "mapping_confidence": "exact",
        "quality_flags": "none",
        "adjusted_p_value": 0.01,
        "tissue": "gill",
        "life_stage": "adult",
        "phenotype": "thermal_tolerance",
        "stressor": "temperature",
        "molecular_direction": "up",
        "species": "Crassostrea gigas",
        "ortholog_reference": None,
    }
    path = tmp_path / "evidence.tsv"
    pd.DataFrame([dict(base, **row) for row in rows]).to_csv(path, sep="\t", index=False)
    return path


def test_meta_analysis_never_pools_different_effect_scales(tmp_path):
    from aree.meta_analysis.random_effects import meta_analysis_table

    evidence = pd.read_csv(_evidence_rows(tmp_path, [
        {"feature_id_standardized": "G1", "study_id": "A", "effect_size": 2.0, "standard_error": 0.2},
        {"feature_id_standardized": "G1", "study_id": "B", "effect_size": 0.1, "standard_error": 0.02,
         "effect_size_type": "methylation_difference"},
    ]), sep="\t")
    table = meta_analysis_table(evidence)
    assert sorted(table["effect_size_type"]) == ["log2_fold_change", "methylation_difference"]
    assert set(table["n_effects"]) == {1}


def test_scoring_summarizes_effects_within_each_scale(tmp_path):
    path = _evidence_rows(tmp_path, [
        {"feature_id_standardized": "G1", "study_id": "A", "effect_size": 2.0, "standard_error": 0.2},
        {"feature_id_standardized": "G1", "study_id": "B", "effect_size": -1.0, "standard_error": 0.2},
        {"feature_id_standardized": "G1", "study_id": "C", "effect_size": 0.1, "standard_error": 0.02,
         "effect_size_type": "methylation_difference"},
    ])
    from aree.prioritize.scoring import _typed_effect_summary

    magnitude, direction = _typed_effect_summary(pd.read_csv(path, sep="\t"))
    assert magnitude == pytest.approx(((1.5 / 2.5) + (0.1 / 2.5)) / 2)
    # The lone methylation effect cannot show consistency, so only the split log2FC pair counts.
    assert direction == pytest.approx(0.5)


def test_one_effect_per_assay_is_not_high_priority(tmp_path):
    path = _evidence_rows(tmp_path, [
        {"feature_id_standardized": "G1", "study_id": "A", "effect_size": 1.0, "standard_error": 0.1},
        {"feature_id_standardized": "G1", "study_id": "B", "effect_size": -0.2, "standard_error": 0.05,
         "feature_type": "genomic_region", "effect_size_type": "methylation_difference"},
    ])
    scores = pd.read_csv(score_candidates(evidence_path=path, output_path=tmp_path / "s.tsv"), sep="\t")
    row = scores.iloc[0]
    assert pd.isna(row["direction_consistency"])
    assert row["consistency_flag"] == "not replicated within an effect-size type"
    assert row["category"] == "Multi-omics convergence candidate"


def test_heterogeneity_penalty_comes_from_scored_evidence(tmp_path):
    from aree.meta_analysis.random_effects import meta_analysis_table

    def score_with_standard_error(se, name):
        path = _evidence_rows(tmp_path, [
            {"feature_id_standardized": "G1", "study_id": "A", "effect_size": 1.0, "standard_error": se},
            {"feature_id_standardized": "G1", "study_id": "B", "effect_size": 3.0, "standard_error": se},
        ])
        i2 = meta_analysis_table(pd.read_csv(path, sep="\t"))["i2_percent"].iloc[0]
        score = pd.read_csv(score_candidates(evidence_path=path, output_path=tmp_path / name), sep="\t")["score"].iloc[0]
        return i2, score

    # Only the standard errors differ, and they feed no other score component,
    # so the score gap is exactly the I2 penalty computed from this evidence.
    i2_wide, score_wide = score_with_standard_error(10.0, "wide.tsv")
    i2_tight, score_tight = score_with_standard_error(0.1, "tight.tsv")
    assert i2_wide == 0.0 and i2_tight > 90.0
    assert score_wide - score_tight == pytest.approx(i2_tight / 100.0 * 0.15, abs=1e-4)


def test_scoring_phenotype_filter_limits_evidence(tmp_path):
    path = _evidence_rows(tmp_path, [
        {"feature_id_standardized": "G1", "study_id": "A", "effect_size": 1.0, "standard_error": 0.1},
        {"feature_id_standardized": "G2", "study_id": "B", "effect_size": 1.0, "standard_error": 0.1,
         "phenotype": "survival"},
    ])
    scores = pd.read_csv(score_candidates(evidence_path=path, output_path=tmp_path / "s.tsv", phenotype="survival"), sep="\t")
    assert list(scores["candidate_id"]) == ["G2"]


def test_filtered_cards_are_scored_from_the_evidence_they_show(tmp_path):
    from aree.prioritize.scoring import score_table

    evidence = harmonize_demo(tmp_path / "evidence.tsv")
    cards = build_evidence_cards(phenotype="survival", evidence_path=evidence, output_dir=tmp_path / "cards")
    table = pd.read_csv(evidence, sep="\t")
    expected = score_table(table[table["phenotype"] == "survival"]).set_index("candidate_id")["score"]
    assert len(cards) == len(expected)
    for card in cards:
        candidate = card.read_text().splitlines()[0].replace("# Evidence Card: ", "")
        assert "- Candidate score: {}".format(expected[candidate]) in card.read_text()


def test_effects_without_standard_errors_are_reported_not_dropped(tmp_path):
    from aree.meta_analysis.random_effects import meta_analysis_table

    evidence = pd.read_csv(_evidence_rows(tmp_path, [
        {"feature_id_standardized": "G1", "study_id": "A", "effect_size": 1.0, "standard_error": 0.2},
        {"feature_id_standardized": "G1", "study_id": "B", "effect_size": 2.0, "standard_error": None},
        {"feature_id_standardized": "G2", "study_id": "B", "effect_size": 2.0, "standard_error": None},
    ]), sep="\t")
    table = meta_analysis_table(evidence).set_index("feature_id_standardized")
    assert table.loc["G1", "pooling_status"] == "single_effect"
    assert table.loc["G1", "n_effects_excluded"] == 1
    assert table.loc["G1", "excluded_study_ids"] == "B"
    assert table.loc["G2", "pooling_status"] == "no_standard_errors"
    assert table.loc["G2", "n_effects"] == 0
    assert pd.isna(table.loc["G2", "pooled_effect"])


def test_unpoolable_evidence_is_scored_and_surfaced(tmp_path):
    from typer.testing import CliRunner

    from aree.cli import app

    path = _evidence_rows(tmp_path, [
        {"feature_id_standardized": "G2", "study_id": "B", "effect_size": 2.0, "standard_error": None},
    ])
    scores = pd.read_csv(score_candidates(evidence_path=path, output_path=tmp_path / "s.tsv"), sep="\t")
    assert scores["score"].notna().all()
    result = CliRunner().invoke(app, ["meta-analyze", "--evidence", str(path), "--output", str(tmp_path / "m.tsv")])
    assert result.exit_code == 0, result.output
    assert "1 effects lack standard errors" in result.output and "studies: B" in result.output
    report = build_demo_report(output_path=tmp_path / "r.md", evidence_path=path, scores_path=tmp_path / "s.tsv")
    assert "| B | 1 | 0 | 1 |" in report.read_text()
    cards = build_evidence_cards(evidence_path=path, output_dir=tmp_path / "cards")
    assert "not pooled in meta-analysis): B" in cards[0].read_text()


def test_card_filenames_are_portable():
    from aree.reporting.evidence_cards import card_filename

    assert card_filename("NCBI:LOC105317001") == "NCBI_LOC105317001.md"
    assert card_filename("a/b\\c") == "a_b_c.md"
    assert card_filename("CON") == "_CON.md"
    assert card_filename("feature. ") == "feature.md"


def test_rebuilding_cards_removes_stale_cards_only(tmp_path):
    cards_dir = tmp_path / "cards"
    cards_dir.mkdir()
    (cards_dir / "OLD_GENE.md").write_text("# Evidence Card: OLD_GENE\n")
    (cards_dir / "README.md").write_text("# Notes kept by a curator\n")
    evidence = harmonize_demo(tmp_path / "evidence.tsv")
    cards = build_evidence_cards(evidence_path=evidence, output_dir=cards_dir)
    names = sorted(path.name for path in cards_dir.glob("*.md"))
    assert names == sorted([path.name for path in cards] + ["README.md"])
    assert not any(":" in name for name in names)


def test_colliding_card_filenames_fail_before_writing(tmp_path):
    path = _evidence_rows(tmp_path, [
        {"feature_id_standardized": "NCBI:G1", "study_id": "A", "effect_size": 1.0, "standard_error": 0.1},
        {"feature_id_standardized": "ncbi_g1", "study_id": "B", "effect_size": 1.0, "standard_error": 0.1},
    ])
    with pytest.raises(ValueError, match="same card file"):
        build_evidence_cards(evidence_path=path, output_dir=tmp_path / "cards")
    assert not list((tmp_path / "cards").glob("*.md"))


def test_bulk_identifier_mapping_matches_single_lookups():
    from aree.harmonize.identifiers import load_mapping, map_identifiers

    mapping = load_mapping(ROOT / "data/mappings/cgigas_cgi_to_ncbi_gene_rs2024_06_v1.tsv")
    ids = list(mapping["feature_id_original"].head(20)) + ["NOT_A_REAL_ID", mapping["feature_id_original"].iloc[0]]
    assert map_identifiers(ids, mapping) == [map_identifier(original_id, mapping) for original_id in ids]


def test_column_groups_matches_groupby_order_and_skips_missing_keys():
    from aree.groups import column_groups

    frame = pd.DataFrame({"key": ["b", "a", None, "b", "a"], "value": [1, 2, 3, 4, 5]})
    groups = [(key, list(columns["value"])) for key, columns in column_groups(frame, "key", ["value"])]
    expected = [(key, list(group["value"])) for key, group in frame.groupby("key")]
    assert groups == expected == [("a", [2, 5]), ("b", [1, 4])]


def test_project_root_resolution(tmp_path, monkeypatch):
    import aree.paths as paths

    project = tmp_path / "project"
    (project / "registry" / "studies").mkdir(parents=True)
    nested = project / "data" / "processed"
    nested.mkdir(parents=True)
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()

    # A source checkout (editable install) wins when AREE_ROOT is unset.
    monkeypatch.delenv("AREE_ROOT", raising=False)
    assert paths.project_root() == ROOT.resolve()

    # A regular install has no project around the package: use the nearest project above cwd.
    monkeypatch.setattr(paths, "_SOURCE_CHECKOUT", tmp_path / "site-packages")
    monkeypatch.chdir(nested)
    assert paths.project_root() == project.resolve()
    monkeypatch.chdir(elsewhere)
    assert paths.project_root() == elsewhere.resolve()

    # AREE_ROOT overrides everything.
    monkeypatch.setenv("AREE_ROOT", str(project))
    assert paths.root_path("registry") == project.resolve() / "registry"


def test_schemas_ship_inside_the_package():
    from aree.paths import package_path

    assert package_path("schemas", "study.schema.json").is_file()
    assert package_path("schemas", "evidence.schema.json").is_file()
