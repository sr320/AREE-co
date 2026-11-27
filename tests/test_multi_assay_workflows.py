"""Tests for validate-processed and Nextflow multi-assay workflows (methylation, proteomics, metabolomics)."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from aree.raw.assay_evidence import validate_processed_table

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = os.environ.get("AREE_NEXTFLOW", "nextflow")


def _works(command):
    try:
        return subprocess.run(command, capture_output=True, timeout=300).returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


REQUIRED = os.environ.get("AREE_REQUIRE_PIPELINE") == "1"


def _requires(available, reason):
    if REQUIRED and not available:
        pytest.fail("AREE_REQUIRE_PIPELINE=1 but " + reason, pytrace=False)
    return pytest.mark.skipif(not available, reason=reason)


needs_nextflow = _requires(_works([NEXTFLOW, "-version"]), "usable Nextflow not found")
needs_limma = _requires(
    _works(["Rscript", "-e", "stopifnot(requireNamespace('limma', quietly=TRUE))"]),
    "R with limma is required for multi-assay differential runs",
)


def test_validate_processed_rejects_invalid_rows(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text(
        "sample_comparison\tfeature_id_original\tfeature_type\tmolecular_direction\t"
        "effect_size\teffect_size_type\tstandard_error\tp_value\tadjusted_p_value\tanalysis_method\n"
        "heat_vs_ctrl\tCPG_1\tcpg_site\tsideways\t1.2\tm_value_diff\t0.2\t0.01\t0.04\tlimma\n"
    )
    with pytest.raises(ValueError, match="invalid molecular_direction"):
        validate_processed_table(bad)


@needs_nextflow
@needs_limma
def test_methylation_workflow_end_to_end(tmp_path):
    data_dir = tmp_path / "data"
    subprocess.run(
        [sys.executable, str(ROOT / "workflows/methylation/tests/make_test_data.py"), "--output", str(data_dir)],
        check=True,
    )
    res = subprocess.run(
        [
            NEXTFLOW, "-q", "run", str(ROOT / "workflows/methylation/main.nf"), "-profile", "test",
            "--study_id", "SYNTH_METH",
            "--coverage_table", str(data_dir / "methylation_coverage.tsv"),
            "--design", str(data_dir / "design.csv"),
            "--test_level", "heat",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert res.returncode == 0, res.stdout + res.stderr
    ev_dir = tmp_path / "results/methylation/SYNTH_METH/evidence"
    proc = pd.read_csv(ev_dir / "SYNTH_METH_methylation.tsv", sep="\t")
    assert set(proc["molecular_direction"]).issubset({"hypermethylated", "hypomethylated", "unknown"})
    sig = proc[proc["adjusted_p_value"] < 0.05]
    assert len(sig) >= 30
    report = json.loads((ev_dir / "validation_report.json").read_text())
    assert report["status"] == "valid"


@needs_nextflow
@needs_limma
def test_proteomics_workflow_end_to_end(tmp_path):
    data_dir = tmp_path / "data"
    subprocess.run(
        [sys.executable, str(ROOT / "workflows/proteomics/tests/make_test_data.py"), "--output", str(data_dir)],
        check=True,
    )
    res = subprocess.run(
        [
            NEXTFLOW, "-q", "run", str(ROOT / "workflows/proteomics/main.nf"), "-profile", "test",
            "--study_id", "SYNTH_PROT",
            "--abundance_table", str(data_dir / "proteomics_abundance.tsv"),
            "--design", str(data_dir / "design.csv"),
            "--test_level", "heat",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert res.returncode == 0, res.stdout + res.stderr
    ev_dir = tmp_path / "results/proteomics/SYNTH_PROT/evidence"
    proc = pd.read_csv(ev_dir / "SYNTH_PROT_proteomics.tsv", sep="\t")
    assert set(proc["feature_type"]) == {"protein"}
    assert set(proc["molecular_direction"]).issubset({"higher", "lower", "unknown"})
    assert (proc["adjusted_p_value"] < 0.05).sum() >= 25
    assert json.loads((ev_dir / "validation_report.json").read_text())["status"] == "valid"


@needs_nextflow
@needs_limma
def test_metabolomics_workflow_end_to_end(tmp_path):
    data_dir = tmp_path / "data"
    subprocess.run(
        [sys.executable, str(ROOT / "workflows/metabolomics/tests/make_test_data.py"), "--output", str(data_dir)],
        check=True,
    )
    res = subprocess.run(
        [
            NEXTFLOW, "-q", "run", str(ROOT / "workflows/metabolomics/main.nf"), "-profile", "test",
            "--study_id", "SYNTH_METAB",
            "--abundance_table", str(data_dir / "metabolomics_abundance.tsv"),
            "--design", str(data_dir / "design.csv"),
            "--test_level", "heat",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert res.returncode == 0, res.stdout + res.stderr
    ev_dir = tmp_path / "results/metabolomics/SYNTH_METAB/evidence"
    proc = pd.read_csv(ev_dir / "SYNTH_METAB_metabolomics.tsv", sep="\t")
    assert set(proc["feature_type"]) == {"metabolite"}
    assert (proc["adjusted_p_value"] < 0.05).sum() >= 16
    assert json.loads((ev_dir / "validation_report.json").read_text())["status"] == "valid"
