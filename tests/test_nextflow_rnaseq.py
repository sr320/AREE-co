"""Run workflows/rnaseq on synthetic data. Skipped unless Nextflow (and, for the full run, the bio tools) exist.

Set AREE_NEXTFLOW to a Nextflow launcher if the one on PATH is not usable, and
AREE_REQUIRE_PIPELINE=1 (as the CI pipeline job does) to fail instead of skipping when tools are missing.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
PIPELINE = ROOT / "workflows" / "rnaseq"
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
needs_bio_tools = _requires(
    bool(shutil.which("salmon"))
    and bool(shutil.which("fastqc"))
    and _works(["Rscript", "-e", "stopifnot(requireNamespace('DESeq2', quietly=TRUE), "
                                 "requireNamespace('tximport', quietly=TRUE))"]),
    "salmon, fastqc and R with DESeq2/tximport are required for the full pipeline run",
)


@pytest.fixture(scope="module")
def study(tmp_path_factory):
    data = tmp_path_factory.mktemp("synthetic_study")
    subprocess.run([sys.executable, str(PIPELINE / "tests" / "make_test_data.py"), "--output", str(data)], check=True)
    return data


def run_pipeline(study, workdir, *extra):
    command = [
        NEXTFLOW, "-q", "run", str(PIPELINE), "-profile", "test",
        "--study_id", "SIM_HEAT", "--test_level", "heat",
        "--samplesheet", str(study / "samplesheet.csv"), "--design", str(study / "design.csv"),
        "--tx2gene", str(study / "tx2gene.tsv"), "--gff", str(study / "genomic.gff.gz"),
        "--transcripts", str(study / "transcripts.fa"), "--quality_flags", "raw_reanalysis,synthetic_test_data",
        *extra,
    ]
    result = subprocess.run(command, cwd=workdir, capture_output=True, text=True, timeout=1800)
    assert result.returncode == 0, result.stdout + result.stderr
    return workdir / "results" / "rnaseq" / "SIM_HEAT"


@needs_nextflow
def test_pipeline_wiring_with_stub_run(study, tmp_path):
    results = run_pipeline(study, tmp_path, "-stub-run")
    for name in ("SIM_HEAT_rnaseq.tsv", "SIM_HEAT_mapping.tsv", "SIM_HEAT_annotations.tsv"):
        assert (results / "evidence" / name).exists()
    assert (results / "deseq2" / "heat_vs_control_deseq2_all_genes.tsv").exists()
    assert len(list((results / "salmon").glob("*/quant.sf"))) == 6


@needs_nextflow
@needs_bio_tools
def test_pipeline_recovers_simulated_expression_changes(study, tmp_path):
    results = run_pipeline(study, tmp_path)
    table = pd.read_csv(results / "deseq2" / "heat_vs_control_deseq2_all_genes.tsv", sep="\t")
    index = table["feature_id_standardized"].str.rsplit(":", n=1).str[-1].astype(int) - 105500000
    significant = table["padj"] < 0.05
    up = set(index[significant & (table["log2FoldChange"] > 0)])
    down = set(index[significant & (table["log2FoldChange"] < 0)])
    # make_test_data.py simulates genes 0-29 at 4x and 30-59 at 0.25x under heat.
    assert len(up & set(range(0, 30))) >= 27
    assert len(down & set(range(30, 60))) >= 27
    assert not (up | down) - set(range(0, 60))

    qc = pd.read_csv(results / "salmon" / "salmon_qc_summary.tsv", sep="\t")
    assert list(qc["sample"]) == ["control_A", "control_B", "control_C", "heat_A", "heat_B", "heat_C"]
    assert (qc["mapping_rate_percent"] > 90).all()

    from aree.harmonize.processed import harmonize_processed

    evidence = harmonize_processed(
        "CGIG_HEAT_RNASEQ_PRJNA516762",
        results / "evidence" / "SIM_HEAT_rnaseq.tsv",
        output_path=tmp_path / "harmonized.tsv",
        mapping_path=results / "evidence" / "SIM_HEAT_mapping.tsv",
    )
    harmonized = pd.read_csv(evidence, sep="\t")
    assert len(harmonized) == len(table)
    assert set(harmonized["mapping_confidence"]) == {"exact"}
    # analysis_method records the tool versions that actually ran, not a hard-coded default.
    method = (results / "deseq2" / "analysis_method.txt").read_text().strip()
    salmon_version = subprocess.run(["salmon", "--version"], capture_output=True, text=True).stdout.split()[-1]
    assert method.startswith("Salmon_{}_tximport_".format(salmon_version)) and "_DESeq2_" in method
    assert set(harmonized["analysis_method"]) == {method}
