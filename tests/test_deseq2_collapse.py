import csv
import shutil
import subprocess

import numpy as np
import pytest


SCRIPT = "scripts/run_salmon_tximport_deseq2.R"


def _has_deseq2():
    if shutil.which("Rscript") is None:
        return False
    probe = subprocess.run(
        ["Rscript", "-e", "stopifnot(requireNamespace('DESeq2', quietly=TRUE), requireNamespace('tximport', quietly=TRUE))"],
        capture_output=True,
    )
    return probe.returncode == 0


pytestmark = pytest.mark.skipif(not _has_deseq2(), reason="Rscript with DESeq2 and tximport is required")


def _write_study(tmp_path):
    """Three pools per arm, each sequenced as three libraries of the same RNA."""
    rng = np.random.default_rng(7)
    genes = 300
    quant_dir = tmp_path / "salmon"
    tx2gene = tmp_path / "tx2gene.tsv"
    tx2gene.write_text("transcript_id\tgene_id\n" + "".join(f"tx{g}\tgene{g}\n" for g in range(genes)))
    design_rows = []
    for condition in ("control", "acidified"):
        for timepoint in ("d7", "d28", "d56"):
            pool_mean = rng.gamma(2.0, 100.0, genes)
            for replicate in (1, 2, 3):
                sample = f"{condition}_{timepoint}_{replicate}"
                counts = rng.poisson(pool_mean)
                sample_dir = quant_dir / sample
                sample_dir.mkdir(parents=True)
                with (sample_dir / "quant.sf").open("w") as handle:
                    handle.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
                    for gene, count in enumerate(counts):
                        handle.write(f"tx{gene}\t1500\t1300.0\t{count / 13.0:.4f}\t{count}\n")
                design_rows.append({
                    "sample": sample, "condition": condition, "replicate": f"{timepoint}_{replicate}",
                    "run_accession": f"SRR{len(design_rows)}", "timepoint": timepoint,
                })
    design = tmp_path / "design.csv"
    with design.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(design_rows[0]))
        writer.writeheader()
        writer.writerows(design_rows)
    return quant_dir, design, tx2gene


def _run(quant_dir, design, tx2gene, output_dir, *options):
    return subprocess.run(
        ["Rscript", SCRIPT, str(quant_dir), str(design), str(tx2gene), str(output_dir), "--test=acidified", *options],
        capture_output=True, text=True,
    )


def test_collapse_sums_technical_replicates_into_one_sample_per_pool(tmp_path):
    quant_dir, design, tx2gene = _write_study(tmp_path)
    output_dir = tmp_path / "deseq2"
    result = _run(quant_dir, design, tx2gene, output_dir,
                  "--replicates=3", "--covariate=timepoint", "--collapse-by=condition+timepoint")
    assert result.returncode == 0, result.stderr

    with (output_dir / "acidified_vs_control_normalized_counts.tsv").open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
    assert header[1:] == ["control_d7", "control_d28", "control_d56", "acidified_d7", "acidified_d28", "acidified_d56"]

    with (output_dir / "collapsed_pools.csv").open(newline="") as handle:
        pools = list(csv.DictReader(handle))
    assert len(pools) == 18
    assert {row["pool"] for row in pools if row["library"].startswith("control_d7_")} == {"control_d7"}

    with (output_dir / "sample_deseq2_qc.tsv").open() as handle:
        assert len(handle.read().strip().splitlines()) == 1 + 6


def test_uncollapsed_replicate_count_is_still_checked_per_library(tmp_path):
    quant_dir, design, tx2gene = _write_study(tmp_path)
    # Without collapsing there are nine libraries per arm, so requiring three must fail.
    result = _run(quant_dir, design, tx2gene, tmp_path / "deseq2", "--replicates=3", "--covariate=timepoint")
    assert result.returncode != 0
    assert "expected 3 control and 3 acidified samples" in result.stderr


def test_collapse_refuses_a_pool_that_spans_covariate_levels(tmp_path):
    # Pooling on condition alone would merge libraries from different sampling days.
    quant_dir, design, tx2gene = _write_study(tmp_path)
    result = _run(quant_dir, design, tx2gene, tmp_path / "deseq2",
                  "--replicates=1", "--covariate=timepoint", "--collapse-by=condition")
    assert result.returncode != 0
    assert "spans more than one condition or covariate level" in result.stderr
