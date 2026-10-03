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


def _add_design_column(design, name, values):
    with design.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row[name] = values(row)
    with design.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_several_covariates_enter_the_model_in_order(tmp_path):
    quant_dir, design, tx2gene = _write_study(tmp_path)
    # A batch that crosses both conditions and all days, like a sequencing lane.
    _add_design_column(design, "batch", lambda row: "lane" + row["replicate"][-1])
    output_dir = tmp_path / "deseq2"
    result = _run(quant_dir, design, tx2gene, output_dir, "--replicates=9", "--covariate=timepoint+batch")
    assert result.returncode == 0, result.stderr
    assert (output_dir / "design_formula.txt").read_text().strip() == "~timepoint + batch + condition"


def test_aliased_covariates_are_refused(tmp_path):
    quant_dir, design, tx2gene = _write_study(tmp_path)
    # A renamed copy of timepoint passes the per-covariate check but duplicates it in the model.
    _add_design_column(design, "day", lambda row: row["timepoint"].upper())
    result = _run(quant_dir, design, tx2gene, tmp_path / "deseq2", "--replicates=9", "--covariate=timepoint+day")
    assert result.returncode != 0
    assert "is not full rank" in result.stderr


def test_covariate_with_some_single_condition_levels_is_fitted(tmp_path):
    quant_dir, design, tx2gene = _write_study(tmp_path)
    # Like an unbalanced multiplex pool: "a" holds control d7 only, "b" acidified d7 only, and "c"
    # carries both arms, so the acidification effect is still estimable within "c".
    def pool(row):
        if row["timepoint"] == "d7":
            return "a" if row["condition"] == "control" else "b"
        return "c"
    _add_design_column(design, "pool", pool)
    result = _run(quant_dir, design, tx2gene, tmp_path / "deseq2", "--replicates=9", "--covariate=pool")
    assert result.returncode == 0, result.stderr


def test_covariate_nested_entirely_within_condition_is_refused(tmp_path):
    quant_dir, design, tx2gene = _write_study(tmp_path)
    _add_design_column(design, "tank", lambda row: row["condition"] + "_tank")
    result = _run(quant_dir, design, tx2gene, tmp_path / "deseq2", "--replicates=9", "--covariate=tank")
    assert result.returncode != 0
    assert "is confounded with condition" in result.stderr


def test_tag_seq_drops_the_transcript_length_offset(tmp_path):
    """With --tag-seq=yes, normalization is one size factor per sample, whatever Salmon's lengths."""
    rng = np.random.default_rng(11)
    genes = 200
    quant_dir = tmp_path / "salmon"
    tx2gene = tmp_path / "tx2gene.tsv"
    tx2gene.write_text("transcript_id\tgene_id\n" + "".join(f"tx{g}\tgene{g}\n" for g in range(genes)))
    rows, raw = [], {}
    for condition in ("control", "acidified"):
        for replicate in (1, 2, 3):
            sample = f"{condition}_{replicate}"
            counts = rng.poisson(rng.gamma(2.0, 100.0, genes)) + 10
            # Effective lengths that differ by sample would enter the tximport offset.
            lengths = rng.uniform(200, 3000, genes)
            sample_dir = quant_dir / sample
            sample_dir.mkdir(parents=True)
            with (sample_dir / "quant.sf").open("w") as handle:
                handle.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
                for gene, (count, length) in enumerate(zip(counts, lengths)):
                    handle.write(f"tx{gene}\t{length + 200:.0f}\t{length:.1f}\t{count / length:.4f}\t{count}\n")
            raw[sample] = counts
            rows.append({"sample": sample, "condition": condition, "replicate": str(replicate),
                         "run_accession": f"SRR{len(rows)}"})
    design = tmp_path / "design.csv"
    with design.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    output_dir = tmp_path / "deseq2"
    result = _run(quant_dir, design, tx2gene, output_dir, "--replicates=3", "--tag-seq=yes")
    assert result.returncode == 0, result.stderr
    with (output_dir / "acidified_vs_control_normalized_counts.tsv").open() as handle:
        table = list(csv.reader(handle, delimiter="\t"))
    header = table[0]
    for column, sample in enumerate(header[1:], start=1):
        ratios = [float(row[column]) / raw[sample][int(row[0].split("gene")[-1])] for row in table[1:]]
        assert max(ratios) / min(ratios) == pytest.approx(1.0, abs=1e-6)


def test_tag_seq_rejects_other_values(tmp_path):
    quant_dir, design, tx2gene = _write_study(tmp_path)
    result = _run(quant_dir, design, tx2gene, tmp_path / "out", "--tag-seq=maybe")
    assert result.returncode != 0
    assert "--tag-seq must be yes or no" in result.stderr
