import csv
import gzip

from aree.raw.fastq import run_mates, write_design_samplesheet
from scripts.prepare_prjeb18614_evidence import prepare
from scripts.prepare_prjeb18614_fastqs import load_manifest


MANIFEST = "data/manifests/CGIG_HAB_RNASEQ_PRJEB18614_runs.tsv"


def test_both_diets_are_sampled_once_at_the_same_thirteen_times():
    rows = load_manifest(MANIFEST)
    times = {condition: sorted(row["timepoint"] for row in rows if row["condition"] == condition)
             for condition in ("control", "alexandrium")}
    assert times["control"] == times["alexandrium"]
    assert len(set(times["control"])) == 13


def test_one_paired_lane_per_library_the_lower_of_its_pair():
    for row in load_manifest(MANIFEST):
        assert run_mates(row) == ("1", "2")
        first, second = (int(lane[1:]) for lane in row["sequencing_pool"].split("-"))
        assert second == first + 1
        assert row["lane"] == "L{:03d}".format(first)
        assert row["library_name"].endswith("_" + row["lane"])
        # The skipped duplicate lane is recorded for provenance.
        assert "not used" in row["assignment_evidence"]


def test_sequencing_pools_leave_the_diet_effect_estimable():
    rows = load_manifest(MANIFEST)
    arms = {}
    for row in rows:
        arms.setdefault(row["sequencing_pool"], set()).add(row["condition"])
    # Unbalanced pools are tolerated only while some pool still carries both diets.
    assert any(conditions == {"control", "alexandrium"} for conditions in arms.values())


def test_design_sheet_carries_timepoint_and_sequencing_pool(tmp_path):
    rows = load_manifest(MANIFEST)
    design = tmp_path / "design.csv"
    write_design_samplesheet(rows, design, extra_columns=("timepoint", "sequencing_pool"))
    with design.open(newline="") as handle:
        design_rows = list(csv.DictReader(handle))
    assert list(design_rows[0]) == ["sample", "condition", "replicate", "run_accession", "timepoint", "sequencing_pool"]
    assert len({row["sample"] for row in design_rows}) == 26


def test_alexandrium_reanalysis_exports_flagged_evidence(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t1.5\t0.25\t0.0001\t0.001\n"
    )
    gff = tmp_path / "reference.gff.gz"
    with gzip.open(gff, "wt") as handle:
        handle.write(
            "chr1\tRefSeq\tgene\t1\t10\t.\t+\t.\t"
            "ID=gene-A;Dbxref=GeneID:1;Name=A;description=alpha;gene=A;gene_biotype=protein_coding\n"
        )
    processed = tmp_path / "processed.tsv"
    prepare(results, gff, processed, tmp_path / "mapping.tsv", tmp_path / "annotations.tsv")
    processed_rows = list(csv.DictReader(processed.open(), delimiter="\t"))
    assert processed_rows[0]["sample_comparison"] == "alexandrium_minutum_vs_heterocapsa_triquetra_diet"
    assert processed_rows[0]["molecular_direction"] == "up"
    flags = processed_rows[0]["quality_flags"]
    assert "time_averaged_over_diurnal_cycle" in flags and "sequencing_pool_partly_confounded" in flags
