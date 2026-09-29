import csv
import gzip

from aree.raw.fastq import write_design_samplesheet
from scripts.prepare_prjna756710_evidence import prepare
from scripts.prepare_prjna756710_fastqs import load_manifest


MANIFEST = "data/manifests/CGIG_SALINITY_RNASEQ_PRJNA756710_runs.tsv"


def test_run_manifest_is_timepoint_balanced_across_both_arms():
    rows = load_manifest(MANIFEST)
    assert all(row["layout"] == "paired" for row in rows)
    # Three libraries per arm at each time, so time since infection is not confounded with salinity.
    for condition in ("control", "low_salinity"):
        for timepoint in ("h12", "h48"):
            cell = [row for row in rows if row["condition"] == condition and row["timepoint"] == timepoint]
            assert len(cell) == 3


def test_shared_zero_hour_baselines_are_not_in_the_contrast():
    runs = {row["run_accession"] for row in load_manifest(MANIFEST)}
    assert not runs & {"SRR15559847", "SRR15559848"}


def test_library_names_agree_with_biosample_assignments():
    for row in load_manifest(MANIFEST):
        # C10 is the 10 ppt arm and C23 the 30 ppt arm; suffixes 1-3 are 12 h and 4-6 are 48 h.
        prefix = "C10" if row["condition"] == "low_salinity" else "C23"
        timepoint, replicate = row["replicate"].split("_")
        offset = 0 if timepoint == "h12" else 3
        assert row["library_name"] == "{}-{}".format(prefix, int(replicate) + offset)
        assert len(row["fastq_1_md5"]) == len(row["fastq_2_md5"]) == 32


def test_design_sheet_carries_timepoint_for_blocking(tmp_path):
    rows = load_manifest(MANIFEST)
    design = tmp_path / "design.csv"
    write_design_samplesheet(rows, design, extra_columns=("timepoint",))
    with design.open(newline="") as handle:
        design_rows = list(csv.DictReader(handle))
    assert list(design_rows[0]) == ["sample", "condition", "replicate", "run_accession", "timepoint"]
    assert {row["timepoint"] for row in design_rows} == {"h12", "h48"}
    assert len({row["sample"] for row in design_rows}) == len(design_rows)


def test_salinity_reanalysis_exports_low_salinity_evidence(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t2.0\t0.25\t0.0001\t0.001\n"
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
    assert processed_rows[0]["sample_comparison"] == "low_salinity_10ppt_vs_30ppt_during_vibrio_infection"
    assert processed_rows[0]["molecular_direction"] == "up"
    # Every library was infected, so the evidence must say the salinity effect is conditional on it.
    assert "infection_background_all_arms" in processed_rows[0]["quality_flags"]
