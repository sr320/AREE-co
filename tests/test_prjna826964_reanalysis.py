import csv
import gzip

from aree.raw.fastq import write_design_samplesheet
from scripts.prepare_prjna826964_evidence import prepare
from scripts.prepare_prjna826964_fastqs import load_locations, load_manifest


MANIFEST = "data/manifests/CGIG_OA_RNASEQ_PRJNA826964_runs.tsv"
LOCATIONS = "data/manifests/CGIG_OA_RNASEQ_PRJNA826964_sra_locations.tsv"


def test_every_run_has_an_sra_location_for_checksum_verification():
    locations = load_locations(LOCATIONS)
    assert {row["run_accession"] for row in load_manifest(MANIFEST)} == set(locations)
    for record in locations.values():
        assert len(record["sra_md5"]) == 32
        assert int(record["sra_size"]) > 0


def test_run_manifest_is_timepoint_balanced_across_both_arms():
    rows = load_manifest(MANIFEST)
    assert all(row["layout"] == "paired" for row in rows)
    # Three libraries per arm at each sampling day, so duration is not confounded with treatment.
    for condition in ("control", "acidified"):
        for timepoint in ("d7", "d28", "d56"):
            cell = [row for row in rows if row["condition"] == condition and row["timepoint"] == timepoint]
            assert len(cell) == 3


def test_library_names_agree_with_biosample_assignments():
    for row in load_manifest(MANIFEST):
        prefix = "C" if row["condition"] == "control" else "O"
        assert row["library_name"] == prefix + row["replicate"][1:]
        assert len(row["fastq_1_md5"]) == len(row["fastq_2_md5"]) == 32


def test_design_sheet_carries_timepoint_for_blocking(tmp_path):
    rows = load_manifest(MANIFEST)
    design = tmp_path / "design.csv"
    write_design_samplesheet(rows, design, extra_columns=("timepoint",))
    with design.open(newline="") as handle:
        design_rows = list(csv.DictReader(handle))
    assert list(design_rows[0]) == ["sample", "condition", "replicate", "run_accession", "timepoint"]
    assert {row["timepoint"] for row in design_rows} == {"d7", "d28", "d56"}
    assert len({row["sample"] for row in design_rows}) == len(design_rows)


def test_oa_reanalysis_exports_acidification_evidence(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t-2.0\t0.25\t0.0001\t0.001\n"
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
    assert processed_rows[0]["sample_comparison"] == "ocean_acidification_vs_ambient_control"
    assert processed_rows[0]["molecular_direction"] == "down"
    assert "ocean_acidification_response" in processed_rows[0]["quality_flags"]
