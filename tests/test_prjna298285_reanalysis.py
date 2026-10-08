import csv
import gzip
import re
from collections import Counter

from aree.raw.fastq import run_mates
from scripts.prepare_prjna298285_evidence import CONTRASTS, prepare
from scripts.prepare_prjna298285_fastqs import DEFAULT_LOCATIONS, load_manifest, write_contrast_sheet
from scripts.prepare_prjna1329250_fastqs import load_locations


MANIFEST = "data/manifests/CGIG_LARVAL_OA_TEMP_RNASEQ_PRJNA298285_runs.tsv"
CONDITIONS = {"control": ("pH 8.1", "20 °C"), "acidified": ("pH 7.9", "20 °C"),
              "warmed": ("pH 8.1", "22 °C"), "acidified_warmed": ("pH 7.9", "22 °C")}


def test_manifest_is_two_libraries_per_treatment_and_stage():
    rows = load_manifest(MANIFEST)
    assert set(Counter((row["condition"], row["stage"]) for row in rows).values()) == {2}
    for row in rows:
        ph, temp = CONDITIONS[row["condition"]]
        assert "{}, {}".format(ph, temp) in row["assignment_evidence"]
        assert run_mates(row) == ("1", "2")


def test_the_unnamed_library_is_the_only_inferred_replicate():
    inferred = [row for row in load_manifest(MANIFEST) if "inferred" in row["assignment_evidence"]]
    assert [(row["run_accession"], row["replicate"]) for row in inferred] == [
        (inferred[0]["run_accession"], "late_veliger_B")]
    assert inferred[0]["library_name"].startswith("C3759")
    for row in load_manifest(MANIFEST):
        if row not in inferred:
            assert re.search(r"(ambi|acid|warm|intr){}$".format(row["replicate"][-1]), row["library_name"])


def test_every_run_has_a_checksummed_ncbi_archive():
    locations = load_locations(DEFAULT_LOCATIONS)
    assert set(locations) == {row["run_accession"] for row in load_manifest(MANIFEST)}


def test_each_contrast_sheet_is_six_v_six(tmp_path):
    rows = load_manifest(MANIFEST)
    for test in ("acidified", "warmed", "acidified_warmed"):
        path = tmp_path / "{}.csv".format(test)
        write_contrast_sheet(rows, test, path)
        assert Counter(row["condition"] for row in csv.DictReader(path.open())) == {"control": 6, test: 6}


def test_records_share_the_ambient_dependence_group(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t1.0\t0.3\t0.0001\t0.001\n"
    )
    gff = tmp_path / "reference.gff.gz"
    with gzip.open(gff, "wt") as handle:
        handle.write(
            "chr1\tRefSeq\tgene\t1\t10\t.\t+\t.\t"
            "ID=gene-A;Dbxref=GeneID:1;Name=A;description=alpha;gene=A;gene_biotype=protein_coding\n"
        )
    # Only acidification is exported; the warming contrasts are confounded with developmental rate.
    assert [spec["study_id"] for spec in CONTRASTS.values()] == ["CGIG_LARVAL_OA_RNASEQ_PRJNA298285"]
    for contrast, spec in CONTRASTS.items():
        processed = tmp_path / "{}.tsv".format(contrast)
        prepare(contrast, results, gff, processed, tmp_path / "mapping.tsv", tmp_path / "annotations.tsv")
        rows = list(csv.DictReader(processed.open(), delimiter="\t"))
        assert rows[0]["sample_comparison"] == spec["sample_comparison"]
        assert "dependence_group=PRJNA298285_ambient_larvae" in rows[0]["quality_flags"].split(";")
