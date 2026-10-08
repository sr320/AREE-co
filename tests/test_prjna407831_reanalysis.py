import csv
import gzip
import re
from collections import Counter

from aree.raw.fastq import run_mates
from scripts.prepare_prjna407831_evidence import CONTRASTS, prepare
from scripts.prepare_prjna407831_fastqs import DEFAULT_LOCATIONS, load_manifest
from scripts.prepare_prjna1329250_fastqs import load_locations


MANIFEST = "data/manifests/CGIG_HEAT_RNASEQ_PRJNA407831_runs.tsv"


def test_manifest_is_four_populations_by_three_tanks_by_three_times():
    rows = load_manifest(MANIFEST)
    assert set(Counter((row["population"], row["condition"]) for row in rows).values()) == {3}
    assert {row["population"] for row in rows} == {"BYQJ", "BYQX", "LTJ", "LTX"}
    for row in rows:
        alias = re.search(r"alias=(\w+);", row["assignment_evidence"]).group(1)
        assert alias == "{}{}{}".format(row["population"], row["condition"][1:], row["replicate"].split("_")[1])
        assert "temp=35_{}".format(row["condition"][1:]) in row["assignment_evidence"]
        assert run_mates(row) == ("1", "2")


def test_every_run_has_a_checksummed_ncbi_archive():
    locations = load_locations(DEFAULT_LOCATIONS)
    assert set(locations) == {row["run_accession"] for row in load_manifest(MANIFEST)}
    assert all(len(location["sra_md5"]) == 32 for location in locations.values())


def test_both_records_share_the_0h_dependence_group(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t5.0\t0.5\t0.0001\t0.001\n"
    )
    gff = tmp_path / "reference.gff.gz"
    with gzip.open(gff, "wt") as handle:
        handle.write(
            "chr1\tRefSeq\tgene\t1\t10\t.\t+\t.\t"
            "ID=gene-A;Dbxref=GeneID:1;Name=A;description=alpha;gene=A;gene_biotype=protein_coding\n"
        )
    for contrast, spec in CONTRASTS.items():
        processed = tmp_path / "{}.tsv".format(contrast)
        prepare(contrast, results, gff, processed, tmp_path / "mapping.tsv", tmp_path / "annotations.tsv")
        rows = list(csv.DictReader(processed.open(), delimiter="\t"))
        assert rows[0]["sample_comparison"] == spec["sample_comparison"]
        assert {"dependence_group=PRJNA407831_0h_pools", "pooled_libraries"} <= set(rows[0]["quality_flags"].split(";"))


def test_the_0h_outlier_is_left_out_of_both_design_sheets(tmp_path):
    from scripts.prepare_prjna407831_fastqs import EXCLUDED_SAMPLES, write_contrast_sheet

    rows = load_manifest(MANIFEST)
    assert EXCLUDED_SAMPLES == {"h0_LTJ_2"}
    for test in ("h6", "h24"):
        path = tmp_path / "{}.csv".format(test)
        write_contrast_sheet(rows, test, path)
        design = list(csv.DictReader(path.open()))
        assert Counter(row["condition"] for row in design) == {"control": 11, test: 12}
        assert "h0_LTJ_2" not in {row["sample"] for row in design}
