import csv
import gzip
import re
from collections import Counter

from aree.raw.fastq import run_mates
from scripts.prepare_prjeb58194_evidence import CONTRASTS as EVIDENCE
from scripts.prepare_prjeb58194_evidence import prepare
from scripts.prepare_prjeb58194_fastqs import CONTRASTS, DEFAULT_LOCATIONS, load_manifest, write_contrast_sheet
from scripts.prepare_prjna1329250_fastqs import load_locations


MANIFEST = "data/manifests/CGIG_PESTO_RNASEQ_PRJEB58194_runs.tsv"


def test_f0_groups_match_their_sample_titles():
    for row in load_manifest(MANIFEST):
        if row["generation"] != "F0":
            continue
        title = re.search(r"sample title=(\w+)", row["assignment_evidence"]).group(1)
        code = {"f0_control": "T", "f0_exposed": "E"}[row["condition"]]
        stage = {"gastrula": "g", "d_larva": "d"}[row["stage"]]
        assert title == "B{}{}{}".format(code, row["replicate"].split("_")[-1 if "ERR" not in row["replicate"] else -2], stage)
        assert run_mates(row) == ("1", "2")


def test_f1_groups_match_their_sample_titles():
    for row in load_manifest(MANIFEST):
        if row["generation"] == "F1":
            stage = {"gastrula": "Gastrula", "pediveliger": "Metamorphosis"}[row["stage"]]
            assert "{}_RNA_F1_{}".format(stage, row["condition"][3:]) in row["assignment_evidence"]


def test_every_run_has_a_checksummed_ncbi_archive_and_spot_count():
    locations = load_locations(DEFAULT_LOCATIONS)
    rows = load_manifest(MANIFEST)
    assert set(locations) == {row["run_accession"] for row in rows}
    assert all(int(row["read_count"]) > 0 for row in rows)


def test_contrast_sheets(tmp_path):
    rows = load_manifest(MANIFEST)
    expected = {"f0_exposed": 18, "f1_direct": 6, "f1_parental": 6}
    for contrast in CONTRASTS:
        path = tmp_path / "{}.csv".format(contrast)
        write_contrast_sheet(rows, contrast, path)
        design = list(csv.DictReader(path.open()))
        assert Counter(row["condition"] for row in design) == {"control": expected[contrast], contrast: expected[contrast]}
        assert len({row["stage"] for row in design}) == 2


def test_only_f1_records_share_the_tt_dependence_group(tmp_path):
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
    for contrast in EVIDENCE:
        processed = tmp_path / "{}.tsv".format(contrast)
        prepare(contrast, results, gff, processed, tmp_path / "mapping.tsv", tmp_path / "annotations.tsv")
        flags = list(csv.DictReader(processed.open(), delimiter="\t"))[0]["quality_flags"].split(";")
        assert ("dependence_group=PRJEB58194_F1_TT" in flags) == contrast.startswith("f1")
