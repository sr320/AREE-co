import csv
import gzip
import re
from collections import Counter

from aree.raw.fastq import run_mates
from scripts.prepare_prjna852380_evidence import CONTRASTS, prepare
from scripts.prepare_prjna852380_fastqs import DEFAULT_LOCATIONS, load_manifest, write_contrast_sheet
from scripts.prepare_prjna1329250_fastqs import load_locations


MANIFEST = "data/manifests/CGIG_VIBRIO_RNASEQ_PRJNA852380_runs.tsv"


def test_groups_match_sample_names():
    rows = load_manifest(MANIFEST)
    assert set(Counter(row["condition"] for row in rows).values()) == {3}
    for row in rows:
        name = re.search(r"sample_name=(\w+);", row["assignment_evidence"]).group(1)
        expected = "control{}".format(row["replicate"]) if row["condition"] == "control" else \
            "V{}_{}".format(row["condition"].split("_")[1], row["replicate"])
        assert name == expected
        assert run_mates(row) == ("1", "2")


def test_every_run_has_a_checksummed_ncbi_archive():
    assert set(load_locations(DEFAULT_LOCATIONS)) == {row["run_accession"] for row in load_manifest(MANIFEST)}


def test_only_12h_and_24h_are_exported_against_shared_controls(tmp_path):
    assert set(CONTRASTS) == {"vibrio_12h", "vibrio_24h"}
    rows = load_manifest(MANIFEST)
    for test in CONTRASTS:
        path = tmp_path / "{}.csv".format(test)
        write_contrast_sheet(rows, test, path)
        assert Counter(row["condition"] for row in csv.DictReader(path.open())) == {"control": 3, test: 3}
    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t4.0\t0.5\t0.0001\t0.001\n"
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
        row = list(csv.DictReader(processed.open(), delimiter="\t"))[0]
        assert row["sample_comparison"] == spec["sample_comparison"]
        assert {"dependence_group=PRJNA852380_controls", "vibrio_species_not_stated"} <= set(row["quality_flags"].split(";"))
