import csv
import gzip
from collections import Counter

from aree.raw.fastq import run_mates
from scripts.prepare_prjna856813_evidence import CONTRASTS, prepare
from scripts.prepare_prjna856813_fastqs import DEFAULT_LOCATIONS, load_manifest
from scripts.prepare_prjna1329250_fastqs import load_locations


MANIFEST = "data/manifests/CGIG_TIRE_RNASEQ_PRJNA856813_runs.tsv"
CODE = {"control": "CTL", "leachate_high": "LEA_H", "tire_particles_high": "MR_H"}


def test_manifest_is_day41_gill_six_per_arm_and_matches_sample_names():
    rows = load_manifest(MANIFEST)
    assert set(Counter(row["condition"] for row in rows).values()) == {6}
    for row in rows:
        assert run_mates(row) == ("1",)
        assert "sample_name=D41_GI_{}_{};".format(CODE[row["condition"]], row["replicate"]) in row["assignment_evidence"]


def test_every_run_has_a_checksummed_ncbi_archive():
    locations = load_locations(DEFAULT_LOCATIONS)
    assert set(locations) == {row["run_accession"] for row in load_manifest(MANIFEST)}
    assert all(len(location["sra_md5"]) == 32 for location in locations.values())


def test_both_records_share_the_gill_control_dependence_group(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t1.5\t0.3\t0.0001\t0.001\n"
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
        assert "dependence_group=PRJNA856813_gill_controls" in rows[0]["quality_flags"].split(";")
