import csv
import gzip

import numpy as np

from aree.raw.fastq import run_mates
from scripts.prepare_prjna735889_evidence import prepare
from scripts.prepare_prjna735889_fastqs import DEFAULT_LOCATIONS, EXCLUDED_RUNS, load_manifest
from scripts.prepare_prjna1329250_fastqs import load_locations
from scripts.score_prjna735889_gametogenesis import module_scores, widest_gap_threshold


MANIFEST = "data/manifests/CGIG_OA_RNASEQ_PRJNA735889_runs.tsv"


def test_bands_exclude_the_tipping_window():
    rows = load_manifest(MANIFEST)
    for row in rows:
        ph = float(row["ph_total"])
        if row["condition"] == "control":
            assert 7.4 <= ph <= 7.8
        else:
            assert 6.5 <= ph <= 6.8
    # One tank per pH level (pHT 6.5 has two), five tanks per band.
    tanks = {(row["condition"], row["tank"]) for row in rows}
    assert len([tank for tank in tanks if tank[0] == "control"]) == 5
    assert len([tank for tank in tanks if tank[0] == "acidified"]) == 5
    assert len({row["tank"] for row in rows}) == 10


def test_biosample_attributes_agree_with_the_assigned_tank_and_ph():
    for row in load_manifest(MANIFEST):
        assert "title={}{};".format(row["tank"][1:], row["replicate"].split("_", 1)[1]) in row["assignment_evidence"]
        assert "treatment=pHT {}".format(row["ph_total"]) in row["assignment_evidence"]
        assert run_mates(row) == ("1", "2")


def test_every_run_has_one_ncbi_archive_and_13ebis_is_excluded():
    rows = load_manifest(MANIFEST)
    locations = load_locations(DEFAULT_LOCATIONS)
    assert set(locations) == {row["run_accession"] for row in rows}
    assert all(len(location["sra_md5"]) == 32 for location in locations.values())
    excluded = [row for row in rows if row["run_accession"] in EXCLUDED_RUNS]
    assert [row["replicate"] for row in excluded] == ["t13_Ebis"]


def test_gametogenesis_module_separates_active_oysters(tmp_path):
    rng = np.random.default_rng(0)
    samples = ["s{}".format(i) for i in range(10)]
    active = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
    base = rng.normal(8, 0.1, (20, 10))
    module = 8 + 3 * active + rng.normal(0, 0.05, (5, 10))
    counts = 2 ** np.vstack([module, base]) - 1
    normalized = tmp_path / "counts.tsv"
    with normalized.open("w") as handle:
        handle.write("feature_id_standardized\t" + "\t".join(samples) + "\n")
        for i, row in enumerate(counts):
            handle.write("g{}\t".format(i) + "\t".join("{:.3f}".format(value) for value in row) + "\n")
    pca = tmp_path / "pca.csv"
    with pca.open("w") as handle:
        handle.write("PC1,name\n")
        for name, flag in zip(samples, active):
            handle.write("{},{}\n".format(-10.0 * flag + rng.normal(0, 0.1), name))
    names, scores, size = module_scores(normalized, pca)
    assert size == 5
    called = scores > widest_gap_threshold(scores)
    assert list(called) == [bool(flag) for flag in active]


def test_ph_band_reanalysis_exports_gametogenesis_blocked_evidence(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t-6.7\t0.9\t0.0001\t0.001\n"
    )
    gff = tmp_path / "reference.gff.gz"
    with gzip.open(gff, "wt") as handle:
        handle.write(
            "chr1\tRefSeq\tgene\t1\t10\t.\t+\t.\t"
            "ID=gene-A;Dbxref=GeneID:1;Name=A;description=alpha;gene=A;gene_biotype=protein_coding\n"
        )
    processed = tmp_path / "processed.tsv"
    prepare(results, gff, processed, tmp_path / "mapping.tsv", tmp_path / "annotations.tsv")
    rows = list(csv.DictReader(processed.open(), delimiter="\t"))
    assert rows[0]["sample_comparison"] == "ph_6.5_6.8_vs_ph_7.4_7.8_23d_juveniles"
    assert rows[0]["molecular_direction"] == "down"
    assert "gametogenesis_blocked_contrast" in rows[0]["quality_flags"]
