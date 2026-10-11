import csv
import re
from collections import Counter

from aree.raw.fastq import run_mates
from scripts.prepare_prjna1244231_fastqs import DEFAULT_LOCATIONS, load_manifest, write_contrast_sheet
from scripts.prepare_prjna1329250_fastqs import load_locations


MANIFEST = "data/manifests/CGIG_ERK_HEAT_RNASEQ_PRJNA1244231_runs.tsv"
CODE = {"sw_0h": "SW-0", "sw_12h": "SW-12", "erk_0h": "ERK-0", "erk_12h": "ERK-12"}


def test_conditions_match_sample_names():
    rows = load_manifest(MANIFEST)
    assert Counter(row["condition"] for row in rows) == {"sw_0h": 3, "sw_12h": 3, "erk_0h": 3, "erk_12h": 3}
    for row in rows:
        name = re.search(r"sample_name=([\w-]+)", row["assignment_evidence"]).group(1)
        assert name == "{}-{}".format(CODE[row["condition"]], row["replicate"])
        assert run_mates(row) == ("1", "2")


def test_every_run_has_a_checksummed_ncbi_archive():
    assert set(load_locations(DEFAULT_LOCATIONS)) == {row["run_accession"] for row in load_manifest(MANIFEST)}


def test_heat_contrast_is_seawater_only_three_v_three(tmp_path):
    path = tmp_path / "sw_12h.csv"
    write_contrast_sheet(load_manifest(MANIFEST), "sw_12h", path)
    rows = list(csv.DictReader(path.open()))
    assert Counter(row["condition"] for row in rows) == {"control": 3, "sw_12h": 3}
    assert all(row["sample"].startswith("sw_") for row in rows)
