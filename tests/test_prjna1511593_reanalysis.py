import csv
import re
from collections import Counter

from aree.raw.fastq import run_mates
from scripts.prepare_prjna1511593_fastqs import DEFAULT_LOCATIONS, load_manifest, write_contrast_sheet
from scripts.prepare_prjna1329250_fastqs import load_locations


MANIFEST = "data/manifests/CGIG_HEAT_RNASEQ_PRJNA1511593_runs.tsv"
CODE = {"control": "H0", "heat_37C_12h": "H37H12", "heat_42C_1h": "H42H1"}


def test_only_the_gigas_line_is_used_and_codes_match_sample_names():
    rows = load_manifest(MANIFEST)
    assert Counter(row["condition"] for row in rows) == {"control": 3, "heat_37C_12h": 3, "heat_42C_1h": 3}
    for row in rows:
        name = re.search(r"sample_name=(\w+)", row["assignment_evidence"]).group(1)
        assert name == "GG{}S{}".format(CODE[row["condition"]], row["replicate"])
        assert "organism Magallana gigas" in row["assignment_evidence"]
        assert run_mates(row) == ("1", "2")


def test_every_run_has_a_checksummed_ncbi_archive():
    assert set(load_locations(DEFAULT_LOCATIONS)) == {row["run_accession"] for row in load_manifest(MANIFEST)}


def test_contrast_sheets_are_three_v_three(tmp_path):
    rows = load_manifest(MANIFEST)
    for test in ("heat_37C_12h", "heat_42C_1h"):
        path = tmp_path / "{}.csv".format(test)
        write_contrast_sheet(rows, test, path)
        assert Counter(row["condition"] for row in csv.DictReader(path.open())) == {"control": 3, test: 3}

