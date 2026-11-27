import csv
import gzip
import hashlib
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

import pytest

from aree.raw.fastq import (
    download, download_once, download_runs, load_run_manifest, run_mates, total_fastq_bytes,
    write_nfcore_samplesheet,
)
from aree.raw.gff import gene_annotations, parse_attributes


class RangeHandler(SimpleHTTPRequestHandler):
    """Static file server that honors single "bytes=N-" Range requests, like ENA's FTP-over-HTTP."""

    def send_head(self):
        requested = self.headers.get("Range")
        if not requested:
            return super().send_head()
        path = self.translate_path(self.path)
        data = open(path, "rb").read()
        start = int(requested.split("=")[1].rstrip("-"))
        self.send_response(206)
        self.send_header("Content-Length", str(len(data) - start))
        self.end_headers()
        self.wfile.write(data[start:])
        return None

    def log_message(self, *args):
        pass


@pytest.fixture
def server(tmp_path):
    served = tmp_path / "served"
    served.mkdir()
    payload = b"@read\nACGT\n+\nIIII\n" * 5000
    (served / "run_1.fastq.gz").write_bytes(payload)
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), partial(RangeHandler, directory=str(served)))
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    url = "http://127.0.0.1:{}/run_1.fastq.gz".format(httpd.server_address[1])
    yield url, payload, hashlib.md5(payload).hexdigest()
    httpd.shutdown()


def test_fresh_download_is_verified_and_promoted(server, tmp_path):
    url, payload, md5 = server
    destination = tmp_path / "out" / "run_1.fastq.gz"
    destination.parent.mkdir()
    download_once(url, destination, len(payload), md5)
    assert destination.read_bytes() == payload
    assert not destination.with_suffix(".gz.part").exists()


def test_partial_download_resumes_with_range_request(server, tmp_path, capsys):
    url, payload, md5 = server
    destination = tmp_path / "run_1.fastq.gz"
    destination.with_suffix(".gz.part").write_bytes(payload[:1000])
    download_once(url, destination, len(payload), md5)
    assert destination.read_bytes() == payload
    assert "resuming" in capsys.readouterr().out


def test_checksum_mismatch_is_never_promoted(server, tmp_path):
    url, payload, _ = server
    destination = tmp_path / "run_1.fastq.gz"
    with pytest.raises(ValueError, match="MD5"):
        download_once(url, destination, len(payload), "0" * 32)
    assert not destination.exists()


def test_existing_file_must_match_checksum(tmp_path):
    destination = tmp_path / "run_1.fastq.gz"
    destination.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="Existing file fails"):
        download_once("http://unused.invalid/", destination, 7, "0" * 32)


def test_network_errors_are_retried_then_raised(tmp_path, monkeypatch):
    sleeps = []
    monkeypatch.setattr("aree.raw.fastq.time.sleep", sleeps.append)
    with pytest.raises(OSError):
        download("http://127.0.0.1:9/unreachable", tmp_path / "x.fastq.gz", 1, "0" * 32, retries=2, timeout=1)
    assert sleeps == [1, 2]


def test_run_manifest_condition_counts_are_enforced(tmp_path):
    manifest = tmp_path / "runs.tsv"
    manifest.write_text("run_accession\tcondition\nSRR1\tcontrol\nSRR2\theat\n")
    assert len(load_run_manifest(manifest, {"control": 1, "heat": 1})) == 2
    with pytest.raises(ValueError, match="Expected runs per condition"):
        load_run_manifest(manifest, {"control": 3, "heat": 3})


def test_gff_attributes_and_gene_annotations(tmp_path):
    assert parse_attributes("ID=gene-1;description=heat%20shock%3B%20protein;flag\n") == {
        "ID": "gene-1",
        "description": "heat shock; protein",
    }
    gff = tmp_path / "ref.gff.gz"
    lines = [
        "##gff-version 3",
        "NC_1\tRefSeq\tgene\t1\t9\t.\t+\t.\tDbxref=GeneID:11;gene=hsp70;gene_biotype=protein_coding",
        "NC_1\tRefSeq\tpseudogene\t1\t9\t.\t+\t.\tDbxref=GeneID:12,Other:3;Name=LOC12",
        "NC_1\tRefSeq\tmRNA\t1\t9\t.\t+\t.\tDbxref=GeneID:13",
    ]
    with gzip.open(gff, "wt") as handle:
        handle.write("\n".join(lines) + "\n")
    annotations = gene_annotations(gff, "REL")
    assert sorted(annotations) == ["NCBI:GeneID:11", "NCBI:GeneID:12"]
    assert annotations["NCBI:GeneID:11"]["gene_symbol"] == "hsp70"
    assert annotations["NCBI:GeneID:12"]["gene_symbol"] == "LOC12"


def test_single_end_runs_download_one_file_and_leave_fastq_2_empty(server, tmp_path):
    url, payload, digest = server
    row = {
        "run_accession": "run", "condition": "control", "replicate": "1", "layout": "single",
        "fastq_1": url, "fastq_1_bytes": str(len(payload)), "fastq_1_md5": digest,
        "fastq_2": "", "fastq_2_bytes": "", "fastq_2_md5": "",
    }
    assert run_mates(row) == ("1",)
    assert total_fastq_bytes([row]) == len(payload)
    download_runs([row], tmp_path / "fastq")
    assert sorted(path.name for path in (tmp_path / "fastq").iterdir()) == ["run_1.fastq.gz"]
    sheet = tmp_path / "samplesheet.csv"
    write_nfcore_samplesheet([row], tmp_path / "fastq", sheet)
    with sheet.open(newline="") as handle:
        written = next(csv.DictReader(handle))
    assert written["fastq_1"].endswith("run_1.fastq.gz")
    assert written["fastq_2"] == ""


def test_paired_end_is_the_default_and_unknown_layouts_are_rejected():
    assert run_mates({"run_accession": "run"}) == ("1", "2")
    assert run_mates({"run_accession": "run", "layout": "paired"}) == ("1", "2")
    with pytest.raises(ValueError, match="Unknown library layout"):
        run_mates({"run_accession": "run", "layout": "SINGLE"})


def test_sra_helpers_verify_and_discard(tmp_path):
    from aree.raw.sra import discard, load_locations, sra_path, verify_sra

    loc_path = tmp_path / "locations.tsv"
    payload = b"sra-bytes"
    md5 = hashlib.md5(payload).hexdigest()
    loc_path.write_text(
        "run_accession\tsra_size\tsra_md5\n"
        f"SRR1\t{len(payload)}\t{md5}\n"
        f"SRR2\t{len(payload)}\t{'0' * 32}\n"
        f"SRR3\t{len(payload)}\t{md5}\n"
    )
    locations = load_locations(loc_path)
    assert set(locations) == {"SRR1", "SRR2", "SRR3"}

    sra_dir = tmp_path / "sra"
    sra_dir.mkdir()
    sra_path(sra_dir, "SRR1").write_bytes(payload)
    sra_path(sra_dir, "SRR2").write_bytes(payload)

    rows = [{"run_accession": r} for r in ("SRR1", "SRR2", "SRR3")]
    assert sorted(verify_sra(rows, locations, sra_dir, workers=2)) == ["SRR2", "SRR3"]
    assert verify_sra([{"run_accession": "SRR1"}], locations, sra_dir, workers=1) == []

    discard([sra_path(sra_dir, "SRR1"), sra_path(sra_dir, "MISSING")])
    assert not sra_path(sra_dir, "SRR1").exists()


def test_merge_annotations_and_clean_r_na_strings(tmp_path):
    from aree.raw.deseq2_evidence import export_deseq2_evidence, merge_annotations

    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    header = "feature_id_standardized\tgene_symbol\tdescription\tgene_biotype\tannotation_release\n"
    a.write_text(header + "NCBI:GeneID:2\tG2\tbeta\tprotein_coding\tREL\n")
    b.write_text(header + "NCBI:GeneID:1\tG1\talpha\tprotein_coding\tREL\n"
                 + "NCBI:GeneID:2\tG2\tbeta_dup\tprotein_coding\tREL\n")
    merged = tmp_path / "nested" / "merged.tsv"
    merge_annotations([a, b], merged)
    lines = merged.read_text().splitlines()
    assert lines == [
        header.strip(),
        "NCBI:GeneID:1\tG1\talpha\tprotein_coding\tREL",
        "NCBI:GeneID:2\tG2\tbeta\tprotein_coding\tREL",
    ]

    bad = tmp_path / "bad.tsv"
    bad.write_text("other_col\n1\n")
    with pytest.raises(ValueError, match="different columns"):
        merge_annotations([a, bad], merged)

    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t0.5\tNA\tNA\tNA\n"
    )
    gff = tmp_path / "ref.gff.gz"
    with gzip.open(gff, "wt") as handle:
        handle.write("NC_1\tRefSeq\tgene\t1\t9\t.\t+\t.\tDbxref=GeneID:1;gene=G1;gene_biotype=protein_coding\n")
    processed = tmp_path / "processed.tsv"
    export_deseq2_evidence(
        results, gff, processed, tmp_path / "map.tsv", tmp_path / "ann.tsv",
        sample_comparison="test_vs_ctrl", quality_flags=["raw_reanalysis"],
    )
    with processed.open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["standard_error"] == "" and row["p_value"] == "" and row["adjusted_p_value"] == ""
    assert "deseq2_significance_unavailable" in row["quality_flags"]


def test_scripts_do_not_import_sibling_scripts():
    from pathlib import Path

    scripts_dir = Path(__file__).resolve().parents[1] / "scripts"
    offenders = []
    for path in sorted(scripts_dir.glob("*.py")):
        for line in path.read_text().splitlines():
            stripped = line.lstrip()
            if stripped.startswith("from scripts.") or stripped.startswith("import scripts."):
                offenders.append(f"{path.name}: {stripped}")
    assert not offenders
