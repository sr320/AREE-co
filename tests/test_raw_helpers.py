import gzip
import hashlib
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

import pytest

from aree.raw.fastq import download, download_once, load_run_manifest
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
