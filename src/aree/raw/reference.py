"""Versioned multi-species RefSeq catalog, verified fetching, tx2gene and gentrome builders."""

import gzip
import hashlib
import http.client
import json
import os
import re
import shutil
import socket
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml

from aree.raw.fastq import download as fastq_download
from aree.raw.fastq import md5sum
from aree.raw.gff import parse_attributes


DEFAULT_CATALOG_PATH = Path(__file__).resolve().parents[3] / "config" / "references.yaml"
NETWORK_ERRORS = (socket.timeout, TimeoutError, ConnectionError, urllib.error.URLError, http.client.IncompleteRead)
_GTF_ATTR_RE = re.compile(r'([A-Za-z0-9_]+)\s+"([^"]*)"')
_GENIID_RE = re.compile(r"GeneID:(\d+)")


def _open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return Path(path).open("r", encoding="utf-8", errors="replace")


def _sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _verify_file(
    path: Path,
    expected_bytes: Optional[int] = None,
    expected_md5: Optional[str] = None,
    expected_sha256: Optional[str] = None,
) -> bool:
    if expected_bytes is not None and path.stat().st_size != int(expected_bytes):
        return False
    if expected_md5 and md5sum(path) != expected_md5:
        return False
    if expected_sha256 and _sha256sum(path) != expected_sha256:
        return False
    return True


def load_reference_catalog(catalog_path: Optional[Path] = None) -> Dict[str, dict]:
    """Load the YAML reference catalog mapping RefSeq assembly accessions to metadata and files."""
    path = Path(catalog_path) if catalog_path else DEFAULT_CATALOG_PATH
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict) or "references" not in data:
        raise ValueError("Invalid reference catalog at {}".format(path))
    return data


def _normalize_key(text: str) -> str:
    return re.sub(r"[\s_]+", " ", text.strip().lower())


def list_references(catalog_path: Optional[Path] = None) -> List[dict]:
    """Return all configured RefSeq entries sorted by assembly accession."""
    catalog = load_reference_catalog(catalog_path)
    return [catalog["references"][acc] for acc in sorted(catalog["references"])]


def get_reference_spec(identifier: Optional[str] = None, catalog_path: Optional[Path] = None) -> dict:
    """Resolve an assembly accession, species name, synonym, or species_code to a catalog entry."""
    catalog = load_reference_catalog(catalog_path)
    refs = catalog["references"]
    query = identifier or catalog.get("default_reference")
    if not query:
        raise ValueError("No reference identifier provided and no default_reference configured")
    if query in refs:
        return dict(refs[query])
    norm = _normalize_key(str(query))
    for acc, spec in refs.items():
        candidates = [acc, spec.get("species", ""), spec.get("species_code", "")] + list(spec.get("synonyms", []))
        if any(_normalize_key(c) == norm for c in candidates if c):
            return dict(spec)
    available = ", ".join("{} ({})".format(k, v.get("species", "")) for k, v in sorted(refs.items()))
    raise KeyError("Unknown reference identifier {!r}. Available: {}".format(identifier, available))


def _download_stream(
    url: str,
    destination: Path,
    expected_bytes: Optional[int] = None,
    expected_md5: Optional[str] = None,
    expected_sha256: Optional[str] = None,
    retries: int = 6,
    timeout: int = 60,
) -> Path:
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if _verify_file(destination, expected_bytes, expected_md5, expected_sha256):
            return destination
        raise ValueError("Existing file fails verification: {}".format(destination))

    partial = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(retries + 1):
        try:
            offset = partial.stat().st_size if partial.exists() else 0
            request = urllib.request.Request(url)
            if offset and not str(url).startswith("file:"):
                request.add_header("Range", "bytes={}-".format(offset))
            with urllib.request.urlopen(request, timeout=timeout) as source:
                append = offset > 0 and getattr(source, "status", None) == 206
                with partial.open("ab" if append else "wb") as target:
                    shutil.copyfileobj(source, target, length=1024 * 1024)
            if not _verify_file(partial, expected_bytes, expected_md5, expected_sha256):
                raise ValueError("Downloaded file fails checksum/size verification: {}".format(partial))
            os.replace(str(partial), str(destination))
            return destination
        except NETWORK_ERRORS:
            if attempt == retries:
                raise
            time.sleep(min(2 ** attempt, 30))
    return destination


def download_artifact(
    artifact: dict,
    destination: Path,
    mirror_dir: Optional[Path] = None,
    retries: int = 6,
    timeout: int = 60,
) -> Path:
    """Fetch and verify a single catalog artifact, checking optional local mirror first."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    expected_bytes = artifact.get("bytes")
    expected_md5 = artifact.get("md5")
    expected_sha256 = artifact.get("sha256")

    if mirror_dir is not None:
        mirror_dir = Path(mirror_dir)
        gz_candidate = mirror_dir / artifact.get("filename", destination.name)
        short_candidate = mirror_dir / destination.name
        if gz_candidate.exists():
            if _verify_file(gz_candidate, expected_bytes, expected_md5, expected_sha256):
                if gz_candidate.resolve() != destination.resolve():
                    shutil.copy2(gz_candidate, destination)
                return destination
        elif short_candidate.exists():
            if _verify_file(short_candidate, expected_bytes, expected_md5, expected_sha256):
                if short_candidate.resolve() != destination.resolve():
                    shutil.copy2(short_candidate, destination)
                return destination

    url = artifact["url"]
    if expected_bytes is not None and expected_md5 and not str(url).startswith("file:"):
        fastq_download(url, destination, int(expected_bytes), expected_md5, retries=retries, timeout=timeout)
        return destination
    return _download_stream(
        url,
        destination,
        expected_bytes=expected_bytes,
        expected_md5=expected_md5,
        expected_sha256=expected_sha256,
        retries=retries,
        timeout=timeout,
    )


def _extract_gtf_attrs(text: str) -> Tuple[Dict[str, str], Optional[str]]:
    attrs = {}  # type: Dict[str, str]
    gene_id_num = None
    for key, val in _GTF_ATTR_RE.findall(text):
        if key == "db_xref" and val.startswith("GeneID:"):
            gene_id_num = val.split(":", 1)[1]
        elif key not in attrs:
            attrs[key] = val
    if gene_id_num is None:
        match = _GENIID_RE.search(text)
        if match:
            gene_id_num = match.group(1)
    return attrs, gene_id_num


def build_tx2gene_from_gff(annotation_path: Path, output_path: Path) -> int:
    """Build sorted ``transcript_id\\tgene_id`` (``NCBI:GeneID:<id>``) from a GFF3 or GTF file."""
    annotation_path = Path(annotation_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    is_gtf_hint = ".gtf" in annotation_path.name.lower()
    gene_by_gff_id = {}  # type: Dict[str, str]
    tx_to_gene = {}  # type: Dict[str, str]

    with _open_text(annotation_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            attr_text = fields[8]
            is_gtf = is_gtf_hint or ('"' in attr_text and "=" not in attr_text.split(";", 1)[0])
            if is_gtf:
                attrs, gene_num = _extract_gtf_attrs(attr_text)
                tx_id = attrs.get("transcript_id", "").strip()
                if tx_id and gene_num:
                    tx_to_gene[tx_id] = "NCBI:GeneID:{}".format(gene_num)
                continue

            attrs = parse_attributes(attr_text)
            dbxref = attrs.get("Dbxref", "")
            gene_match = _GENIID_RE.search(dbxref)
            gene_id = "NCBI:GeneID:{}".format(gene_match.group(1)) if gene_match else None
            if feature == "gene" and attrs.get("ID") and gene_id:
                gene_by_gff_id[attrs["ID"]] = gene_id
                continue
            tx_id = attrs.get("transcript_id", "").strip()
            if not tx_id:
                continue
            resolved_gene = gene_id
            if not resolved_gene:
                for parent in attrs.get("Parent", "").split(","):
                    if parent in gene_by_gff_id:
                        resolved_gene = gene_by_gff_id[parent]
                        break
            if resolved_gene:
                tx_to_gene[tx_id] = resolved_gene

    pairs = sorted(tx_to_gene.items(), key=lambda item: (item[1], item[0]))
    with output_path.open("w", encoding="utf-8", newline="") as out:
        out.write("transcript_id\tgene_id\n")
        for tx_id, gene_id in pairs:
            out.write("{}\t{}\n".format(tx_id, gene_id))
    return len(pairs)


def write_minimal_gff_from_gtf(gtf_path: Path, gff_out_path: Path) -> int:
    """Convert RefSeq GTF gene features to a compact GFF3 suitable for ``aree.raw.gff.gene_annotations``."""
    gtf_path = Path(gtf_path)
    gff_out_path = Path(gff_out_path)
    gff_out_path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    seen_genes = set()
    opener = gzip.open if str(gff_out_path).endswith(".gz") else Path.open
    with _open_text(gtf_path) as src, opener(gff_out_path, "wt", encoding="utf-8") as dst:
        dst.write("##gff-version 3\n")
        for line in src:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs, gene_num = _extract_gtf_attrs(fields[8])
            if not gene_num or gene_num in seen_genes:
                continue
            seen_genes.add(gene_num)
            symbol = attrs.get("gene") or attrs.get("gene_id") or "LOC{}".format(gene_num)
            desc = attrs.get("description", "")
            biotype = attrs.get("gene_biotype", "")
            gff_attrs = [
                "ID=gene-{}".format(urllib.parse.quote(symbol, safe="-_.")),
                "Dbxref=GeneID:{}".format(gene_num),
                "Name={}".format(urllib.parse.quote(symbol, safe="-_.")),
            ]
            if desc:
                gff_attrs.append("description={}".format(urllib.parse.quote(desc, safe="-_. ")))
            if biotype:
                gff_attrs.append("gene_biotype={}".format(urllib.parse.quote(biotype, safe="-_.")))
            row = fields[:8] + [";".join(gff_attrs)]
            dst.write("\t".join(row) + "\n")
            written += 1
    return written


def build_gentrome(
    rna_fna: Path,
    genomic_fna: Path,
    gentrome_out: Path,
    decoys_out: Path,
) -> Tuple[Path, Path]:
    """Concatenate transcriptome FASTA and genome FASTA into gentrome + decoys.txt for Salmon."""
    rna_fna = Path(rna_fna)
    genomic_fna = Path(genomic_fna)
    gentrome_out = Path(gentrome_out)
    decoys_out = Path(decoys_out)
    gentrome_out.parent.mkdir(parents=True, exist_ok=True)
    decoys_out.parent.mkdir(parents=True, exist_ok=True)

    out_opener = gzip.open if str(gentrome_out).endswith(".gz") else Path.open
    with out_opener(gentrome_out, "wt", encoding="utf-8") as gentrome, \
            decoys_out.open("w", encoding="utf-8") as decoys:
        with _open_text(rna_fna) as rna_handle:
            shutil.copyfileobj(rna_handle, gentrome)
        with _open_text(genomic_fna) as genome_handle:
            for line in genome_handle:
                gentrome.write(line)
                if line.startswith(">"):
                    seq_id = line[1:].strip().split()[0]
                    if seq_id:
                        decoys.write("{}\n".format(seq_id))
    return gentrome_out, decoys_out


def fetch_reference(
    identifier: Optional[str] = None,
    output_dir: Optional[Path] = None,
    include_genome: bool = False,
    catalog_path: Optional[Path] = None,
    mirror_dir: Optional[Path] = None,
) -> dict:
    """Download/stage a cataloged RefSeq bundle, build ``tx2gene.tsv`` (and gentrome if requested), and write manifest."""
    spec = get_reference_spec(identifier, catalog_path=catalog_path)
    accession = spec["accession"]
    dest_dir = Path(output_dir) if output_dir else Path("data/reference") / accession
    dest_dir.mkdir(parents=True, exist_ok=True)

    files = spec["files"]
    rna_path = download_artifact(files["rna_fna"], dest_dir / "rna.fna.gz", mirror_dir=mirror_dir)

    annotation_path = None
    gff_path = dest_dir / "genomic.gff.gz"
    gtf_path = None

    # Prefer genomic_gff if explicitly verified or local mirror holds GFF; otherwise use verified genomic_gtf.
    gff_meta = files.get("genomic_gff")
    gtf_meta = files.get("genomic_gtf")

    if gtf_meta and (gtf_meta.get("md5") or gtf_meta.get("sha256")):
        # Check if mirror has plain genomic.gtf (common in prepared /Volumes/ugreen-2TB/aree/reference bundles)
        plain_mirror_gtf = Path(mirror_dir) / "genomic.gtf" if mirror_dir else None
        if plain_mirror_gtf and plain_mirror_gtf.exists():
            gtf_path = plain_mirror_gtf
        else:
            gtf_path = download_artifact(gtf_meta, dest_dir / "genomic.gtf.gz", mirror_dir=mirror_dir)
        annotation_path = gtf_path

    if gff_meta and (gff_meta.get("md5") or (mirror_dir and (Path(mirror_dir) / gff_meta["filename"]).exists()) or annotation_path is None):
        gff_path = download_artifact(gff_meta, dest_dir / "genomic.gff.gz", mirror_dir=mirror_dir)
        annotation_path = gff_path
    elif not gff_path.exists() and gtf_path is not None:
        write_minimal_gff_from_gtf(gtf_path, gff_path)

    tx2gene_path = dest_dir / "tx2gene.tsv"
    n_tx = build_tx2gene_from_gff(annotation_path, tx2gene_path)

    gentrome_path = None
    decoys_path = None
    if include_genome and "genomic_fna" in files:
        genome_path = download_artifact(files["genomic_fna"], dest_dir / "genomic.fna.gz", mirror_dir=mirror_dir)
        gentrome_path, decoys_path = build_gentrome(
            rna_path,
            genome_path,
            dest_dir / "gentrome.fa.gz",
            dest_dir / "decoys.txt",
        )

    manifest = {
        "accession": accession,
        "species": spec["species"],
        "species_code": spec.get("species_code"),
        "ncbi_taxid": spec.get("ncbi_taxid"),
        "assembly_id": spec.get("assembly_id"),
        "annotation_release": spec["annotation_release"],
        "transcripts": str((gentrome_path or rna_path).resolve()),
        "rna_fna": str(rna_path.resolve()),
        "gff": str(gff_path.resolve()),
        "gtf": str(gtf_path.resolve()) if gtf_path else None,
        "tx2gene": str(tx2gene_path.resolve()),
        "decoys": str(decoys_path.resolve()) if decoys_path else None,
        "tx2gene_count": n_tx,
    }
    manifest_path = dest_dir / "reference_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    manifest["manifest_path"] = str(manifest_path.resolve())
    return manifest
