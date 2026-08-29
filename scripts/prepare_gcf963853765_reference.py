"""Download and validate the current Pacific oyster RefSeq reference."""

import argparse
import csv
import gzip
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import unquote

try:
    from scripts.prepare_prjna694496_fastqs import download, nearest_existing_ancestor
except ModuleNotFoundError:  # Direct execution places scripts/ on sys.path.
    from prepare_prjna694496_fastqs import download, nearest_existing_ancestor


ACCESSION = "GCF_963853765.1"
ASSEMBLY = "GCF_963853765.1_xbMagGiga1.1"
BASE_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/853/765/{}/".format(ASSEMBLY)
DEFAULT_OUTPUT = Path("data/reference") / ACCESSION
ARTIFACTS = {
    "{}_genomic.fna.gz".format(ASSEMBLY): (174082374, "e1d6d4f1a8f9c5cd87bb006dd5fd9f98"),
    "{}_genomic.gff.gz".format(ASSEMBLY): (15657002, "f358436291fd77299046268947f8d6bb"),
    "{}_rna.fna.gz".format(ASSEMBLY): (36305575, "0321e363d5f15f9ac7d13e69ff20963a"),
}


def parse_attributes(text):
    attributes = {}
    for item in text.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = unquote(value)
    return attributes


def build_tx2gene(gff_path, output_path):
    mappings = {}
    with gzip.open(gff_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attributes = parse_attributes(fields[8])
            transcript_id = attributes.get("transcript_id")
            dbxref = attributes.get("Dbxref", "")
            gene_ids = [part.split(":", 1)[1] for part in dbxref.split(",") if part.startswith("GeneID:")]
            if not transcript_id or len(gene_ids) != 1:
                continue
            mapping = ("NCBI:GeneID:{}".format(gene_ids[0]), attributes.get("gene", ""))
            if transcript_id in mappings and mappings[transcript_id] != mapping:
                raise ValueError("Conflicting gene assignments for {}".format(transcript_id))
            mappings[transcript_id] = mapping
    if not mappings:
        raise ValueError("No transcript-to-gene mappings found in {}".format(gff_path))
    with Path(output_path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["transcript_id", "gene_id", "gene_symbol"])
        for transcript_id, (gene_id, gene_symbol) in sorted(mappings.items()):
            writer.writerow([transcript_id, gene_id, gene_symbol])
    print("wrote {} transcript-to-gene mappings to {}".format(len(mappings), output_path))


def build_gentrome(transcriptome_path, genome_path, gentrome_path, decoys_path):
    with Path(gentrome_path).open("wb") as target:
        for source_path in (transcriptome_path, genome_path):
            with gzip.open(source_path, "rb") as source:
                shutil.copyfileobj(source, target, length=1024 * 1024)
    decoys = []
    with gzip.open(genome_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                decoys.append(line[1:].split(None, 1)[0])
    if not decoys:
        raise ValueError("No genome decoy sequences found in {}".format(genome_path))
    with Path(decoys_path).open("w", encoding="utf-8") as handle:
        handle.write("\n".join(decoys) + "\n")
    print("wrote decoy-aware gentrome {} with {} genome decoys".format(gentrome_path, len(decoys)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--retries", type=int, default=12)
    parser.add_argument("--build-derived-only", action="store_true")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    if args.retries < 0:
        parser.error("--retries cannot be negative")

    required_bytes = sum(size for size, _ in ARTIFACTS.values()) * 3
    probe_dir = nearest_existing_ancestor(args.output_dir)
    free_bytes = shutil.disk_usage(str(probe_dir)).free
    print("reference capacity checked at: {}".format(probe_dir))
    print("required free space: {:.2f} GB".format(required_bytes / 1e9))
    print("available free space: {:.2f} GB".format(free_bytes / 1e9))
    if free_bytes < required_bytes:
        raise OSError("Insufficient reference scratch space")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    if not args.build_derived_only:
        jobs = []
        for filename, (expected_bytes, expected_md5) in ARTIFACTS.items():
            jobs.append(
                (
                    BASE_URL + filename,
                    args.output_dir / filename,
                    expected_bytes,
                    expected_md5,
                    args.retries,
                )
            )
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(download, *job) for job in jobs]
            for future in as_completed(futures):
                future.result()

    gff_path = args.output_dir / "{}_genomic.gff.gz".format(ASSEMBLY)
    genome_path = args.output_dir / "{}_genomic.fna.gz".format(ASSEMBLY)
    transcriptome_path = args.output_dir / "{}_rna.fna.gz".format(ASSEMBLY)
    for required_path in (gff_path, genome_path, transcriptome_path):
        if not required_path.exists():
            raise FileNotFoundError("Missing reference artifact: {}".format(required_path))
    build_tx2gene(gff_path, args.output_dir / "{}_tx2gene.tsv".format(ACCESSION))
    build_gentrome(
        transcriptome_path,
        genome_path,
        args.output_dir / "{}_gentrome.fa".format(ACCESSION),
        args.output_dir / "{}_decoys.txt".format(ACCESSION),
    )


if __name__ == "__main__":
    main()
