"""Validate clean FASTQ replacements and quarantine failed resumed files."""

import argparse
import csv
import os
from pathlib import Path

from prepare_prjna694496_fastqs import DEFAULT_MANIFEST, md5sum


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--fastq-dir", type=Path, required=True)
    args = parser.parse_args()

    with args.manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    promoted = 0
    for row in rows:
        run = row["run_accession"]
        for mate in ("1", "2"):
            destination = args.fastq_dir / "{}_{}.fastq.gz".format(run, mate)
            replacement = destination.with_suffix(destination.suffix + ".fresh")
            failed_resume = destination.with_suffix(destination.suffix + ".part")
            quarantine = destination.with_suffix(destination.suffix + ".md5-failed")
            if not replacement.exists():
                continue
            expected_bytes = int(row["fastq_{}_bytes".format(mate)])
            expected_md5 = row["fastq_{}_md5".format(mate)]
            if replacement.stat().st_size != expected_bytes:
                raise ValueError("Replacement has wrong byte count: {}".format(replacement))
            actual_md5 = md5sum(replacement)
            if actual_md5 != expected_md5:
                raise ValueError("Replacement fails ENA MD5 validation: {}".format(replacement))
            if destination.exists():
                raise FileExistsError("Refusing to overwrite final FASTQ: {}".format(destination))
            if failed_resume.exists():
                if quarantine.exists():
                    raise FileExistsError("Quarantine path already exists: {}".format(quarantine))
                os.replace(str(failed_resume), str(quarantine))
                print("quarantined failed resume {}".format(quarantine))
            os.replace(str(replacement), str(destination))
            print("promoted verified replacement {}".format(destination))
            promoted += 1
    if promoted == 0:
        raise ValueError("No .fresh replacement files found")
    print("promoted {} checksum-verified replacements".format(promoted))


if __name__ == "__main__":
    main()
