"""Validate an SRA-recovered mate against its extraction and ENA mate 2.

Headers may differ in their descriptions; accession.spot IDs, sequences and
qualities must agree. This is not a claim of ENA compressed-file MD5 identity.
"""

import argparse
import gzip
import hashlib
import io
import json
from contextlib import ExitStack
from pathlib import Path


def record(handle):
    header = handle.readline()
    if not header:
        return None
    sequence, plus, quality = (handle.readline().rstrip(b"\r\n") for _ in range(3))
    if not header.startswith(b"@") or not plus.startswith(b"+") or not sequence or len(sequence) != len(quality):
        raise ValueError("Malformed FASTQ record")
    return header.split()[0], sequence, quality


def validate(recovered, extracted1, extracted2, ena2, run, expected):
    paths = [recovered, extracted1, extracted2, ena2]
    digests = [hashlib.sha256() for _ in paths]
    count = 0
    with ExitStack() as stack:
        # Large buffers avoid alternating tiny reads across four external-drive files.
        handles = [stack.enter_context(io.BufferedReader(gzip.open(p, "rb"), 4 * 1024 * 1024)
                   if ".gz" in p.name else p.open("rb", buffering=4 * 1024 * 1024)) for p in paths]
        while True:
            rows = [record(h) for h in handles]
            if all(r is None for r in rows):
                break
            if any(r is None for r in rows):
                raise ValueError("FASTQs have unequal record counts")
            count += 1
            if len({r[0] for r in rows}) != 1 or not rows[0][0].startswith(("@" + run + ".").encode()):
                raise ValueError(f"Mismatched spot IDs at record {count}")
            if rows[0] != rows[1] or rows[2] != rows[3]:
                raise ValueError(f"Recovered/extracted/ENA content mismatch at record {count}")
            for digest, row in zip(digests, rows):
                digest.update(b"\n".join(row) + b"\n")
            if count % 1000000 == 0:
                print(f"Validated {count:,} / {expected:,} pairs", flush=True)
    if count != expected:
        raise ValueError(f"Expected {expected} pairs, observed {count}")
    return {"run_accession": run, "validated_pairs": count, "status": "pass",
            "validation": "All spot IDs paired; recovered mate 1 matches SRA extraction; extracted mate 2 matches ENA sequences and qualities.",
            "ena_compressed_md5_equivalent": False,
            "normalized_sha256": {str(p): d.hexdigest() for p, d in zip(paths, digests)}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("recovered", "extracted1", "extracted2", "ena2", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--expected", type=int, required=True)
    args = parser.parse_args()
    result = validate(args.recovered, args.extracted1, args.extracted2, args.ena2, args.run, args.expected)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
