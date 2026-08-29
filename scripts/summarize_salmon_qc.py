#!/usr/bin/env python3
"""Create a sample-level QC table from Salmon metadata."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


FIELDS = [
    "sample",
    "condition",
    "replicate",
    "run_accession",
    "salmon_version",
    "library_type",
    "processed_fragments",
    "mapped_fragments",
    "mapping_rate_percent",
    "decoy_fragments",
    "decoy_rate_percent",
    "fragment_length_mean",
    "fragment_length_sd",
    "sequence_bias_corrected",
    "gc_bias_corrected",
]


def summarize_salmon_qc(
    quant_dir: Path,
    design_sheet: Path,
    output_path: Path,
    *,
    allow_incomplete: bool = False,
) -> Path:
    """Summarize Salmon ``meta_info.json`` files in design-sheet order."""
    with design_sheet.open(newline="") as handle:
        design_rows = list(csv.DictReader(handle))

    rows = []
    missing = []
    for design in design_rows:
        sample = design["sample"]
        metadata_path = quant_dir / sample / "aux_info" / "meta_info.json"
        if not metadata_path.is_file():
            missing.append(sample)
            continue
        metadata = json.loads(metadata_path.read_text())
        processed = int(metadata["num_processed"])
        decoy_fragments = int(metadata.get("num_decoy_fragments", 0))
        library_types = metadata.get("library_types", [])
        rows.append(
            {
                "sample": sample,
                "condition": design["condition"],
                "replicate": design["replicate"],
                "run_accession": design["run_accession"],
                "salmon_version": metadata["salmon_version"],
                "library_type": ";".join(library_types),
                "processed_fragments": processed,
                "mapped_fragments": int(metadata["num_mapped"]),
                "mapping_rate_percent": f'{float(metadata["percent_mapped"]):.6f}',
                "decoy_fragments": decoy_fragments,
                "decoy_rate_percent": f"{100 * decoy_fragments / processed:.6f}",
                "fragment_length_mean": metadata.get("frag_length_mean", ""),
                "fragment_length_sd": metadata.get("frag_length_sd", ""),
                "sequence_bias_corrected": metadata.get("seq_bias_correct", ""),
                "gc_bias_corrected": metadata.get("gc_bias_correct", ""),
            }
        )

    if missing and not allow_incomplete:
        raise FileNotFoundError(
            "missing Salmon metadata for expected samples: " + ", ".join(missing)
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quant-dir", type=Path, required=True)
    parser.add_argument("--design-sheet", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()
    output = summarize_salmon_qc(
        args.quant_dir,
        args.design_sheet,
        args.output,
        allow_incomplete=args.allow_incomplete,
    )
    print(f"wrote Salmon QC summary: {output}")


if __name__ == "__main__":
    main()
