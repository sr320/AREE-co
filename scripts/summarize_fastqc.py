"""Summarize FastQC archives into a compact study-level TSV."""

import argparse
import csv
import zipfile
from pathlib import Path


def parse_fastqc_archive(path):
    with zipfile.ZipFile(path) as archive:
        member = next(name for name in archive.namelist() if name.endswith("/fastqc_data.txt"))
        lines = archive.read(member).decode("utf-8").splitlines()
    metrics = {}
    modules = {}
    current_module = None
    for line in lines:
        if line.startswith(">>END_MODULE"):
            current_module = None
        elif line.startswith(">>"):
            module, status = line[2:].split("\t", 1)
            current_module = module
            modules[module] = status
        elif current_module == "Basic Statistics" and line and not line.startswith("#"):
            key, value = line.split("\t", 1)
            metrics[key] = value
    return metrics, modules


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastqc_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    archives = sorted(args.fastqc_dir.glob("*_fastqc.zip"))
    if not archives:
        raise FileNotFoundError("No FastQC ZIP archives found in {}".format(args.fastqc_dir))

    module_names = set()
    parsed = []
    for archive in archives:
        metrics, modules = parse_fastqc_archive(archive)
        module_names.update(modules)
        parsed.append((archive, metrics, modules))
    ordered_modules = sorted(module_names)
    fieldnames = [
        "sample_file",
        "total_sequences",
        "poor_quality_sequences",
        "sequence_length",
        "gc_percent",
    ] + ["status_{}".format(name.lower().replace(" ", "_")) for name in ordered_modules]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for archive, metrics, modules in parsed:
            row = {
                "sample_file": metrics.get("Filename", archive.stem),
                "total_sequences": metrics.get("Total Sequences", ""),
                "poor_quality_sequences": metrics.get("Sequences flagged as poor quality", ""),
                "sequence_length": metrics.get("Sequence length", ""),
                "gc_percent": metrics.get("%GC", ""),
            }
            for module in ordered_modules:
                row["status_{}".format(module.lower().replace(" ", "_"))] = modules.get(module, "")
            writer.writerow(row)
    print("wrote {} FastQC summaries to {}".format(len(parsed), args.output))


if __name__ == "__main__":
    main()
