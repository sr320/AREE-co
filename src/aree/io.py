from pathlib import Path

import pandas as pd


def read_tsv(path, **kwargs):
    # round_trip parsing reads floats back exactly, so rewritten tables do not drift
    # (e.g. 0.018 -> 0.018000000000000002) and match across pandas versions.
    return pd.read_csv(path, sep="\t", float_precision="round_trip", **kwargs)


def evidence_files(path):
    """The TSVs behind an evidence path: the file itself, or every per-study ``*.tsv`` in a directory."""
    path = Path(path)
    if not path.is_dir():
        return [path]
    files = sorted(path.glob("*.tsv"))
    if not files:
        raise ValueError("No evidence TSVs found in {}".format(path))
    return files


def read_evidence(path, **kwargs):
    """Read harmonized evidence from one TSV or from a directory of per-study TSVs (concatenated)."""
    tables = [read_tsv(file, **kwargs) for file in evidence_files(path)]
    return tables[0] if len(tables) == 1 else pd.concat(tables, ignore_index=True)
