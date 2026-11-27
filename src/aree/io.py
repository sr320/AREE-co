from pathlib import Path

import pandas as pd


_PANDAS_MAJOR = int(pd.__version__.split(".")[0])


def read_tsv(path, **kwargs):
    # Pandas < 3.0 needs float_precision="round_trip" so rewritten tables do not drift
    # (e.g. 0.018 -> 0.018000000000000002). Pandas 3.0+ uses round-trip parsing by default
    # and deprecates the argument.
    if _PANDAS_MAJOR < 3 and "float_precision" not in kwargs:
        kwargs["float_precision"] = "round_trip"
    return pd.read_csv(path, sep="\t", **kwargs)


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
