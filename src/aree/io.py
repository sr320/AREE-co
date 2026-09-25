import pandas as pd


def read_tsv(path, **kwargs):
    # round_trip parsing reads floats back exactly, so rewritten tables do not drift
    # (e.g. 0.018 -> 0.018000000000000002) and match across pandas versions.
    return pd.read_csv(path, sep="\t", float_precision="round_trip", **kwargs)
