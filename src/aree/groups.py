import math


def is_missing(value):
    return value is None or (isinstance(value, float) and math.isnan(value))


def present(values):
    """Non-missing values, matching pandas' skipna/dropna behavior."""
    return [value for value in values if not is_missing(value)]


def column_groups(frame, by, columns):
    """Yield ``(key, {column: array})`` per group, in sorted key order like ``groupby`` iteration.

    Rows keep their original order within a group and groups with a missing key are skipped, as
    in ``groupby``. Working on plain arrays avoids the per-group pandas overhead that dominated
    runtime on real studies (tens of thousands of groups with a handful of rows each).
    """
    arrays = {column: frame[column].to_numpy() for column in columns}
    groups = frame.groupby(by).indices
    for key in sorted(groups, key=lambda k: k if isinstance(k, tuple) else (k,)):
        if any(is_missing(part) for part in (key if isinstance(key, tuple) else (key,))):
            continue
        positions = groups[key]
        yield key, {column: values[positions] for column, values in arrays.items()}
