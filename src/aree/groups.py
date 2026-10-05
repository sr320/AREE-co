import math


def dependence_group(study_id, flags):
    """Explicit replication groups travel with evidence outside the checkout."""
    groups = {flag.split("=", 1)[1] for flag in str(flags).split(";")
              if flag.startswith("dependence_group=")}
    if len(groups) > 1 or "" in groups:
        raise ValueError("Invalid dependence group for {}".format(study_id))
    return "group:" + next(iter(groups)) if groups else "study:" + str(study_id)


def dependence_groups(evidence):
    flags = evidence["quality_flags"] if "quality_flags" in evidence else [None] * len(evidence)
    groups, seen = [], {}
    for study, flag in zip(evidence["study_id"], flags):
        group = dependence_group(study, flag)
        if study in seen and seen[study] != group:
            raise ValueError("Inconsistent dependence groups within a study")
        seen[study] = group
        groups.append(group)
    return groups


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
