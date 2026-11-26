import gzip
import io

import numpy as np
import pytest

from scripts.check_replicate_genotypes import compare, pileup


def test_pileup_handles_terminal_nh_tag_and_lowercase_reads(tmp_path):
    targets = tmp_path / "targets.txt"
    targets.write_text("tx1\n")
    transcripts = tmp_path / "transcripts.fa.gz"
    with gzip.open(transcripts, "wt") as handle:
        handle.write(">tx1\nACGTACGTACGT\n")

    sam = io.StringIO(
        "@HD\tVN:1.6\n"
        # NH:i:1 as the last tag (no trailing tab) -> kept
        "r1\t0\ttx1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tAS:i:0\tNH:i:1\n"
        # Lowercase bases with terminal NH:i:1 -> kept and mapped to A/C/G/T
        "r2\t0\ttx1\t1\t60\t4M\t*\t0\t0\tacgt\tIIII\tNH:i:1\n"
        # NH:i:10 (multi-mapper) -> skipped
        "r3\t0\ttx1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tNH:i:10\n"
        # Missing NH tag -> skipped
        "r4\t0\ttx1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tAS:i:0\n"
        # Fewer than 12 fields (no optional tags) -> skipped
        "r5\t0\ttx1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
    )
    output = tmp_path / "sample.npz"
    pileup(sam, targets, transcripts, output)

    counts = np.load(output)["tx1"]
    # Positions 0..3 should each have count 2 (r1 + r2) on A, C, G, T respectively.
    assert counts[:4].sum(axis=1).tolist() == [2, 2, 2, 2]
    assert np.diag(counts[:4]).tolist() == [2, 2, 2, 2]


def test_compare_raises_clear_error_when_no_sites_pass_filters(tmp_path):
    paths = []
    for name in ("g_1.npz", "g_2.npz"):
        path = tmp_path / name
        counts = np.zeros((10, 4), dtype=np.uint32)
        counts[:, 0] = 5  # below min_depth=30 and monomorphic
        np.savez_compressed(path, tx1=counts)
        paths.append(str(path))

    with pytest.raises(ValueError, match="No sites met min_depth=30"):
        compare(paths)
