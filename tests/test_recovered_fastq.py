import gzip

import pytest

from scripts.validate_recovered_fastq import validate


def inputs(tmp_path):
    paths = [tmp_path / "recovered.gz", tmp_path / "extract1",
             tmp_path / "extract2", tmp_path / "ena2.gz"]
    for index, path in enumerate(paths):
        description = b" length=4" if index != 3 else b" 1/2"
        data = b"@SRR1.1" + description + b"\nACGT\n+\nFFFF\n"
        if path.suffix == ".gz":
            with gzip.open(path, "wb") as handle:
                handle.write(data)
        else:
            path.write_bytes(data)
    return paths


def test_description_differences_are_allowed(tmp_path):
    result = validate(*inputs(tmp_path), "SRR1", 1)
    assert result["validated_pairs"] == 1
    assert result["ena_compressed_md5_equivalent"] is False


@pytest.mark.parametrize("replacement", [
    b"@SRR1.2\nACGT\n+\nFFFF\n",  # wrong spot
    b"@SRR1.1\nTCGT\n+\nFFFF\n",  # changed sequence
    b"@SRR1.1\nACGT\n+\nFFFE\n",  # changed quality
    b"@SRR1.1\nACGT\n+\nFFF\n",   # truncated quality
    b"",                              # missing mate
])
def test_invalid_recovery_rejected(tmp_path, replacement):
    paths = inputs(tmp_path)
    paths[1].write_bytes(replacement)
    with pytest.raises(ValueError):
        validate(*paths, "SRR1", 1)


def test_manifest_count_enforced(tmp_path):
    with pytest.raises(ValueError, match="Expected 2 pairs"):
        validate(*inputs(tmp_path), "SRR1", 2)
