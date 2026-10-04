import csv
import gzip
from collections import Counter

from aree.raw.fastq import run_mates
from scripts.prepare_prjna913164_evidence import CONTRASTS, prepare
from scripts.prepare_prjna913164_fastqs import EXCLUDED_SAMPLES, load_manifest


MANIFEST = "data/manifests/CGIG_HEATWAVE_RNASEQ_PRJNA913164_runs.tsv"


def test_manifest_is_balanced_over_treatment_and_ploidy():
    rows = load_manifest(MANIFEST)
    cells = Counter((row["condition"], row["ploidy"]) for row in rows)
    assert set(cells.values()) == {12}
    assert {condition for condition, _ in cells} == {"control", "heat", "heat_emersion"}
    assert {ploidy for _, ploidy in cells} == {"diploid", "triploid"}


def test_every_library_is_single_end_tag_seq_with_a_cited_label():
    for row in load_manifest(MANIFEST):
        assert run_mates(row) == ("1",)
        assert row["read_length_bp"] == "101"
        assert len(row["fastq_1_md5"]) == 32
        assert "attributes.tsv description=" in row["assignment_evidence"]
        assert "'{}, ".format(row["ploidy"]) in row["assignment_evidence"]


def test_qc_exclusions_keep_every_cell_replicated():
    rows = load_manifest(MANIFEST)
    assert EXCLUDED_SAMPLES <= {row["replicate"] for row in rows}
    # The authors' outlier D54 passes AREE's QC and is kept.
    assert "D54" not in EXCLUDED_SAMPLES
    kept = Counter((row["condition"], row["ploidy"]) for row in rows if row["replicate"] not in EXCLUDED_SAMPLES)
    assert min(kept.values()) >= 9


def test_each_stressor_is_its_own_record_sharing_the_controls(tmp_path):
    assert {spec["study_id"] for spec in CONTRASTS.values()} == {
        "CGIG_HEAT_RNASEQ_PRJNA913164", "CGIG_MULTISTRESS_RNASEQ_PRJNA913164",
    }
    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t8.5\t0.5\t0.0001\t0.001\n"
    )
    gff = tmp_path / "reference.gff.gz"
    with gzip.open(gff, "wt") as handle:
        handle.write(
            "chr1\tRefSeq\tgene\t1\t10\t.\t+\t.\t"
            "ID=gene-A;Dbxref=GeneID:1;Name=A;description=alpha;gene=A;gene_biotype=protein_coding\n"
        )
    for contrast, spec in CONTRASTS.items():
        processed = tmp_path / "{}.tsv".format(contrast)
        prepare(contrast, results, gff, processed, tmp_path / "mapping.tsv", tmp_path / "annotations.tsv")
        rows = list(csv.DictReader(processed.open(), delimiter="\t"))
        assert rows[0]["sample_comparison"] == spec["sample_comparison"]
        assert rows[0]["molecular_direction"] == "up"
        flags = rows[0]["quality_flags"].split(";")
        assert {"shared_control_group", "quantseq_3prime_tag_seq", "ploidy_blocked_contrast",
                "dependence_group=PRJNA913164_20C_controls"} <= set(flags)
