import csv
import gzip

from aree.raw.fastq import run_mates
from scripts.prepare_prjna593309_evidence import CONTRASTS as EVIDENCE
from scripts.prepare_prjna593309_evidence import prepare
from scripts.prepare_prjna593309_fastqs import CONTRASTS, contrast_rows, load_manifest


MANIFEST = "data/manifests/CGIG_OSHV1_TEMP_RNASEQ_PRJNA593309_runs.tsv"
ALIAS = {"ctr_21C_0h": "Cg_CTR_21_0hpc", "oshv_21C_48h": "Cg_OsHV_21_48hpc",
         "oshv_21C_12h": "Cg_OsHV_21_12hpc", "oshv_29C_12h": "Cg_OsHV_29_12hpc"}


def test_manifest_groups_agree_with_the_sample_aliases():
    for row in load_manifest(MANIFEST):
        assert "alias={}_{};".format(ALIAS[row["condition"]], row["replicate"]) in row["assignment_evidence"]
        assert run_mates(row) == ("1", "2")
        assert len(row["fastq_1_md5"]) == len(row["fastq_2_md5"]) == 32


def test_each_record_is_three_tanks_per_arm():
    rows = load_manifest(MANIFEST)
    for contrast, labels in CONTRASTS.items():
        design = contrast_rows(rows, contrast)
        assert sorted(row["condition"] for row in design) == sorted([*labels] * 3)
        assert len({row["tank"] for row in design}) == (3 if contrast == "oshv1" else 6)
    # The OsHV-1 record's 0 h and 48 h oysters share the three 21 C tanks, so tank can pair them.
    oshv = contrast_rows(rows, "oshv1")
    assert {row["tank"] for row in oshv if row["condition"] == "control"} == \
        {row["tank"] for row in oshv if row["condition"] == "oshv1"}


def test_records_export_with_the_shared_tank_dependence_group(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t2.0\t0.4\t0.0001\t0.001\n"
    )
    gff = tmp_path / "reference.gff.gz"
    with gzip.open(gff, "wt") as handle:
        handle.write(
            "chr1\tRefSeq\tgene\t1\t10\t.\t+\t.\t"
            "ID=gene-A;Dbxref=GeneID:1;Name=A;description=alpha;gene=A;gene_biotype=protein_coding\n"
        )
    for contrast, spec in EVIDENCE.items():
        processed = tmp_path / "{}.tsv".format(contrast)
        prepare(contrast, results, gff, processed, tmp_path / "mapping.tsv", tmp_path / "annotations.tsv")
        rows = list(csv.DictReader(processed.open(), delimiter="\t"))
        assert rows[0]["sample_comparison"] == spec["sample_comparison"]
        assert "dependence_group=PRJNA593309_21C_tanks" in rows[0]["quality_flags"].split(";")


def test_the_duplicated_29c_oyster_is_summed_into_one_animal():
    from scripts.prepare_prjna593309_fastqs import SAME_ANIMAL

    design = contrast_rows(load_manifest(MANIFEST), "heat")
    animals = {row["animal"] for row in design if row["condition"] == "heat"}
    assert SAME_ANIMAL == {"oshv_29C_12h_2": "oshv_29C_12h_1"}
    assert len(animals) == 2
    assert len({row["animal"] for row in design if row["condition"] == "control"}) == 3
