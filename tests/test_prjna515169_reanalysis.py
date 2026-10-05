import csv
import gzip

from aree.raw.fastq import run_mates
from scripts.prepare_prjna515169_evidence import CONTRASTS, prepare
from scripts.prepare_prjna515169_fastqs import load_manifest


MANIFEST = "data/manifests/CGIG_VIBRIO_RNASEQ_PRJNA515169_runs.tsv"
ALIAS = {"control": "RNAseq_anesthesis", "vcrassostreae_J2_9": "RNAseq_Vcrass_J2-9",
         "vtasmaniensis_LGP32": "RNAseq_Vtasma_LGP32"}


def test_manifest_has_one_pool_per_condition_per_experiment():
    rows = load_manifest(MANIFEST)
    for condition in ALIAS:
        assert sorted(row["experiment"] for row in rows if row["condition"] == condition) == ["exp1", "exp2", "exp3"]
    for row in rows:
        assert "alias={}_replicate{} ".format(ALIAS[row["condition"]], row["replicate"]) in row["assignment_evidence"]
        assert row["experiment"] == "exp" + row["replicate"]
        assert run_mates(row) == ("1", "2")
        assert len(row["fastq_1_md5"]) == len(row["fastq_2_md5"]) == 32


def test_both_records_share_the_mock_pool_dependence_group(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t3.0\t0.5\t0.0001\t0.001\n"
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
        flags = rows[0]["quality_flags"].split(";")
        assert {"dependence_group=PRJNA515169_mock_pools", "pooled_libraries"} <= set(flags)
