import csv
import gzip

import pandas as pd

from scripts.compare_prjna516762_prjna694496 import compare
from scripts.prepare_prjna516762_evidence import prepare
from scripts.prepare_prjna516762_fastqs import load_manifest


def test_run_manifest_has_three_controls_and_three_heat_libraries():
    rows = load_manifest("data/manifests/CGIG_HEAT_RNASEQ_PRJNA516762_runs.tsv")
    assert {row["condition"] for row in rows} == {"control", "heat"}
    assert {row["library_name"] for row in rows} == {"T1", "T2", "T3", "T7", "T8", "T9"}
    assert all(row["layout"] == "paired" for row in rows)


def test_heat_reanalysis_exports_exact_current_reference_evidence(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t2.0\t0.25\t0.0001\t0.001\n"
        "NCBI:GeneID:2\t-1.0\t0.4\t0.2\t\n"
    )
    gff = tmp_path / "reference.gff.gz"
    with gzip.open(gff, "wt") as handle:
        handle.write(
            "chr1\tRefSeq\tgene\t1\t10\t.\t+\t.\t"
            "ID=gene-A;Dbxref=GeneID:1;Name=A;description=alpha;gene=A;gene_biotype=protein_coding\n"
            "chr1\tRefSeq\tgene\t20\t30\t.\t+\t.\t"
            "ID=gene-B;Dbxref=GeneID:2;Name=B;description=beta;gene=B;gene_biotype=protein_coding\n"
        )
    processed = tmp_path / "processed.tsv"
    mapping = tmp_path / "mapping.tsv"
    annotations = tmp_path / "annotations.tsv"
    prepare(results, gff, processed, mapping, annotations)
    processed_rows = list(csv.DictReader(processed.open(), delimiter="\t"))
    mapping_rows = list(csv.DictReader(mapping.open(), delimiter="\t"))
    assert processed_rows[0]["sample_comparison"] == "heat_35C_2h_vs_control_12C"
    assert processed_rows[0]["molecular_direction"] == "up"
    assert "acute_heat_response" in processed_rows[0]["quality_flags"]
    assert "deseq2_significance_unavailable" in processed_rows[1]["quality_flags"]
    assert all(row["mapping_confidence"] == "exact" for row in mapping_rows)


def test_cross_context_comparison_separates_concordant_and_opposite_results(tmp_path):
    columns = ["feature_id_standardized", "gene_symbol", "log2FoldChange", "lfcSE",
               "stat", "pvalue", "padj"]
    heat = pd.DataFrame([
        ["NCBI:GeneID:1", "A", 2.0, 0.2, 10.0, 1e-8, 1e-7],
        ["NCBI:GeneID:2", "B", -1.0, 0.2, -5.0, 1e-5, 1e-4],
    ], columns=columns)
    selection = pd.DataFrame([
        ["NCBI:GeneID:1", "A", 1.0, 0.2, 5.0, 1e-5, 1e-4],
        ["NCBI:GeneID:2", "B", 1.0, 0.2, 5.0, 1e-5, 1e-4],
    ], columns=columns)
    heat_path, selection_path = tmp_path / "heat.tsv", tmp_path / "selection.tsv"
    heat.to_csv(heat_path, sep="\t", index=False)
    selection.to_csv(selection_path, sep="\t", index=False)
    result = compare(heat_path, selection_path, tmp_path / "out")
    classes = dict(zip(result.feature_id_standardized, result.cross_context_class))
    assert classes["NCBI:GeneID:1"] == "significant_both_concordant"
    assert classes["NCBI:GeneID:2"] == "significant_both_opposite"
