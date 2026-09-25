import csv
import gzip

import pytest

from scripts.prepare_prjna694496_evidence import prepare
from scripts.summarize_refseq_annotations import (
    benjamini_hochberg,
    hypergeometric_upper_tail,
)


def test_gene_and_pseudogene_annotations_export_as_exact_mappings(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t1.5\t0.2\t0.001\t0.01\n"
        "NCBI:GeneID:2\t-0.5\t0.3\t0.2\t\n"
    )
    gff = tmp_path / "reference.gff.gz"
    with gzip.open(gff, "wt") as handle:
        handle.write(
            "chr1\tRefSeq\tgene\t1\t10\t.\t+\t.\t"
            "ID=gene-A;Dbxref=GeneID:1;Name=A;description=alpha;gene=A;gene_biotype=protein_coding\n"
            "chr1\tRefSeq\tpseudogene\t20\t30\t.\t+\t.\t"
            "ID=gene-B;Dbxref=GeneID:2;Name=B;description=beta;gene=B;gene_biotype=transcribed_pseudogene\n"
        )
    processed = tmp_path / "processed.tsv"
    mapping = tmp_path / "mapping.tsv"
    annotations = tmp_path / "annotations.tsv"
    prepare(results, gff, processed, mapping, annotations)

    processed_rows = list(csv.DictReader(processed.open(), delimiter="\t"))
    mapping_rows = list(csv.DictReader(mapping.open(), delimiter="\t"))
    annotation_rows = list(csv.DictReader(annotations.open(), delimiter="\t"))
    assert [row["molecular_direction"] for row in processed_rows] == ["up", "down"]
    assert all(row["mapping_confidence"] == "exact" for row in mapping_rows)
    assert annotation_rows[1]["gene_biotype"] == "transcribed_pseudogene"
    assert "deseq2_significance_unavailable" in processed_rows[1]["quality_flags"]


def test_hypergeometric_tail_and_bh_adjustment():
    assert hypergeometric_upper_tail(10, 5, 5, 5) == pytest.approx(1 / 252)
    adjusted = benjamini_hochberg([0.01, 0.04, 0.03])
    assert adjusted == pytest.approx([0.03, 0.04, 0.04])
