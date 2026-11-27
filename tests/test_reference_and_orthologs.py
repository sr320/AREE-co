import gzip
import hashlib
from pathlib import Path

import pandas as pd
import pytest
import yaml

from aree.harmonize.identifiers import load_mapping, map_identifier
from aree.harmonize.orthologs import (
    build_all_multispecies_maps,
    build_mgigas_orthogroup_index,
    build_species_to_mgigas_map,
    load_single_copy_orthogroups,
)
from aree.raw.reference import (
    build_tx2gene_from_gff,
    fetch_reference,
    get_reference_spec,
    list_references,
)

ROOT = Path(__file__).resolve().parents[1]


def test_reference_catalog_resolves_all_eight_bivalves():
    refs = list_references()
    assert len(refs) == 8
    assert get_reference_spec("mgigas")["accession"] == "GCF_963853765.1"
    assert get_reference_spec("Crassostrea gigas")["accession"] == "GCF_963853765.1"
    assert get_reference_spec("crassostrea_virginica")["accession"] == "GCF_053477285.1"
    assert get_reference_spec("GCF_026571515.1")["species_code"] == "rphil"
    with pytest.raises(KeyError):
        get_reference_spec("UNKNOWN_SP")


def test_build_tx2gene_supports_both_gff3_and_gtf(tmp_path):
    gff_path = tmp_path / "sample.gff.gz"
    with gzip.open(gff_path, "wt", encoding="utf-8") as h:
        h.write("##gff-version 3\n")
        h.write("NC_1\tRefSeq\tgene\t1\t500\t.\t+\t.\tID=gene-A;Dbxref=GeneID:105300001;Name=GENE_A\n")
        h.write("NC_1\tRefSeq\tmRNA\t1\t500\t.\t+\t.\tID=rna-1;Parent=gene-A;transcript_id=XM_0001.1\n")
    tx_gff = tmp_path / "tx_gff.tsv"
    assert build_tx2gene_from_gff(gff_path, tx_gff) == 1
    df_gff = pd.read_csv(tx_gff, sep="\t")
    assert df_gff.to_dict(orient="records") == [{"transcript_id": "XM_0001.1", "gene_id": "NCBI:GeneID:105300001"}]

    gtf_path = tmp_path / "sample.gtf.gz"
    with gzip.open(gtf_path, "wt", encoding="utf-8") as h:
        h.write('NC_1\tGnomon\tgene\t1\t500\t.\t+\t.\tgene_id "GENE_B"; db_xref "GeneID:111100002"; description "test gene";\n')
        h.write('NC_1\tGnomon\ttranscript\t1\t500\t.\t+\t.\tgene_id "GENE_B"; transcript_id "XM_0002.1"; db_xref "GeneID:111100002";\n')
    tx_gtf = tmp_path / "tx_gtf.tsv"
    assert build_tx2gene_from_gff(gtf_path, tx_gtf) == 1
    df_gtf = pd.read_csv(tx_gtf, sep="\t")
    assert df_gtf.to_dict(orient="records") == [{"transcript_id": "XM_0002.1", "gene_id": "NCBI:GeneID:111100002"}]


def test_fetch_reference_local_file_urls_with_genome_decoys(tmp_path):
    src = tmp_path / "upstream"
    src.mkdir()
    rna = src / "rna.fna.gz"
    gtf = src / "genomic.gtf.gz"
    fna = src / "genomic.fna.gz"

    with gzip.open(rna, "wt") as h:
        h.write(">XM_1.1\nATGCATGCATGC\n")
    with gzip.open(gtf, "wt") as h:
        h.write('NC_1\tGnomon\tgene\t1\t500\t.\t+\t.\tgene_id "HSP70"; db_xref "GeneID:105312345"; gene "HSP70"; gene_biotype "protein_coding";\n')
        h.write('NC_1\tGnomon\ttranscript\t1\t500\t.\t+\t.\tgene_id "HSP70"; transcript_id "XM_1.1"; db_xref "GeneID:105312345";\n')
    with gzip.open(fna, "wt") as h:
        h.write(">NC_1 chromosome 1\nGGGGCCCCAAAATTTT\n")

    def md5_hex(p):
        return hashlib.md5(p.read_bytes()).hexdigest()

    cat_path = tmp_path / "catalog.yaml"
    cat_path.write_text(yaml.safe_dump({
        "version": 1,
        "default_reference": "GCF_TEST.1",
        "references": {
            "GCF_TEST.1": {
                "accession": "GCF_TEST.1",
                "species": "Crassostrea virginica",
                "species_code": "cvirg",
                "ncbi_taxid": 6565,
                "assembly_id": "TESTv1",
                "annotation_release": "GCF_TEST.1-RS_2026_01",
                "files": {
                    "rna_fna": {"filename": rna.name, "url": rna.resolve().as_uri(), "md5": md5_hex(rna)},
                    "genomic_gtf": {"filename": gtf.name, "url": gtf.resolve().as_uri(), "md5": md5_hex(gtf)},
                    "genomic_fna": {"filename": fna.name, "url": fna.resolve().as_uri(), "md5": md5_hex(fna)},
                },
            }
        },
    }))

    dest = tmp_path / "out_ref"
    manifest = fetch_reference("cvirg", output_dir=dest, include_genome=True, catalog_path=cat_path)
    assert manifest["accession"] == "GCF_TEST.1"
    assert manifest["tx2gene_count"] == 1
    assert Path(manifest["tx2gene"]).read_text().strip().endswith("XM_1.1\tNCBI:GeneID:105312345")
    assert Path(manifest["decoys"]).read_text().strip() == "NC_1"
    assert Path(manifest["gff"]).exists()
    assert Path(manifest["transcripts"]).name == "gentrome.fa.gz"


def test_ortholog_builders_cardinality_and_identifier_lookup(tmp_path):
    of_dir = tmp_path / "Results"
    (of_dir / "Orthogroups").mkdir(parents=True)
    (of_dir / "Orthologues" / "Orthologues_mgigas").mkdir(parents=True)

    (of_dir / "Orthogroups" / "Orthogroups_SingleCopyOrthologues.txt").write_text("OG0001\n")
    (of_dir / "Orthogroups" / "Orthogroups.tsv").write_text(
        "Orthogroup\tcvirg\tmcali\tmedul\tmgall\tmgigas\tmmerc\tmtros\trphil\n"
        "OG0001\t201\t301\t401\t501\t101\t601\t701\t801\n"
        "OG0002\t202, 203\t\t\t\t102\t\t\t\n"
        "OG0003\t204\t\t\t\t103, 104\t\t\t\n"
        "OG0004\t205, 206\t\t\t\t105, 106\t\t\t\n"
    )
    pair_path = of_dir / "Orthologues" / "Orthologues_mgigas" / "mgigas__v__cvirg.tsv"
    pair_path.write_text(
        "Orthogroup\tmgigas\tcvirg\n"
        "OG0001\t101\t201\n"
        "OG0002\t102\t202, 203\n"
        "OG0003\t103, 104\t204\n"
        "OG0004\t105, 106\t205, 206\n"
    )

    sc = load_single_copy_orthogroups(of_dir / "Orthogroups" / "Orthogroups_SingleCopyOrthologues.txt")
    map_out = tmp_path / "cvirg_map.tsv"
    summary = build_species_to_mgigas_map(pair_path, "cvirg", single_copy_ogs=sc, output_path=map_out)
    assert summary["one_to_one"] == 1
    assert summary["many_to_one"] == 2
    assert summary["one_to_many"] == 1
    assert summary["many_to_many_inferred"] == 2

    mapping = load_mapping(map_out)
    assert len(mapping) == 6
    assert map_identifier("NCBI:GeneID:201", mapping)["feature_id_standardized"] == "NCBI:GeneID:101"
    assert map_identifier("NCBI:GeneID:201", mapping)["mapping_confidence"] == "one-to-one ortholog"
    assert map_identifier("NCBI:GeneID:202", mapping)["mapping_confidence"] == "many-to-one ortholog"
    assert map_identifier("NCBI:GeneID:204", mapping)["mapping_confidence"] == "one-to-many ortholog"
    assert map_identifier("NCBI:GeneID:205", mapping)["mapping_confidence"] == "inferred"

    idx_out = tmp_path / "idx.tsv"
    assert build_mgigas_orthogroup_index(of_dir / "Orthogroups" / "Orthogroups.tsv", sc, idx_out) == 6

    all_out = tmp_path / "release"
    summaries = build_all_multispecies_maps(of_dir, all_out, species_codes=["cvirg"])
    assert len(summaries) == 1
    assert (all_out / "release_summary.tsv").exists()


def test_committed_ortholog_release_integrity():
    ortho_dir = ROOT / "data" / "mappings" / "orthologs"
    summary = pd.read_csv(ortho_dir / "release_summary.tsv", sep="\t")
    assert set(summary["species_code"]) == {"cvirg", "medul", "mgall", "mcali", "mtros", "mmerc", "rphil"}
    cvirg = pd.read_csv(ortho_dir / "cvirg_to_mgigas_orthofinder_v1.tsv", sep="\t", usecols=[
        "feature_id_original", "feature_id_standardized", "mapping_confidence", "orthogroup_single_copy_8sp",
    ])
    assert cvirg["feature_id_original"].is_unique
    assert set(cvirg["mapping_confidence"]).issubset({
        "one-to-one ortholog", "many-to-one ortholog", "one-to-many ortholog", "inferred"
    })
    assert (ortho_dir / "mgigas_orthogroups_v1.tsv").exists()
