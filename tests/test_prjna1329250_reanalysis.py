import csv
import gzip

from aree.raw.fastq import write_design_samplesheet
from scripts.prepare_prjna1329250_evidence import prepare
from scripts.prepare_prjna1329250_fastqs import load_locations, load_manifest


MANIFEST = "data/manifests/CGIG_OSHV1_RNASEQ_PRJNA1329250_runs.tsv"
LOCATIONS = "data/manifests/CGIG_OSHV1_RNASEQ_PRJNA1329250_sra_locations.tsv"


def test_run_manifest_is_lineage_balanced_across_both_arms():
    rows = load_manifest(MANIFEST)
    assert {row["condition"] for row in rows} == {"control", "oshv1_uvar"}
    assert all(row["layout"] == "paired" for row in rows)
    # Both lineages contribute equally to each arm, so infection is not confounded with lineage.
    for condition, expected in (("control", 6), ("oshv1_uvar", 15)):
        arm = [row for row in rows if row["condition"] == condition]
        assert sum(row["population"] == "Midori" for row in arm) == expected
        assert sum(row["population"] == "Miyagi" for row in arm) == expected


def test_every_run_has_an_sra_location_for_checksum_verification():
    rows = load_manifest(MANIFEST)
    locations = load_locations(LOCATIONS)
    assert {row["run_accession"] for row in rows} == set(locations)
    for record in locations.values():
        assert len(record["sra_md5"]) == 32
        assert int(record["sra_size"]) > 0


def test_design_sheet_carries_population_for_blocking(tmp_path):
    rows = load_manifest(MANIFEST)
    design = tmp_path / "design.csv"
    write_design_samplesheet(rows, design, extra_columns=("population",))
    with design.open(newline="") as handle:
        design_rows = list(csv.DictReader(handle))
    assert list(design_rows[0]) == ["sample", "condition", "replicate", "run_accession", "population"]
    assert {row["population"] for row in design_rows} == {"Midori", "Miyagi"}
    assert len({row["sample"] for row in design_rows}) == len(design_rows)


def test_oshv1_reanalysis_exports_pathogen_challenge_evidence(tmp_path):
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
    assert processed_rows[0]["sample_comparison"] == "oshv1_microvariant_challenge_vs_unchallenged_control"
    assert processed_rows[0]["molecular_direction"] == "up"
    assert "pathogen_challenge_response" in processed_rows[0]["quality_flags"]
    assert "deseq2_significance_unavailable" in processed_rows[1]["quality_flags"]
    # The reanalysis must not inherit the module's default Salmon/DESeq2 version string.
    assert processed_rows[0]["analysis_method"] == "Salmon_2.3.4_tximport_1.30.0_DESeq2_1.42.1_unshrunk_effect"
