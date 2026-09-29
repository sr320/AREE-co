import csv
import gzip

from aree.raw.fastq import run_mates, write_design_samplesheet, write_nfcore_samplesheet
from scripts.prepare_prjeb86618_evidence import prepare
from scripts.prepare_prjeb86618_fastqs import load_manifest


MANIFEST = "data/manifests/CGIG_STARVATION_RNASEQ_PRJEB86618_runs.tsv"


def test_manifest_is_the_family_balanced_zero_hour_contrast():
    rows = load_manifest(MANIFEST)
    assert {row["timepoint"] for row in rows} == {"h0"}
    # Three pools per diet in each family, so family is not confounded with diet.
    for condition in ("control", "starved"):
        for family in ("F11N", "F14R"):
            assert len([row for row in rows if row["condition"] == condition and row["family"] == family]) == 3


def test_every_run_is_single_end_with_one_checksummed_file():
    for row in load_manifest(MANIFEST):
        assert run_mates(row) == ("1",)
        assert row["fastq_2"] == row["fastq_2_bytes"] == row["fastq_2_md5"] == ""
        assert len(row["fastq_1_md5"]) == 32


def test_ena_aliases_agree_with_the_assigned_design():
    for row in load_manifest(MANIFEST):
        # sam_P_<family>_<3 fed | 5 starved>_<hours>_<replicate>_<ctrl | nut>-rnaseq
        code, tag = ("3", "ctrl") if row["condition"] == "control" else ("5", "nut")
        replicate = row["replicate"].split("_")[1]
        alias = "sam_P_{}_{}_0_{}_{}-rnaseq".format(row["family"], code, replicate, tag)
        assert "alias={};".format(alias) in row["assignment_evidence"]
        # Fed controls come from PRJEB86525 and starved pools from PRJEB86618.
        project = "PRJEB86525" if row["condition"] == "control" else "PRJEB86618"
        assert row["assignment_evidence"].startswith("ENA {} ".format(project))


def test_sample_sheets_are_single_end_and_carry_family(tmp_path):
    rows = load_manifest(MANIFEST)
    nfcore = tmp_path / "nfcore.csv"
    write_nfcore_samplesheet(rows, tmp_path / "fastq", nfcore)
    with nfcore.open(newline="") as handle:
        assert all(row["fastq_2"] == "" for row in csv.DictReader(handle))
    design = tmp_path / "design.csv"
    write_design_samplesheet(rows, design, extra_columns=("family",))
    with design.open(newline="") as handle:
        design_rows = list(csv.DictReader(handle))
    assert list(design_rows[0]) == ["sample", "condition", "replicate", "run_accession", "family"]
    assert len({row["sample"] for row in design_rows}) == len(design_rows)


def test_starvation_reanalysis_exports_pre_infection_evidence(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "feature_id_standardized\tlog2FoldChange\tlfcSE\tpvalue\tpadj\n"
        "NCBI:GeneID:1\t-2.0\t0.25\t0.0001\t0.001\n"
    )
    gff = tmp_path / "reference.gff.gz"
    with gzip.open(gff, "wt") as handle:
        handle.write(
            "chr1\tRefSeq\tgene\t1\t10\t.\t+\t.\t"
            "ID=gene-A;Dbxref=GeneID:1;Name=A;description=alpha;gene=A;gene_biotype=protein_coding\n"
        )
    processed = tmp_path / "processed.tsv"
    prepare(results, gff, processed, tmp_path / "mapping.tsv", tmp_path / "annotations.tsv")
    processed_rows = list(csv.DictReader(processed.open(), delimiter="\t"))
    assert processed_rows[0]["sample_comparison"] == "starved_vs_fed_ad_libitum_16d_before_oshv1"
    assert processed_rows[0]["molecular_direction"] == "down"
    assert "pre_infection_baseline" in processed_rows[0]["quality_flags"]
