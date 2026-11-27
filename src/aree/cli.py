from pathlib import Path
from typing import List, Optional

import typer

from aree.attribution import (
    add_contributor as add_contributor_entry,
    load_contributors,
    render_bibtex,
    render_citation_cff,
    render_plain_citation,
    sync_attribution as sync_attribution_files,
)
from aree.export.database import export_jsonld_catalog, export_parquet_bundle, export_schema_sql, export_sqlite
from aree.harmonize.processed import harmonize_demo as harmonize_demo_data
from aree.harmonize.processed import harmonize_processed
from aree.io import read_tsv
from aree.intake.registry import register_study as register_study_file
from aree.meta_analysis.random_effects import run_meta_analysis
from aree.meta_analysis.synthesis import run_real_synthesis
from aree.prioritize.scoring import score_candidates
from aree.raw.assay_evidence import export_assay_evidence as export_assay_table
from aree.raw.assay_evidence import validate_processed_table
from aree.raw.deseq2_evidence import ANALYSIS_METHOD, REFSEQ_RELEASE
from aree.raw.deseq2_evidence import export_deseq2_evidence as export_evidence
from aree.raw.reference import fetch_reference as fetch_ref_bundle
from aree.raw.reference import list_references as list_ref_catalog
from aree.reporting.demo_report import build_demo_report as build_report
from aree.reporting.evidence_cards import build_evidence_cards as build_cards
from aree.validation.schemas import validate_study_file


app = typer.Typer(help="Aquaculture Resilience Evidence Companion (AREC) CLI")


@app.command("validate-study")
def validate_study(path: str):
    study = validate_study_file(path)
    typer.echo("valid study registration: {}".format(study["study_id"]))


@app.command("register-study")
def register_study(path: str):
    try:
        registry = register_study_file(path)
        typer.echo("registered study in {}".format(registry))
    except ValueError as exc:
        if "Duplicate study_id" in str(exc):
            typer.echo("{}; registry unchanged".format(exc))
        else:
            raise


EVIDENCE_HELP = "Harmonized evidence TSV, or a directory of per-study TSVs (default: synthetic demo table)."


@app.command("harmonize")
def harmonize(
    study: str = typer.Option(...),
    input: str = typer.Option(...),
    output: str = typer.Option(
        None, "--output",
        help="Evidence TSV to update, or a directory that gets one STUDY_ID.tsv per study. "
        "Required for studies that are not simulated demos.",
    ),
    mapping: str = typer.Option(None, "--mapping"),
):
    try:
        output = harmonize_processed(study, input, output_path=output, mapping_path=mapping)
    except ValueError as exc:
        hint = "--output" if "explicit output path" in str(exc) else None
        raise typer.BadParameter(str(exc), param_hint=hint) from exc
    typer.echo("harmonized evidence written to {}".format(output))


@app.command("export-deseq2-evidence")
def export_deseq2_evidence(
    results: Path = typer.Option(..., help="DESeq2 all-genes TSV (feature_id_standardized, log2FoldChange, lfcSE, ...)."),
    gff: Path = typer.Option(..., help="Gzipped RefSeq GFF of the reference the reads were quantified against."),
    processed: Path = typer.Option(..., help="Processed evidence TSV to write (input to `aree harmonize`)."),
    mapping: Path = typer.Option(..., help="Exact identifier mapping TSV to write."),
    annotations: Path = typer.Option(..., help="RefSeq gene annotation TSV to write."),
    sample_comparison: str = typer.Option(..., help="Contrast label, e.g. heat_vs_control."),
    quality_flag: List[str] = typer.Option(["raw_reanalysis"], help="Study-level quality flag; repeat for several."),
    reference_release: str = typer.Option(REFSEQ_RELEASE, help="Reference assembly and annotation release."),
    analysis_method: str = typer.Option(ANALYSIS_METHOD, help="Tool chain recorded on every evidence row."),
):
    """Turn one DESeq2 contrast on a RefSeq reference into AREE processed evidence and mappings."""
    export_evidence(
        results, gff, processed, mapping, annotations,
        sample_comparison=sample_comparison, quality_flags=quality_flag,
        reference_release=reference_release, analysis_method=analysis_method,
    )


@app.command("harmonize-demo")
def harmonize_demo(output: str = typer.Option(None, "--output")):
    output = harmonize_demo_data(output_path=output)
    typer.echo("demo harmonized evidence written to {}".format(output))


@app.command("meta-analyze")
def meta_analyze(
    phenotype: str = typer.Option(None),
    feature_type: str = typer.Option(None, "--feature-type"),
    evidence: str = typer.Option(None, "--evidence", help=EVIDENCE_HELP),
    output: str = typer.Option(None, "--output"),
    method: str = typer.Option("dl", "--method", help="Random-effects estimator: dl or reml."),
    shared_control_rho: Optional[float] = typer.Option(
        None, "--shared-control-rho",
        help="Sampling correlation in [0, 1) for contrasts sharing a dependence_group (default: fail closed).",
    ),
):
    output = run_meta_analysis(
        phenotype=phenotype, feature_type=feature_type, evidence_path=evidence, output_path=output,
        method=method, shared_control_rho=shared_control_rho,
    )
    typer.echo("meta-analysis written to {}".format(output))
    meta = read_tsv(output)
    excluded = int(meta["n_effects_excluded"].sum())
    if excluded:
        studies = sorted({s for ids in meta["excluded_study_ids"].dropna() for s in ids.split(";") if s})
        no_se_groups = int((meta["pooling_status"] == "no_standard_errors").sum())
        dependent_groups = int((meta["pooling_status"] == "dependent_effects_require_covariance").sum())
        if dependent_groups and no_se_groups == 0:
            typer.echo(
                "note: {} effects were excluded from pooling (studies: {}); "
                "{} feature groups require covariance modeling (--shared-control-rho)".format(
                    excluded, ", ".join(studies), dependent_groups
                )
            )
        else:
            typer.echo(
                "note: {} effects lack standard errors and were not pooled (studies: {}); "
                "{} feature groups have nothing poolable".format(
                    excluded, ", ".join(studies), no_se_groups
                )
            )


@app.command("synthesize-real")
def synthesize_real(
    evidence: str = typer.Option(None, "--evidence", help=EVIDENCE_HELP),
    config: str = typer.Option(None, "--config", help="Path to synthesis_strata.yaml."),
    data_dir: str = typer.Option(None, "--data-dir", help="Output directory for per-stratum pooled TSVs."),
    reports_dir: str = typer.Option(None, "--reports-dir", help="Output directory for synthesis reports and summaries."),
    method: str = typer.Option(None, "--method", help="Random-effects estimator override (dl or reml)."),
    shared_control_rho: Optional[float] = typer.Option(
        None, "--shared-control-rho", help="Override shared-control correlation in [0, 1).",
    ),
):
    """Run comparability-reviewed real-study synthesis across curated strata."""
    report = run_real_synthesis(
        evidence_path=evidence, config_path=config, data_dir=data_dir,
        reports_dir=reports_dir, method=method, shared_control_rho=shared_control_rho,
    )
    typer.echo("real synthesis report written to {}".format(report))


@app.command("score-candidates")
def score(
    evidence: str = typer.Option(None, "--evidence", help=EVIDENCE_HELP),
    phenotype: str = typer.Option(None, help="Score only evidence for this phenotype."),
    stressor: str = typer.Option(None, help="Score only evidence for this stressor."),
    output: str = typer.Option(None, "--output"),
    method: str = typer.Option("dl", "--method", help="Heterogeneity estimator: dl or reml."),
    shared_control_rho: Optional[float] = typer.Option(
        None, "--shared-control-rho",
        help="Sampling correlation in [0, 1) for contrasts sharing a dependence_group.",
    ),
):
    output = score_candidates(
        evidence_path=evidence, output_path=output, phenotype=phenotype, stressor=stressor,
        method=method, shared_control_rho=shared_control_rho,
    )
    typer.echo("candidate scores written to {}".format(output))


@app.command("build-evidence-cards")
def build_evidence_cards(
    phenotype: str = typer.Option(None),
    evidence: str = typer.Option(None, "--evidence", help=EVIDENCE_HELP),
    scores: str = typer.Option(
        None, "--scores", help="Also write the scores behind the cards to this TSV (default: score in memory)."
    ),
    output_dir: str = typer.Option(None, "--output-dir"),
    plots: bool = typer.Option(False, "--plots", help="Render standalone SVG forest plots alongside evidence cards."),
):
    # Cards are scored from the same phenotype-filtered evidence they show; the default
    # candidate_scores.tsv is only written by score-candidates.
    if scores:
        scores = score_candidates(evidence_path=evidence, output_path=scores, phenotype=phenotype)
    paths = build_cards(
        phenotype=phenotype,
        evidence_path=evidence,
        scores_path=scores,
        output_dir=output_dir,
        include_plots=plots,
    )
    typer.echo("wrote {} evidence cards".format(len(paths)))


@app.command("build-demo-report")
def build_demo_report(
    evidence: str = typer.Option(None, "--evidence", help=EVIDENCE_HELP),
    scores: str = typer.Option(None, "--scores"),
    registry: str = typer.Option(None, "--registry"),
    output: str = typer.Option(None, "--output"),
):
    output = build_report(output_path=output, registry_path=registry, evidence_path=evidence, scores_path=scores)
    typer.echo("demo report written to {}".format(output))


@app.command("fetch-reference")
def fetch_reference(
    reference: Optional[str] = typer.Option(
        None, "--reference", help="RefSeq accession, species name, or code (default: catalog default)."
    ),
    output_dir: Optional[Path] = typer.Option(
        None, "--output-dir", help="Destination directory (default: data/reference/<accession>)."
    ),
    catalog: Optional[Path] = typer.Option(None, "--catalog", help="Path to references.yaml."),
    mirror_dir: Optional[Path] = typer.Option(
        None, "--mirror-dir", help="Optional local mirror directory holding verified source files."
    ),
    include_genome: bool = typer.Option(False, "--include-genome", help="Also fetch genome FASTA and build gentrome/decoys."),
    list_only: bool = typer.Option(False, "--list", help="List cataloged references and exit."),
):
    """Fetch and verify a versioned RefSeq reference bundle from the multi-species catalog."""
    if list_only:
        for spec in list_ref_catalog(catalog):
            typer.echo("{}\t{}\t{}\t{}".format(
                spec["accession"], spec.get("species_code", ""), spec["species"], spec["annotation_release"]
            ))
        return
    manifest = fetch_ref_bundle(
        identifier=reference,
        output_dir=output_dir,
        include_genome=include_genome,
        catalog_path=catalog,
        mirror_dir=mirror_dir,
    )
    typer.echo("fetched {} ({}) -> {} ({} tx2gene mappings)".format(
        manifest["accession"], manifest["species"], manifest["manifest_path"], manifest["tx2gene_count"]
    ))


@app.command("export-assay-evidence")
def export_assay_evidence_cmd(
    results: Path = typer.Option(..., "--results", exists=True, help="Differential/quantitative assay TSV."),
    processed: Path = typer.Option(..., "--processed", help="Output processed evidence TSV path."),
    mapping: Path = typer.Option(..., "--mapping", help="Output identifier mapping TSV path."),
    annotations: Path = typer.Option(..., "--annotations", help="Output feature annotations TSV path."),
    sample_comparison: str = typer.Option(..., "--sample-comparison", help="Contrast label."),
    feature_type: str = typer.Option(..., "--feature-type", help="Assay feature type (cpg_site, dmr, protein, metabolite, etc.)."),
    effect_size_type: str = typer.Option(..., "--effect-size-type", help="Effect size unit label."),
    analysis_method: str = typer.Option(..., "--analysis-method", help="Analysis method string with tool versions."),
    reference_release: str = typer.Option("assay_feature_catalog_v1", "--reference-release"),
    direction_mode: str = typer.Option("abundance", "--direction-mode", help="Direction convention: abundance, expression, or methylation."),
    ortholog_reference: str = typer.Option("none", "--ortholog-reference"),
    mapping_release: str = typer.Option("direct_assay_id_v1", "--mapping-release"),
    mapping_evidence: str = typer.Option("Assay-native standardized feature identifier", "--mapping-evidence"),
    quality_flag: Optional[List[str]] = typer.Option(None, "--quality-flag", help="Study-level quality flag; repeat as needed."),
):
    """Export WGBS/methylation, proteomics, or metabolomics assay tables to AREE evidence tables."""
    report = export_assay_table(
        results_path=results,
        processed_path=processed,
        mapping_path=mapping,
        annotation_path=annotations,
        sample_comparison=sample_comparison,
        feature_type=feature_type,
        effect_size_type=effect_size_type,
        analysis_method=analysis_method,
        quality_flags=quality_flag or [],
        reference_release=reference_release,
        direction_mode=direction_mode,
        ortholog_reference=ortholog_reference,
        mapping_release=mapping_release,
        mapping_evidence=mapping_evidence,
    )
    typer.echo("exported {} {} features ({} FDR < 0.05) -> {}".format(
        report["n_rows"], feature_type, report["n_fdr_lt_0_05"], processed
    ))


@app.command("validate-processed")
def validate_processed_cmd(
    processed: Path = typer.Option(..., "--processed", exists=True, help="Processed evidence TSV to validate."),
    mapping: Optional[Path] = typer.Option(None, "--mapping", exists=True, help="Optional identifier mapping TSV."),
    annotation: Optional[Path] = typer.Option(None, "--annotation", exists=True, help="Optional annotation TSV."),
    out: Optional[Path] = typer.Option(None, "--out", help="Optional JSON report path."),
):
    """Validate exported processed evidence, identifier mappings, and annotations."""
    report = validate_processed_table(
        processed_path=processed,
        mapping_path=mapping,
        annotation_path=annotation,
        out_path=out,
    )
    typer.echo("validated {} rows across {} unique features ({} FDR < 0.05)".format(
        report["n_rows"], report["n_unique_features"], report["n_fdr_lt_0_05"]
    ))


@app.command("export-db")
def export_db_cmd(
    format: str = typer.Option("sqlite", "--format", help="Export target: sqlite, parquet, or jsonld."),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        help="Destination file (sqlite/jsonld) or directory (parquet). Defaults under reports/exports/.",
    ),
    schema_out: Optional[Path] = typer.Option(None, "--schema-out", help="Optional path to also write schema.sql."),
    include_simulated: bool = typer.Option(False, "--include-simulated", help="Include simulated study rows."),
):
    """Export curated studies, references, evidence, synthesis, and ortholog summaries to DB formats."""
    fmt = format.strip().lower()
    if schema_out is not None:
        export_schema_sql(schema_out)
    if fmt == "sqlite":
        dest = output or Path("results/aree_evidence.sqlite")
        info = export_sqlite(dest, include_simulated=include_simulated)
        typer.echo("exported SQLite DB to {} ({} bytes)".format(info["output_path"], info["bytes"]))
    elif fmt == "parquet":
        dest = output or Path("results/parquet_bundle")
        info = export_parquet_bundle(dest, include_simulated=include_simulated)
        typer.echo("exported Parquet Data Package to {}".format(info["datapackage_path"]))
    elif fmt == "jsonld":
        dest = output or Path("reports/exports/aree_catalog.jsonld")
        info = export_jsonld_catalog(dest)
        typer.echo("exported JSON-LD catalog to {} ({} datasets)".format(info["output_path"], info["n_datasets"]))
    else:
        raise typer.BadParameter("Unsupported --format {!r}; expected sqlite, parquet, or jsonld".format(format))


@app.command("add-contributor")
def add_contributor_cmd(
    given_name: str = typer.Option(..., "--given-name", help="Contributor's given (first) name."),
    family_name: str = typer.Option(..., "--family-name", help="Contributor's family (last) name."),
    affiliation: Optional[str] = typer.Option(None, "--affiliation", help="Institution or organization."),
    orcid: Optional[str] = typer.Option(None, "--orcid", help="ORCID iD (e.g. 0000-0000-0000-0000)."),
    github: Optional[str] = typer.Option(None, "--github", help="GitHub username (with or without @)."),
    email: Optional[str] = typer.Option(None, "--email", help="Contact email address."),
    role: Optional[List[str]] = typer.Option(
        None,
        "--role",
        help="Contributor role (CRediT role recommended; repeat --role for multiple).",
    ),
    notes: Optional[str] = typer.Option(None, "--notes", help="Optional note on specific contributions."),
    acknowledge_only: bool = typer.Option(
        False,
        "--acknowledge-only",
        help="List in CONTRIBUTORS.md without adding to CITATION.cff / BibTeX author list.",
    ),
    update: bool = typer.Option(
        False,
        "--update",
        help="Update an existing contributor entry matched by ORCID, GitHub, email, or name.",
    ),
    contributors_file: Optional[Path] = typer.Option(
        None, "--contributors-file", help="Path to contributors.yaml (defaults to project root)."
    ),
    sync: bool = typer.Option(
        True, "--sync/--no-sync", help="Regenerate CITATION.cff and CONTRIBUTORS.md after updating contributors.yaml."
    ),
):
    """Add (or update) a contributor in contributors.yaml and sync citation files."""
    try:
        saved_path, action = add_contributor_entry(
            given_name=given_name,
            family_name=family_name,
            affiliation=affiliation,
            orcid=orcid,
            github=github,
            email=email,
            roles=role or [],
            include_in_citation=not acknowledge_only,
            notes=notes,
            update=update,
            path=contributors_file,
        )
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc

    typer.echo("{} {} {} in {}".format(action, given_name.strip(), family_name.strip(), saved_path))
    if sync:
        info = sync_attribution_files(contributors_path=saved_path)
        typer.echo(
            "synced {} contributor(s) ({} cited) -> {}, {}".format(
                info["n_contributors"], info["n_cited_authors"], info["citation_path"], info["markdown_path"]
            )
        )


@app.command("sync-attribution")
def sync_attribution_cmd(
    contributors_file: Optional[Path] = typer.Option(
        None, "--contributors-file", help="Path to contributors.yaml (defaults to project root)."
    ),
    citation_out: Optional[Path] = typer.Option(None, "--citation-out", help="Output path for CITATION.cff."),
    markdown_out: Optional[Path] = typer.Option(None, "--markdown-out", help="Output path for CONTRIBUTORS.md."),
):
    """Validate contributors.yaml and regenerate CITATION.cff and CONTRIBUTORS.md."""
    info = sync_attribution_files(
        contributors_path=contributors_file,
        citation_path=citation_out,
        markdown_path=markdown_out,
    )
    typer.echo(
        "synced {} contributor(s) ({} cited) -> {}, {}".format(
            info["n_contributors"], info["n_cited_authors"], info["citation_path"], info["markdown_path"]
        )
    )


@app.command("cite")
def cite_cmd(
    format: str = typer.Option("text", "--format", help="Citation output format: text, bibtex, or cff."),
    contributors_file: Optional[Path] = typer.Option(
        None, "--contributors-file", help="Path to contributors.yaml (defaults to project root)."
    ),
):
    """Print the current AREC citation in plain text, BibTeX, or CFF format."""
    data = load_contributors(contributors_file)
    fmt = format.strip().lower()
    if fmt == "text":
        typer.echo(render_plain_citation(data))
    elif fmt == "bibtex":
        typer.echo(render_bibtex(data).rstrip())
    elif fmt == "cff":
        typer.echo(render_citation_cff(data).rstrip())
    else:
        raise typer.BadParameter("Unsupported --format {!r}; expected text, bibtex, or cff".format(format))


if __name__ == "__main__":
    app()
