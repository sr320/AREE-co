from pathlib import Path
from typing import List, Optional

import typer

from aree.harmonize.processed import harmonize_demo as harmonize_demo_data
from aree.harmonize.processed import harmonize_processed
from aree.io import read_tsv
from aree.intake.registry import register_study as register_study_file
from aree.meta_analysis.random_effects import run_meta_analysis
from aree.meta_analysis.synthesis import run_real_synthesis
from aree.prioritize.scoring import score_candidates
from aree.raw.deseq2_evidence import ANALYSIS_METHOD, REFSEQ_RELEASE
from aree.raw.deseq2_evidence import export_deseq2_evidence as export_evidence
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
):
    # Cards are scored from the same phenotype-filtered evidence they show; the default
    # candidate_scores.tsv is only written by score-candidates.
    if scores:
        scores = score_candidates(evidence_path=evidence, output_path=scores, phenotype=phenotype)
    paths = build_cards(phenotype=phenotype, evidence_path=evidence, scores_path=scores, output_dir=output_dir)
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


if __name__ == "__main__":
    app()
