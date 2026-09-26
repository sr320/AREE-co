import typer

from aree.harmonize.processed import harmonize_demo as harmonize_demo_data
from aree.harmonize.processed import harmonize_processed
from aree.io import read_tsv
from aree.intake.registry import register_study as register_study_file
from aree.meta_analysis.random_effects import run_meta_analysis
from aree.prioritize.scoring import score_candidates
from aree.reporting.demo_report import build_demo_report as build_report
from aree.reporting.evidence_cards import build_evidence_cards as build_cards
from aree.validation.schemas import validate_study_file


app = typer.Typer(help="Aquaculture Resilience Evidence Engine CLI")


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


EVIDENCE_HELP = "Harmonized evidence TSV (default: synthetic demo table)."


@app.command("harmonize")
def harmonize(
    study: str = typer.Option(...),
    input: str = typer.Option(...),
    output: str = typer.Option(
        None, "--output", help="Evidence TSV to update. Required for studies that are not simulated demos."
    ),
    mapping: str = typer.Option(None, "--mapping"),
):
    try:
        output = harmonize_processed(study, input, output_path=output, mapping_path=mapping)
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--output") from exc
    typer.echo("harmonized evidence written to {}".format(output))


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
):
    output = run_meta_analysis(
        phenotype=phenotype, feature_type=feature_type, evidence_path=evidence, output_path=output
    )
    typer.echo("meta-analysis written to {}".format(output))
    meta = read_tsv(output)
    excluded = int(meta["n_effects_excluded"].sum())
    if excluded:
        studies = sorted({s for ids in meta["excluded_study_ids"].dropna() for s in ids.split(";") if s})
        typer.echo(
            "note: {} effects lack standard errors and were not pooled (studies: {}); "
            "{} feature groups have nothing poolable".format(
                excluded, ", ".join(studies), int((meta["pooling_status"] == "no_standard_errors").sum())
            )
        )


@app.command("score-candidates")
def score(
    evidence: str = typer.Option(None, "--evidence", help=EVIDENCE_HELP),
    phenotype: str = typer.Option(None, help="Score only evidence for this phenotype."),
    stressor: str = typer.Option(None, help="Score only evidence for this stressor."),
    output: str = typer.Option(None, "--output"),
):
    output = score_candidates(evidence_path=evidence, output_path=output, phenotype=phenotype, stressor=stressor)
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
