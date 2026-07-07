import typer

from aree.harmonize.processed import harmonize_demo as harmonize_demo_data
from aree.harmonize.processed import harmonize_processed
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


@app.command("harmonize")
def harmonize(study: str = typer.Option(...), input: str = typer.Option(...)):
    output = harmonize_processed(study, input)
    typer.echo("harmonized evidence written to {}".format(output))


@app.command("harmonize-demo")
def harmonize_demo():
    output = harmonize_demo_data()
    typer.echo("demo harmonized evidence written to {}".format(output))


@app.command("meta-analyze")
def meta_analyze(phenotype: str = typer.Option(None), feature_type: str = typer.Option(None, "--feature-type")):
    output = run_meta_analysis(phenotype=phenotype, feature_type=feature_type)
    typer.echo("meta-analysis written to {}".format(output))


@app.command("score-candidates")
def score():
    output = score_candidates()
    typer.echo("candidate scores written to {}".format(output))


@app.command("build-evidence-cards")
def build_evidence_cards(phenotype: str = typer.Option(None)):
    score_candidates()
    paths = build_cards(phenotype=phenotype)
    typer.echo("wrote {} evidence cards".format(len(paths)))


@app.command("build-demo-report")
def build_demo_report():
    output = build_report()
    typer.echo("demo report written to {}".format(output))


if __name__ == "__main__":
    app()
