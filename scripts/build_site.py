"""Assemble the AREE progress website from committed registry, reports, and docs.

Writes a Quarto website project to ``_site_src/``; render it with
``quarto render _site_src`` (output lands in ``_site/``). The site is rebuilt
from whatever is committed, so new studies, reanalyses, and evidence cards
appear automatically once they reach ``main``.
"""

import argparse
import datetime
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import yaml

from aree.reporting.site_findings import load_findings, write_real_findings
from aree.reporting.evidence_cards import card_filename
from aree.intake.registry import flatten_for_registry, REGISTRY_COLUMNS
from aree.validation.schemas import validate_study_file

ROOT = Path(__file__).resolve().parents[1]
REPO_URL = "https://github.com/sr320/AREE-co"

STYLES = """
.stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(10rem, 1fr)); gap: .75rem; margin: 1.25rem 0 1.5rem; }
.stat { border: 1px solid rgba(127,127,127,.3); border-left: 4px solid var(--bs-primary); border-radius: .4rem; padding: .7rem .9rem; }
.stat-value { font-size: 1.9rem; font-weight: 600; line-height: 1.1; font-variant-numeric: tabular-nums; }
.stat-label { font-size: .85rem; opacity: .75; }
.progress-table td:last-child { white-space: nowrap; }
main table { display: block; max-width: 100%; overflow-x: auto; }
"""

DOC_PAGES = [
    ("why-this-matters.md", "Why this matters"),
    ("methods.md", "Methods"),
    ("interpreting-evidence.md", "Interpreting evidence"),
    ("adding-a-study.md", "Adding a study"),
    ("species-and-genomes.md", "Species and genomes"),
    ("architecture.md", "Architecture"),
    ("design.md", "Design"),
    ("governance.md", "Governance"),
    ("roadmap.md", "Roadmap"),
]


def safe_name(value):
    """File-safe stem; ':' in a relative link would be parsed as a URL scheme."""
    return str(value).replace(":", "_").replace("/", "_")


def md_escape(value):
    if pd.isna(value):
        return ""
    return str(value).replace("|", "\\|").replace("\n", " ")


def md_table(df):
    header = "| " + " | ".join(df.columns) + " |"
    rule = "| " + " | ".join("---" for _ in df.columns) + " |"
    rows = ["| " + " | ".join(md_escape(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join([header, rule] + rows)


def git(*args):
    try:
        return subprocess.run(
            ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def is_simulated(status):
    return str(status).startswith("simulated")


def load_registry():
    registry = pd.read_csv(ROOT / "registry" / "study_registry.csv", dtype=str)
    extra = {}
    for path in sorted((ROOT / "registry" / "studies").glob("*.yaml")):
        record = validate_study_file(path)
        extra[record.get("study_id", path.stem)] = record
    expected = pd.DataFrame([flatten_for_registry(study) for study in extra.values()], columns=REGISTRY_COLUMNS)
    def normalized(frame):
        return frame[REGISTRY_COLUMNS].fillna("").astype(str).sort_values("study_id").reset_index(drop=True)
    if not set(REGISTRY_COLUMNS).issubset(registry.columns) or not normalized(expected).equals(normalized(registry)):
        raise ValueError("Registry CSV is stale or inconsistent with study YAML files; regenerate it before building the site")
    return registry, extra


def analysis_dirs():
    base = ROOT / "reports" / "analysis"
    return sorted(p for p in base.iterdir() if p.is_dir()) if base.exists() else []


def bioproject_for(study_id, extra):
    return (extra.get(study_id, {}).get("accessions") or {}).get("bioproject", "")


def with_title(text):
    """Promote a leading '# Heading' to front matter so Quarto uses it as the page title."""
    lines = text.splitlines()
    if lines and lines[0].startswith("# ") and not text.startswith("---"):
        title = lines[0][2:].strip().replace('"', '\\"')
        return '---\ntitle: "{}"\n---\n'.format(title) + "\n".join(lines[1:]) + "\n"
    return text


def copy_md(src, dst):
    write(dst, with_title(src.read_text()))


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def quarto_config(analyses, cards):
    analysis_menu = [
        {"text": d.name, "href": "analysis/{}/analysis_report.md".format(d.name)}
        for d in analyses
        if (d / "analysis_report.md").exists()
    ]
    config = {
        "project": {
            "type": "website",
            "output-dir": "../_site",
            "resources": ["analysis/**", "downloads/**"],
        },
        "website": {
            "title": "AREE",
            "site-url": "https://sr320.github.io/AREE-co/",
            "repo-url": REPO_URL,
            # Pages are generated outside the checkout, so Quarto's default source URLs
            # point at files that do not exist in the repository. Keep the repository icon
            # and issue action; cards carry explicit links to their committed inputs.
            "repo-actions": ["issue"],
            "page-footer": "Aquaculture Resilience Evidence Engine · association evidence only, not validated biomarkers",
            "navbar": {
                "left": [
                    {"text": "Progress", "href": "index.qmd"},
                    {"text": "Studies", "href": "studies.qmd"},
                    {"text": "Real-study findings", "href": "candidates.qmd"},
                    {"text": "Reanalyses", "menu": [{"text": "All reanalyses", "href": "analyses.qmd"}] + analysis_menu},
                    {"text": "Real evidence cards ({})".format(len(cards)), "href": "evidence_cards/index.qmd"},
                    {"text": "Synthetic demo", "href": "demo/index.qmd"},
                    {
                        "text": "Docs",
                        "menu": [{"text": title, "href": "docs/" + name} for name, title in DOC_PAGES],
                    },
                ],
                "right": [{"icon": "github", "href": REPO_URL}],
            },
        },
        "format": {"html": {"theme": {"light": "cosmo", "dark": "darkly"}, "css": "styles.css", "toc": True}},
    }
    return yaml.safe_dump(config, sort_keys=False)


def progress_page(registry, extra, summaries, cards, highlights, evidence):
    real = registry[~registry["data_status"].map(is_simulated)]
    reanalyzed = real[real["data_status"].str.contains("reanalysis_complete", na=False)]
    sha = git("rev-parse", "--short", "HEAD")
    commit_date = git("log", "-1", "--format=%cd", "--date=format:%Y-%m-%d %H:%M %Z")
    built = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    tiles = [
        ("Real public BioProjects", real["study_id"].map(lambda study: bioproject_for(study, extra)).nunique()),
        ("Registered real-study contrasts", len(real)),
        ("BioProjects reanalysed", reanalyzed["study_id"].map(lambda study: bioproject_for(study, extra)).nunique()),
        ("Contrasts with harmonized evidence", int(summaries["records"].gt(0).sum())),
        ("Real harmonized evidence records", len(evidence)),
        ("Real evidence cards", len(cards)),
    ]
    tile_md = "\n".join(
        '<div class="stat"><div class="stat-value">{:,}</div><div class="stat-label">{}</div></div>'.format(value, label)
        for label, value in tiles
    )

    rows = []
    for rec in real.itertuples(index=False):
        bp = bioproject_for(rec.study_id, extra)
        analysis = ROOT / "reports" / "analysis" / str(bp) / "analysis_report.md"
        intake = ROOT / "reports" / "intake" / "{}.md".format(rec.study_id)
        links = []
        if isinstance(rec.doi, str) and rec.doi:
            links.append("[paper](https://doi.org/{})".format(rec.doi))
        if intake.exists():
            links.append("[intake](intake/{}.md)".format(rec.study_id))
        if bp and analysis.exists():
            links.append("[reanalysis](analysis/{}/analysis_report.md)".format(bp))
        links.insert(0, "[design and limitations](studies/{}.qmd)".format(rec.study_id))
        rows.append([rec.study_id, bp, "{} / {}".format(rec.stressor_class, rec.phenotype), str(rec.analysis_status).replace("_", " "), " · ".join(links)])
    real_table = md_table(pd.DataFrame(rows, columns=[
        "Study", "BioProject", "Stressor / phenotype", "Analysis status", "Links"
    ])) if rows else "_No real studies registered yet._"

    log = git("log", "-15", "--no-merges", "--format=%cd\t%h\t%s", "--date=short")
    activity = []
    for line in filter(None, log.splitlines()):
        date, short, subject = line.split("\t", 2)
        activity.append("- {} · [`{}`]({}/commit/{}) {}".format(date, short, REPO_URL, short, subject))

    examples = []
    for rec in highlights.groupby("study_id", sort=True).head(1).itertuples(index=False):
        examples.append([bioproject_for(rec.study_id, extra),
                         "[{}](evidence_cards/{})".format(rec.feature_id_standardized, card_filename(rec.feature_id_standardized)),
                         rec.sample_comparison, rec.description, rec.inference_status])
    top = md_table(pd.DataFrame(examples, columns=["BioProject", "Feature / real evidence card", "Contrast", "RefSeq description", "Interpretation"])) if examples else "_No real findings available yet._"

    return """---
title: "Aquaculture Resilience Evidence Engine"
subtitle: "Results to date"
toc: false
page-layout: full
---

AREE converts public aquaculture omics studies into harmonized, comparable evidence for resilience-associated biomarker discovery. This page is rebuilt automatically every time results are pushed to `main`.

::: {{.callout-tip appearance="simple"}}
Last updated from commit [`{sha}`]({repo}/commit/{sha}) ({commit_date}); site built {built}.
:::

```{{=html}}
<div class="stats">
{tiles}
</div>
```

[Browse real-study findings](candidates.qmd) · [Browse real evidence cards](evidence_cards/index.qmd)

Findings retain study-specific contrasts, gene descriptions, uncertainty and design limitations. They describe molecular associations, not validated predictors of resilience. The shared DECICOMP controls are not independent replication, and PRJNA735889's tank-dependent analysis remains exploratory.

Counts distinguish public BioProjects from registered contrasts. PRJNA913164 contributes two contrasts with shared controls; neither the contrast count nor the BioProject count establishes independent replication.

[How do expression findings relate to resilience?](docs/interpreting-evidence.md#how-expression-findings-relate-to-resilience) Learn what the comparison supports, why higher expression does not necessarily mean greater resilience, and what validation is needed.

## Real-study progress

::: {{.progress-table}}
{real_table}
:::

## One displayed finding per harmonized study

{top}

Each row is the first finding from its own study's selection; rows are ordered by study ID, not a global biomarker score. See [all displayed findings and full-study downloads](candidates.qmd).

## Workflow demonstration

The [synthetic demo](demo/index.qmd) contains simulated studies, scores, meta-analysis and cards for testing the workflow. It is separate from the real findings above.

## Recent activity

{activity}
""".format(
        sha=sha or "unknown",
        repo=REPO_URL,
        commit_date=commit_date or "unknown",
        built=built,
        tiles=tile_md,
        real_table=real_table,
        top=top,
        activity="\n".join(activity) or "_No history available._",
    )


def studies_page(registry):
    table = registry[[
        "study_id", "species", "assay_type", "tissue", "life_stage", "stressor_class",
        "phenotype", "resilience_classification", "sample_size", "data_status", "analysis_status",
    ]].copy()
    table["study_id"] = [
        "[{0}](studies/{0}.qmd)".format(s)
        for s in table["study_id"]
    ]
    return '---\ntitle: "Registered studies"\n---\n\nFrom `registry/study_registry.csv`. Studies whose data status begins with `simulated` are MVP demo data.\n\n' + md_table(table) + "\n"


def candidates_page(scores, cards):
    table = scores.copy()
    table["candidate_id"] = [
        "[{}](evidence_cards/{}.md)".format(c, safe_name(c)) if safe_name(c) in cards else c
        for c in table["candidate_id"]
    ]
    return (
        '---\ntitle: "Simulated demo candidate scores"\n---\n\n'
        "**Simulated data:** These scores demonstrate the workflow; they are not results from real public studies.\n\n"
        "Transparent candidate prioritization from `data/demo/candidate_scores.tsv` "
        "([download]({}/raw/main/data/demo/candidate_scores.tsv)). Scores rank evidence convergence; "
        "they are not validation.\n\n".format(REPO_URL)
        + md_table(table) + "\n"
    )


def study_detail_page(study):
    parts = ['---\ntitle: "{}"\n---\n'.format(study["study_id"])]
    if is_simulated(study["data_availability"]["status"]):
        parts.append("**Simulated demo study:** All observations are synthetic.\n")
    parts.append("**Primary-analysis sample count:** {}. The counted unit and independent experimental units "
                 "are described under Biological replication; counts may represent animals, pools or libraries.\n"
                 .format(study["sample_size"]))
    if "tank_clustering_unmodeled" in study["analysis_status"]:
        parts.append("::: {.callout-warning}\n**Exploratory analysis: tank clustering is not modeled.** "
                     "The 50 oysters are subsamples from ten tanks (five per pH band). "
                     "Reported p-values, FDR values and standard errors assume independent oysters "
                     "and may overstate precision. A tank-aware reanalysis is required before "
                     "inferential use. These p-values receive no significance reward in candidate scoring.\n:::\n")
    for label, field in [("Treatment", "experimental_treatment"), ("Control", "control_condition"),
                         ("Biological replication", "biological_replication"),
                         ("Phenotype interpretation", "phenotype_direction"),
                         ("Analysis status", "analysis_status"), ("Quality control and results", "quality_control_status")]:
        parts.append("## {}\n\n{}\n".format(label, study[field]))
    parts.append("## Limitations\n")
    parts.extend("- " + limitation for limitation in study["limitations"])
    parts.append("\n## Sources\n")
    parts.extend("- [Source {}]({})".format(i, url)
                 for i, url in enumerate(study["provenance"]["source_links"], 1))
    return "\n".join(parts) + "\n"


def meta_page(meta):
    return (
        '---\ntitle: "Simulated demo meta-analysis"\n---\n\n'
        "**Simulated data:** These effects are synthetic; they are not results from real public studies.\n\n"
        "Random-effects meta-analysis results from `data/demo/meta_analysis.tsv` "
        "([download]({}/raw/main/data/demo/meta_analysis.tsv)).\n\n".format(REPO_URL)
        + md_table(meta) + "\n"
    )


def analyses_page(analyses):
    parts = ['---\ntitle: "Raw-data reanalyses"\n---\n']
    if not analyses:
        parts.append("_No reanalyses committed yet._")
    for d in analyses:
        parts.append("## {}\n".format(d.name))
        reports = sorted(d.rglob("*.md"))
        for r in reports:
            parts.append("- [{}](analysis/{}/{})".format(r.stem.replace("_", " "), d.name, r.relative_to(d)))
        files = sorted(p for p in d.rglob("*") if p.is_file() and p.suffix != ".md")
        if files:
            parts.append("\n**Data files**\n")
            for f in files:
                rel = f.relative_to(d)
                parts.append("- [`{}`](analysis/{}/{}) ({:,} KB)".format(rel, d.name, rel, max(1, f.stat().st_size // 1024)))
        parts.append("")
    return "\n".join(parts) + "\n"


def build(out):
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    registry, extra = load_registry()
    analyses = analysis_dirs()
    scores = pd.read_csv(ROOT / "data" / "demo" / "candidate_scores.tsv", sep="\t", dtype=str)
    meta = pd.read_csv(ROOT / "data" / "demo" / "meta_analysis.tsv", sep="\t", dtype=str)
    real_evidence, summaries, highlights = load_findings(ROOT, extra)
    real_cards = write_real_findings(out, real_evidence, summaries, highlights, extra, md_table, write)

    cards = {}
    for card in sorted((ROOT / "reports" / "evidence_cards").glob("*.md")):
        name = safe_name(card.stem)
        cards[name] = card.stem
        text = with_title(card.read_text())
        # Put the simulation label after front matter so it is visible on each card.
        title_end = text.index("---", 3) + 3
        text = text[:title_end] + "\n\n**Simulated demo evidence:** This card does not describe real-study findings.\n" + text[title_end:]
        write(out / "demo" / "evidence_cards" / (name + ".md"), text)
        # Keep previously published demo-card URLs usable after the move.
        write(out / "evidence_cards" / (name + ".qmd"),
              '---\ntitle: "Synthetic demo card moved"\n---\n\n'
              'This is a synthetic workflow example. [View the demo card](../demo/evidence_cards/{}.md) '
              'or [browse real evidence cards](index.qmd).\n'.format(name))
    card_list = "\n".join("- [{}]({}.md)".format(orig, name) for name, orig in cards.items())
    write(out / "demo" / "evidence_cards" / "index.qmd",
          '---\ntitle: "Simulated demo evidence cards"\n---\n\nThese cards use synthetic observations, not real-study findings. Association evidence only, not validated biomarkers.\n\n'
          + (card_list or "_No evidence cards yet._") + "\n")

    for d in analyses:
        shutil.copytree(d, out / "analysis" / d.name)
    for md in (out / "analysis").rglob("*.md"):
        md.write_text(with_title(md.read_text()))
    for intake in sorted((ROOT / "reports" / "intake").glob("*.md")):
        copy_md(intake, out / "intake" / intake.name)
    copy_md(ROOT / "reports" / "demo_report.md", out / "demo" / "demo_report.md")
    write(out / "demo_report.md", '---\ntitle: "Synthetic demo report moved"\n---\n\n'
          '[View the synthetic demo report](demo/demo_report.md) or [browse real-study findings](candidates.qmd).\n')
    for name, _ in DOC_PAGES:
        if (ROOT / "docs" / name).exists():
            copy_md(ROOT / "docs" / name, out / "docs" / name)
    for study_id, study in extra.items():
        base = out / "demo" if is_simulated(study["data_availability"]["status"]) else out
        write(base / "studies" / (study_id + ".qmd"), study_detail_page(study))

    write(out / "_quarto.yml", quarto_config(analyses, real_cards))
    write(out / "styles.css", STYLES)
    write(out / "index.qmd", progress_page(registry, extra, summaries, real_cards, highlights, real_evidence))
    write(out / "studies.qmd", studies_page(registry[~registry["data_status"].map(is_simulated)]))
    write(out / "demo/studies.qmd", studies_page(registry[registry["data_status"].map(is_simulated)]))
    write(out / "demo/candidates.qmd", candidates_page(scores, cards))
    write(out / "demo/meta-analysis.qmd", meta_page(meta))
    write(out / "demo/index.qmd", '---\ntitle: "Synthetic workflow demonstration"\n---\n\n'
          '**All evidence, scores, effects and cards in this section are simulated.** They are not findings from the real public studies.\n\n'
          '- [Simulated studies](studies.qmd)\n- [Demo candidate scores](candidates.qmd)\n'
          '- [Demo meta-analysis](meta-analysis.qmd)\n- [Demo evidence cards](evidence_cards/index.qmd)\n'
          '- [Demo report](demo_report.md)\n\n[Return to real-study findings](../candidates.qmd)\n')
    write(out / "meta-analysis.qmd", '---\ntitle: "Cross-study comparison boundaries"\n---\n\n'
          'The primary site reports [real study-level findings](candidates.qmd) and [real evidence cards](evidence_cards/index.qmd). '
          'No cross-study pooled effects are presented here: acute exposures, population-selection contrasts, tissues and life stages require comparability review. '
          'Shared controls and exploratory tank-dependent inference also constrain pooling. '
          'The earlier synthetic results are available in [demo meta-analysis](demo/meta-analysis.qmd).\n')
    write(out / "analyses.qmd", analyses_page(analyses))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(ROOT / "_site_src"))
    args = parser.parse_args()
    print("site source written to {}".format(build(Path(args.out))))


if __name__ == "__main__":
    main()
