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
        record = yaml.safe_load(path.read_text()) or {}
        extra[record.get("study_id", path.stem)] = record
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
            "resources": ["analysis/**"],
        },
        "website": {
            "title": "AREE",
            "site-url": "https://sr320.github.io/AREE-co/",
            "repo-url": REPO_URL,
            "repo-actions": ["source", "issue"],
            "page-footer": "Aquaculture Resilience Evidence Engine · association evidence only, not validated biomarkers",
            "navbar": {
                "left": [
                    {"text": "Progress", "href": "index.qmd"},
                    {"text": "Studies", "href": "studies.qmd"},
                    {
                        "text": "Results",
                        "menu": [
                            {"text": "Candidate scores", "href": "candidates.qmd"},
                            {"text": "Meta-analysis", "href": "meta-analysis.qmd"},
                            {"text": "Demo report", "href": "demo_report.md"},
                        ],
                    },
                    {"text": "Reanalyses", "menu": [{"text": "All reanalyses", "href": "analyses.qmd"}] + analysis_menu},
                    {"text": "Evidence cards ({})".format(len(cards)), "href": "evidence_cards/index.qmd"},
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


def progress_page(registry, extra, analyses, cards, scores, evidence):
    real = registry[~registry["data_status"].map(is_simulated)]
    reanalyzed = registry[registry["data_status"].str.contains("reanalysis_complete", na=False)]
    sha = git("rev-parse", "--short", "HEAD")
    commit_date = git("log", "-1", "--format=%cd", "--date=format:%Y-%m-%d %H:%M %Z")
    built = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    tiles = [
        ("Registered studies", len(registry)),
        ("Real public studies", len(real)),
        ("Raw-data reanalyses complete", len(reanalyzed)),
        ("Harmonized evidence records", len(evidence)),
        ("Scored candidates", len(scores)),
        ("Evidence cards", len(cards)),
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
        rows.append([rec.study_id, bp, "{} / {}".format(rec.stressor_class, rec.phenotype), str(rec.analysis_status).replace("_", " "), " · ".join(links)])
    real_table = md_table(pd.DataFrame(rows, columns=[
        "Study", "BioProject", "Stressor / phenotype", "Analysis status", "Links"
    ])) if rows else "_No real studies registered yet._"

    log = git("log", "-15", "--no-merges", "--format=%cd\t%h\t%s", "--date=short")
    activity = []
    for line in filter(None, log.splitlines()):
        date, short, subject = line.split("\t", 2)
        activity.append("- {} · [`{}`]({}/commit/{}) {}".format(date, short, REPO_URL, short, subject))

    top = scores.head(5)[["candidate_id", "score", "category", "n_studies"]].copy() if len(scores) else pd.DataFrame()
    if len(top):
        top["candidate_id"] = [
            "[{}](evidence_cards/{}.md)".format(c, safe_name(c)) if safe_name(c) in cards else c
            for c in top["candidate_id"]
        ]

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

## Real-study progress

::: {{.progress-table}}
{real_table}
:::

Simulated demo studies ({n_sim}) are listed on the [Studies](studies.qmd) page.

## Top candidates

{top}

See [all candidate scores](candidates.qmd). AREE reports associations and evidence convergence; single-study significant features are not validated biomarkers.

## Recent activity

{activity}
""".format(
        sha=sha or "unknown",
        repo=REPO_URL,
        commit_date=commit_date or "unknown",
        built=built,
        tiles=tile_md,
        real_table=real_table,
        n_sim=len(registry) - len(real),
        top=md_table(top) if len(top) else "_No candidates scored yet._",
        activity="\n".join(activity) or "_No history available._",
    )


def studies_page(registry):
    table = registry[[
        "study_id", "species", "assay_type", "tissue", "life_stage", "stressor_class",
        "phenotype", "resilience_classification", "sample_size", "data_status", "analysis_status",
    ]].copy()
    table["study_id"] = [
        "[{0}](intake/{0}.md)".format(s) if (ROOT / "reports" / "intake" / "{}.md".format(s)).exists() else s
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
        '---\ntitle: "Candidate scores"\n---\n\n'
        "Transparent candidate prioritization from `data/demo/candidate_scores.tsv` "
        "([download]({}/raw/main/data/demo/candidate_scores.tsv)). Scores rank evidence convergence; "
        "they are not validation.\n\n".format(REPO_URL)
        + md_table(table) + "\n"
    )


def meta_page(meta):
    return (
        '---\ntitle: "Meta-analysis"\n---\n\n'
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
    evidence = pd.read_csv(ROOT / "data" / "demo" / "harmonized_evidence.tsv", sep="\t", dtype=str)

    cards = {}
    for card in sorted((ROOT / "reports" / "evidence_cards").glob("*.md")):
        name = safe_name(card.stem)
        cards[name] = card.stem
        copy_md(card, out / "evidence_cards" / (name + ".md"))
    card_list = "\n".join("- [{}]({}.md)".format(orig, name) for name, orig in cards.items())
    write(out / "evidence_cards" / "index.qmd",
          '---\ntitle: "Evidence cards"\n---\n\nOne card per candidate: association evidence only, not validated biomarkers.\n\n'
          + (card_list or "_No evidence cards yet._") + "\n")

    for d in analyses:
        shutil.copytree(d, out / "analysis" / d.name)
    for md in (out / "analysis").rglob("*.md"):
        md.write_text(with_title(md.read_text()))
    for intake in sorted((ROOT / "reports" / "intake").glob("*.md")):
        copy_md(intake, out / "intake" / intake.name)
    copy_md(ROOT / "reports" / "demo_report.md", out / "demo_report.md")
    for name, _ in DOC_PAGES:
        if (ROOT / "docs" / name).exists():
            copy_md(ROOT / "docs" / name, out / "docs" / name)

    write(out / "_quarto.yml", quarto_config(analyses, cards))
    write(out / "styles.css", STYLES)
    write(out / "index.qmd", progress_page(registry, extra, analyses, cards, scores, evidence))
    write(out / "studies.qmd", studies_page(registry))
    write(out / "candidates.qmd", candidates_page(scores, cards))
    write(out / "meta-analysis.qmd", meta_page(meta))
    write(out / "analyses.qmd", analyses_page(analyses))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(ROOT / "_site_src"))
    args = parser.parse_args()
    print("site source written to {}".format(build(Path(args.out))))


if __name__ == "__main__":
    main()
