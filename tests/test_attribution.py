from pathlib import Path

import jsonschema
import pytest
import yaml
from typer.testing import CliRunner

from aree.attribution import (
    add_contributor,
    default_citation_path,
    default_contributors_path,
    default_markdown_path,
    load_contributors,
    normalize_orcid,
    render_bibtex,
    render_citation_cff,
    render_contributors_markdown,
    render_plain_citation,
    save_contributors,
    sync_attribution,
)
from aree.cli import app


runner = CliRunner()


def _sample_data():
    return {
        "project": {
            "title": "Aquaculture Resilience Evidence Companion",
            "abbreviation": "AREC",
            "version": "0.1.0",
            "date_released": "2026-11-26",
            "repository_code": "https://github.com/sr320/AREE-co",
            "url": "https://sr320.github.io/AREE-co/",
            "license": "MIT",
            "message": "If you use AREC, please cite this repository.",
        },
        "contributors": [
            {
                "given_name": "Steven",
                "family_name": "Roberts",
                "affiliation": "University of Washington",
                "github": "sr320",
                "roles": ["Conceptualization", "Software"],
                "include_in_citation": True,
            }
        ],
    }


def test_normalize_orcid():
    assert normalize_orcid("0000-0001-8302-1138") == "https://orcid.org/0000-0001-8302-1138"
    assert normalize_orcid("https://orcid.org/0000-0001-8302-113X") == "https://orcid.org/0000-0001-8302-113X"
    assert normalize_orcid("") is None
    assert normalize_orcid(None) is None
    with pytest.raises(ValueError, match="Invalid ORCID"):
        normalize_orcid("not-an-orcid")


def test_add_update_and_sync_roundtrip(tmp_path):
    contrib_path = tmp_path / "contributors.yaml"
    cff_path = tmp_path / "CITATION.cff"
    md_path = tmp_path / "CONTRIBUTORS.md"
    save_contributors(_sample_data(), contrib_path)

    _, action = add_contributor(
        given_name="Jane",
        family_name="Doe",
        affiliation="Marine Lab",
        orcid="0000-0002-1825-0097",
        github="@janedoe",
        email="jane@example.org",
        roles=["Data curation", "Validation"],
        notes="Curated thermal tolerance studies.",
        path=contrib_path,
    )
    assert action == "added"

    # Duplicate without update=True raises
    with pytest.raises(ValueError, match="already exists"):
        add_contributor(given_name="Jane", family_name="Doe", path=contrib_path)

    # Updating merges new roles and preserves existing ones
    _, action = add_contributor(
        given_name="Jane",
        family_name="Doe",
        roles=["Software"],
        update=True,
        path=contrib_path,
    )
    assert action == "updated"

    # Acknowledge-only contributor is excluded from citation outputs
    add_contributor(
        given_name="Alex",
        family_name="Helper",
        roles=["Resources"],
        include_in_citation=False,
        path=contrib_path,
    )

    info = sync_attribution(contrib_path, cff_path, md_path)
    assert info["n_contributors"] == 3
    assert info["n_cited_authors"] == 2

    cff = yaml.safe_load(cff_path.read_text())
    assert cff["cff-version"] == "1.2.0"
    assert [a["family-names"] for a in cff["authors"]] == ["Roberts", "Doe"]
    assert cff["authors"][1]["orcid"] == "https://orcid.org/0000-0002-1825-0097"

    data = load_contributors(contrib_path)
    assert data["contributors"][1]["roles"] == ["Data curation", "Validation", "Software"]

    plain = render_plain_citation(data)
    assert "Roberts, S. & Doe, J. (2026)." in plain
    assert "Helper" not in plain

    bib = render_bibtex(data)
    assert "@software{arec_2026," in bib
    assert "Roberts, Steven and Doe, Jane" in bib
    assert "Helper" not in bib

    md = md_path.read_text()
    assert "Jane Doe" in md
    assert "Alex Helper (acknowledged)" in md
    assert "[0000-0002-1825-0097](https://orcid.org/0000-0002-1825-0097)" in md
    assert "[@janedoe](https://github.com/janedoe)" in md
    assert "Curated thermal tolerance studies." in md


def test_schema_rejects_bad_contributor(tmp_path):
    bad = _sample_data()
    bad["contributors"].append({"given_name": "MissingFamily"})
    with pytest.raises(jsonschema.ValidationError):
        save_contributors(bad, tmp_path / "contributors.yaml")


def test_cli_add_sync_and_cite(tmp_path):
    contrib_path = tmp_path / "contributors.yaml"
    cff_path = tmp_path / "CITATION.cff"
    md_path = tmp_path / "CONTRIBUTORS.md"
    save_contributors(_sample_data(), contrib_path)

    res = runner.invoke(
        app,
        [
            "add-contributor",
            "--given-name", "Sam",
            "--family-name", "Lee",
            "--affiliation", "Shellfish Institute",
            "--role", "Formal analysis",
            "--contributors-file", str(contrib_path),
            "--no-sync",
        ],
    )
    assert res.exit_code == 0, res.output
    assert "added Sam Lee" in res.output

    res_sync = runner.invoke(
        app,
        [
            "sync-attribution",
            "--contributors-file", str(contrib_path),
            "--citation-out", str(cff_path),
            "--markdown-out", str(md_path),
        ],
    )
    assert res_sync.exit_code == 0, res_sync.output
    assert cff_path.exists() and md_path.exists()

    for fmt, needle in [("text", "Roberts, S. & Lee, S."), ("bibtex", "@software{arec_2026,"), ("cff", "cff-version: 1.2.0")]:
        res_cite = runner.invoke(app, ["cite", "--format", fmt, "--contributors-file", str(contrib_path)])
        assert res_cite.exit_code == 0, res_cite.output
        assert needle in res_cite.output


def test_committed_attribution_files_are_valid_and_in_sync():
    data = load_contributors(default_contributors_path())
    assert default_citation_path().read_text() == render_citation_cff(data)
    assert default_markdown_path().read_text() == render_contributors_markdown(data)


def test_website_build_generates_citation_page(tmp_path):
    from scripts import build_site

    out = build_site.build(tmp_path / "site_src")
    citation_qmd = Path(out) / "citation.qmd"
    assert citation_qmd.exists()
    content = citation_qmd.read_text()
    assert 'title: "Citation & Contributors"' in content
    assert "Steven Roberts" in content
    assert "[Studies](studies.qmd)" in content
    assert "## Primary Study Citations" in content
    assert "PRJNA" in content

    index_qmd = (Path(out) / "index.qmd").read_text()
    assert "## How to cite AREC" in index_qmd
    assert "[Citation & contributors](citation.qmd)" in index_qmd

    quarto_yml = yaml.safe_load((Path(out) / "_quarto.yml").read_text())
    footer_right = quarto_yml["website"]["page-footer"]["right"]
    assert any(item.get("href") == "citation.qmd" for item in footer_right)
