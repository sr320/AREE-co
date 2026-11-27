"""Single-source-of-truth contributor attribution and citation management.

Reads ``contributors.yaml`` at the repository root, validates it against
``contributors.schema.json``, and renders ``CITATION.cff``, ``CONTRIBUTORS.md``,
BibTeX, plain-text citations, and website attribution pages.
"""

import datetime
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import jsonschema
import yaml

from aree.paths import root_path
from aree.validation.schemas import load_schema


CREDIT_ROLES = (
    "Conceptualization",
    "Data curation",
    "Formal analysis",
    "Funding acquisition",
    "Investigation",
    "Methodology",
    "Project administration",
    "Resources",
    "Software",
    "Supervision",
    "Validation",
    "Visualization",
    "Writing - original draft",
    "Writing - review & editing",
)

ORCID_RE = re.compile(r"^(?:https?://orcid\.org/)?(\d{4}-\d{4}-\d{4}-\d{3}[\dX])$")


def default_contributors_path() -> Path:
    return root_path("contributors.yaml")


def default_citation_path() -> Path:
    return root_path("CITATION.cff")


def default_markdown_path() -> Path:
    return root_path("CONTRIBUTORS.md")


def normalize_orcid(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    match = ORCID_RE.match(text)
    if not match:
        raise ValueError(
            "Invalid ORCID {!r}; expected 0000-0000-0000-0000 or https://orcid.org/0000-0000-0000-0000".format(value)
        )
    return "https://orcid.org/" + match.group(1)


def normalize_github(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip().lstrip("@")
    return text or None


def _coerce_dates(data: Dict[str, Any]) -> Dict[str, Any]:
    """Ensure YAML-parsed date objects are ISO strings before schema validation."""
    if isinstance(data, dict):
        project = data.get("project")
        if isinstance(project, dict) and isinstance(project.get("date_released"), (datetime.date, datetime.datetime)):
            project["date_released"] = project["date_released"].strftime("%Y-%m-%d")
    return data


def validate_contributors_data(data: Dict[str, Any]) -> Dict[str, Any]:
    _coerce_dates(data)
    jsonschema.Draft7Validator(load_schema("contributors.schema.json")).validate(data)
    return data


def load_contributors(path: Optional[Path] = None) -> Dict[str, Any]:
    target = Path(path) if path is not None else default_contributors_path()
    with target.open() as handle:
        data = yaml.safe_load(handle)
    return validate_contributors_data(data)


def save_contributors(data: Dict[str, Any], path: Optional[Path] = None) -> Path:
    validated = validate_contributors_data(data)
    target = Path(path) if path is not None else default_contributors_path()
    target.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "# Single source of truth for AREC contributors and citation metadata.\n"
        "# Add contributors with `aree add-contributor` or edit this file directly,\n"
        "# then run `aree sync-attribution` to update CITATION.cff and CONTRIBUTORS.md.\n"
    )
    body = yaml.safe_dump(validated, sort_keys=False, allow_unicode=True)
    target.write_text(header + body)
    return target


def _contributor_match_index(contributors: Sequence[Dict[str, Any]], candidate: Dict[str, Any]) -> Optional[int]:
    cand_orcid = normalize_orcid(candidate.get("orcid"))
    cand_github = (normalize_github(candidate.get("github")) or "").lower()
    cand_email = (candidate.get("email") or "").strip().lower()
    cand_name = (
        candidate["given_name"].strip().lower(),
        candidate["family_name"].strip().lower(),
    )

    for idx, entry in enumerate(contributors):
        if cand_orcid and normalize_orcid(entry.get("orcid")) == cand_orcid:
            return idx
        if cand_github and (normalize_github(entry.get("github")) or "").lower() == cand_github:
            return idx
        if cand_email and (entry.get("email") or "").strip().lower() == cand_email:
            return idx
        entry_name = (
            entry["given_name"].strip().lower(),
            entry["family_name"].strip().lower(),
        )
        if entry_name == cand_name:
            return idx
    return None


def _clean_contributor(
    given_name: str,
    family_name: str,
    affiliation: Optional[str] = None,
    orcid: Optional[str] = None,
    github: Optional[str] = None,
    email: Optional[str] = None,
    roles: Optional[Sequence[str]] = None,
    include_in_citation: bool = True,
    notes: Optional[str] = None,
) -> Dict[str, Any]:
    entry: Dict[str, Any] = {
        "given_name": given_name.strip(),
        "family_name": family_name.strip(),
    }
    if affiliation and affiliation.strip():
        entry["affiliation"] = affiliation.strip()
    norm_orcid = normalize_orcid(orcid)
    if norm_orcid:
        entry["orcid"] = norm_orcid
    norm_github = normalize_github(github)
    if norm_github:
        entry["github"] = norm_github
    if email and email.strip():
        entry["email"] = email.strip()
    cleaned_roles = [r.strip() for r in (roles or []) if r and r.strip()]
    if cleaned_roles:
        # Preserve order while removing duplicates
        seen = set()
        unique_roles = []
        for role in cleaned_roles:
            if role not in seen:
                seen.add(role)
                unique_roles.append(role)
        entry["roles"] = unique_roles
    entry["include_in_citation"] = bool(include_in_citation)
    if notes and notes.strip():
        entry["notes"] = notes.strip()
    return entry


def add_contributor(
    given_name: str,
    family_name: str,
    affiliation: Optional[str] = None,
    orcid: Optional[str] = None,
    github: Optional[str] = None,
    email: Optional[str] = None,
    roles: Optional[Sequence[str]] = None,
    include_in_citation: bool = True,
    notes: Optional[str] = None,
    update: bool = False,
    path: Optional[Path] = None,
) -> Tuple[Path, str]:
    """Add or update a contributor in ``contributors.yaml``.

    Returns ``(contributors_path, action)`` where ``action`` is ``'added'`` or ``'updated'``.
    """
    data = load_contributors(path)
    candidate = _clean_contributor(
        given_name=given_name,
        family_name=family_name,
        affiliation=affiliation,
        orcid=orcid,
        github=github,
        email=email,
        roles=roles,
        include_in_citation=include_in_citation,
        notes=notes,
    )
    contributors = data.setdefault("contributors", [])
    match_idx = _contributor_match_index(contributors, candidate)

    if match_idx is not None:
        if not update:
            existing = contributors[match_idx]
            raise ValueError(
                "Contributor already exists: {} {} (pass --update to update their entry)".format(
                    existing["given_name"], existing["family_name"]
                )
            )
        merged = dict(contributors[match_idx])
        for key, value in candidate.items():
            if key == "roles":
                existing_roles = list(merged.get("roles") or [])
                for role in value:
                    if role not in existing_roles:
                        existing_roles.append(role)
                merged["roles"] = existing_roles
            else:
                merged[key] = value
        contributors[match_idx] = merged
        saved = save_contributors(data, path)
        return saved, "updated"

    contributors.append(candidate)
    saved = save_contributors(data, path)
    return saved, "added"


def cited_contributors(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    return [c for c in data.get("contributors", []) if c.get("include_in_citation", True)]


def render_citation_cff(data: Dict[str, Any]) -> str:
    project = data["project"]
    authors = []
    for person in cited_contributors(data):
        author: Dict[str, Any] = {
            "family-names": person["family_name"],
            "given-names": person["given_name"],
        }
        if person.get("affiliation"):
            author["affiliation"] = person["affiliation"]
        if person.get("orcid"):
            author["orcid"] = normalize_orcid(person["orcid"])
        if person.get("email"):
            author["email"] = person["email"]
        authors.append(author)

    if not authors:
        authors = [{"name": "AREC contributors"}]

    cff: Dict[str, Any] = {
        "cff-version": "1.2.0",
        "title": project["title"],
        "message": project["message"],
        "type": "software",
        "authors": authors,
        "version": project["version"],
        "date-released": project["date_released"],
        "repository-code": project["repository_code"],
    }
    if project.get("url"):
        cff["url"] = project["url"]
    if project.get("doi"):
        cff["doi"] = project["doi"]
    cff["license"] = project["license"]

    return yaml.safe_dump(cff, sort_keys=False, allow_unicode=True)


def _initials(given_name: str) -> str:
    parts = [p for p in re.split(r"[\s-]+", given_name.strip()) if p]
    return " ".join(p[0].upper() + "." for p in parts)


def render_plain_citation(data: Dict[str, Any]) -> str:
    project = data["project"]
    year = str(project["date_released"])[:4]
    people = cited_contributors(data)
    if people:
        formatted = ["{}, {}".format(p["family_name"], _initials(p["given_name"])) for p in people]
        if len(formatted) == 1:
            author_str = formatted[0]
        elif len(formatted) == 2:
            author_str = "{} & {}".format(formatted[0], formatted[1])
        else:
            author_str = ", ".join(formatted[:-1]) + ", & " + formatted[-1]
    else:
        author_str = "AREC contributors"

    location = "https://doi.org/{}".format(project["doi"]) if project.get("doi") else project["repository_code"]
    return "{} ({}). {} (Version {}) [Computer software]. {}".format(
        author_str, year, project["title"], project["version"], location
    )


def render_bibtex(data: Dict[str, Any]) -> str:
    project = data["project"]
    year = str(project["date_released"])[:4]
    abbrev = (project.get("abbreviation") or "arec").lower()
    key = "{}_{}".format(re.sub(r"[^a-z0-9]+", "", abbrev) or "arec", year)
    people = cited_contributors(data)
    if people:
        author_str = " and ".join("{}, {}".format(p["family_name"], p["given_name"]) for p in people)
    else:
        author_str = "{AREC contributors}"

    lines = [
        "@software{" + key + ",",
        "  author  = {" + author_str + "},",
        "  title   = {{" + project["title"] + "}},",
        "  version = {" + project["version"] + "},",
        "  year    = {" + year + "},",
        "  url     = {" + project["repository_code"] + "},",
    ]
    if project.get("doi"):
        lines.append("  doi     = {" + project["doi"] + "},")
    lines.append("}")
    return "\n".join(lines) + "\n"


def _md_escape(text: Optional[str]) -> str:
    if not text:
        return ""
    return str(text).replace("|", "\\|").replace("\n", " ")


def render_contributors_markdown(
    data: Dict[str, Any],
    include_front_matter: bool = False,
    primary_studies: Optional[Sequence[Dict[str, Any]]] = None,
) -> str:
    project = data["project"]
    contributors = data.get("contributors", [])
    lines: List[str] = []

    if include_front_matter:
        lines.extend(["---", 'title: "Citation & Contributors"', "---", ""])
    else:
        lines.extend(["# Contributors & Citation", ""])

    lines.append(
        "This file is generated from [`contributors.yaml`](contributors.yaml). "
        "To add or update a contributor, see [Adding a contributor](#adding-a-contributor) below."
        if not include_front_matter
        else "Contributor and citation metadata are maintained in `contributors.yaml` in the repository."
    )
    lines.append("")

    lines.extend(["## How to Cite", ""])
    lines.append(project["message"])
    lines.append("")
    lines.append("**Plain text**")
    lines.append("")
    lines.append("> " + render_plain_citation(data))
    lines.append("")
    lines.append("**BibTeX**")
    lines.append("")
    lines.append("```bibtex")
    lines.append(render_bibtex(data).rstrip())
    lines.append("```")
    lines.append("")
    studies_ref = "[Studies](studies.qmd)" if include_front_matter else "`registry/study_registry.csv`"
    lines.append(
        "Machine-readable citation metadata is also provided in `CITATION.cff` "
        "(used automatically by GitHub's **Cite this repository** button). "
        "When reporting results derived from AREC, please also cite the primary public studies "
        "listed in {}.".format(studies_ref)
    )
    lines.append("")

    lines.extend(["## Contributors", ""])
    if contributors:
        lines.append("| Name | Affiliation | Roles | ORCID | GitHub |")
        lines.append("|---|---|---|---|---|")
        for person in contributors:
            name = "{} {}".format(person["given_name"], person["family_name"])
            if not person.get("include_in_citation", True):
                name += " (acknowledged)"
            affiliation = _md_escape(person.get("affiliation") or "")
            roles = _md_escape(", ".join(person.get("roles") or []))
            orcid_raw = person.get("orcid")
            if orcid_raw:
                orcid_url = normalize_orcid(orcid_raw) or ""
                orcid_id = orcid_url.rsplit("/", 1)[-1]
                orcid_cell = "[{}]({})".format(orcid_id, orcid_url)
            else:
                orcid_cell = ""
            gh = normalize_github(person.get("github"))
            gh_cell = "[@{}](https://github.com/{})".format(gh, gh) if gh else ""
            lines.append("| {} | {} | {} | {} | {} |".format(
                _md_escape(name), affiliation, roles, orcid_cell, gh_cell
            ))
    else:
        lines.append("_No individual contributors listed yet._")
    lines.append("")

    noted = [p for p in contributors if p.get("notes")]
    if noted:
        lines.extend(["### Contribution Notes", ""])
        for person in noted:
            lines.append("- **{} {}**: {}".format(person["given_name"], person["family_name"], person["notes"]))
        lines.append("")

    acknowledgments = project.get("acknowledgments") or []
    if acknowledgments:
        lines.extend(["## Acknowledgments", ""])
        for item in acknowledgments:
            lines.append("- {}".format(item))
        lines.append("")

    if primary_studies:
        lines.extend([
            "## Primary Study Citations",
            "",
            "AREC synthesizes and reanalyzes public datasets generated by the research community. "
            "If you use findings from specific studies, please cite the corresponding primary publications or BioProjects:",
            "",
            "| BioProject | Study contrasts | Citation |",
            "|---|---|---|",
        ])
        for item in primary_studies:
            bp = _md_escape(item.get("bioproject") or "")
            study_links = ", ".join(
                "[{}](studies/{}.qmd)".format(sid, sid) for sid in item.get("study_ids", [])
            )
            citation_text = _md_escape(item.get("citation") or "")
            doi = (item.get("doi") or "").strip()
            if doi:
                citation_text = "{} ([doi:{}](https://doi.org/{}))".format(citation_text, doi, doi)
            lines.append("| {} | {} | {} |".format(bp, study_links, citation_text))
        lines.append("")

    lines.extend([
        "## Adding a Contributor",
        "",
        "Add a contributor from the command line (this updates `contributors.yaml`, `CITATION.cff`, and `CONTRIBUTORS.md` together):",
        "",
        "```bash",
        'aree add-contributor --given-name "Jane" --family-name "Doe" \\',
        '  --affiliation "University of Example" \\',
        '  --orcid "0000-0000-0000-0000" \\',
        '  --github "janedoe" \\',
        '  --role "Data curation" --role "Software"',
        "```",
        "",
        "Or edit `contributors.yaml` directly and run:",
        "",
        "```bash",
        "aree sync-attribution",
        "```",
        "",
        "Standard [CRediT](https://credit.niso.org/) roles are recommended "
        "(e.g., `Conceptualization`, `Data curation`, `Formal analysis`, `Methodology`, "
        "`Software`, `Validation`, `Visualization`, `Writing - review & editing`), "
        "though custom role descriptions are also accepted.",
        "",
    ])
    return "\n".join(lines)


def sync_attribution(
    contributors_path: Optional[Path] = None,
    citation_path: Optional[Path] = None,
    markdown_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Validate ``contributors.yaml`` and regenerate ``CITATION.cff`` and ``CONTRIBUTORS.md``."""
    data = load_contributors(contributors_path)
    cff_out = Path(citation_path) if citation_path is not None else default_citation_path()
    md_out = Path(markdown_path) if markdown_path is not None else default_markdown_path()
    cff_out.write_text(render_citation_cff(data))
    md_out.write_text(render_contributors_markdown(data))
    return {
        "n_contributors": len(data.get("contributors", [])),
        "n_cited_authors": len(cited_contributors(data)),
        "citation_path": str(cff_out),
        "markdown_path": str(md_out),
    }
