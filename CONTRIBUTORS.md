# Contributors & Citation

This file is generated from [`contributors.yaml`](contributors.yaml). To add or update a contributor, see [Adding a contributor](#adding-a-contributor) below.

## How to Cite

If you use AREC, please cite this repository and the underlying public datasets.

**Plain text**

> Roberts, S. (2026). Aquaculture Resilience Evidence Companion (Version 0.1.0) [Computer software]. https://github.com/sr320/AREE-co

**BibTeX**

```bibtex
@software{arec_2026,
  author  = {Roberts, Steven},
  title   = {{Aquaculture Resilience Evidence Companion}},
  version = {0.1.0},
  year    = {2026},
  url     = {https://github.com/sr320/AREE-co},
}
```

Machine-readable citation metadata is also provided in `CITATION.cff` (used automatically by GitHub's **Cite this repository** button). When reporting results derived from AREC, please also cite the primary public studies listed in `registry/study_registry.csv`.

## Contributors

| Name | Affiliation | Roles | ORCID | GitHub |
|---|---|---|---|---|
| Steven Roberts | University of Washington | Conceptualization, Methodology, Software, Data curation, Project administration |  | [@sr320](https://github.com/sr320) |

## Acknowledgments

- AI-assisted software development: Claude (Anthropic), Antigravity / Gemini (Google), and GitHub Copilot were used to assist with code generation, refactoring, testing, documentation, and workflow/website scaffolding. All code, statistical analyses, and scientific interpretations were directed, reviewed, and verified by human contributors, who take full responsibility for the repository and its outputs.

## Adding a Contributor

Add a contributor from the command line (this updates `contributors.yaml`, `CITATION.cff`, and `CONTRIBUTORS.md` together):

```bash
aree add-contributor --given-name "Jane" --family-name "Doe" \
  --affiliation "University of Example" \
  --orcid "0000-0000-0000-0000" \
  --github "janedoe" \
  --role "Data curation" --role "Software"
```

Or edit `contributors.yaml` directly and run:

```bash
aree sync-attribution
```

Standard [CRediT](https://credit.niso.org/) roles are recommended (e.g., `Conceptualization`, `Data curation`, `Formal analysis`, `Methodology`, `Software`, `Validation`, `Visualization`, `Writing - review & editing`), though custom role descriptions are also accepted.
