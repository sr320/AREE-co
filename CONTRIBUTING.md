# Contributing

AREE welcomes code, documentation, workflow modules, and curated public-study metadata.

## Checks

Run these before opening a pull request; CI runs the same checks.

```bash
ruff check .
pytest
```

CI also reruns the demo pipeline and fails if any committed file under `data/demo/` or `reports/` changes, so commit regenerated outputs together with the code that changes them.

## Add a Study

1. Copy `registry/studies/CGIG_HEAT_RNASEQ_001.yaml`.
2. Replace identifiers, accessions, study context, assay metadata, and caveats.
3. Preserve original treatment labels and map them separately to controlled terms.
4. Run `aree validate-study registry/studies/YOUR_STUDY.yaml`.
5. Run `aree register-study registry/studies/YOUR_STUDY.yaml`.

## Evidence Claims

AREE reports associations and evidence convergence. Do not label a candidate as validated unless validation evidence exists outside this resource and is cited.

## Data Governance

Only register public or appropriately shareable datasets. Preserve provenance links, raw/processed availability, limitations, and manual curation decisions.

## Attribution and Adding Contributors

Contributor and citation metadata live in `contributors.yaml` (validated by `src/aree/schemas/contributors.schema.json`). Both `CITATION.cff` and `CONTRIBUTORS.md`—as well as the website's **Citation & contributors** page—are generated from it.

Add yourself or another contributor with:

```bash
aree add-contributor \
  --given-name "Jane" \
  --family-name "Doe" \
  --affiliation "University of Example" \
  --orcid "0000-0000-0000-0000" \
  --github "janedoe" \
  --role "Data curation" --role "Software"
```

This appends the entry to `contributors.yaml` and immediately regenerates `CITATION.cff` and `CONTRIBUTORS.md`. Use `--update` to update an existing entry, or `--acknowledge-only` to list someone in `CONTRIBUTORS.md` without adding them to the citation author list.

If you edit `contributors.yaml` by hand, regenerate the derived files with:

```bash
aree sync-attribution
```

Print the current citation in plain text, BibTeX, or CFF format with `aree cite --format text|bibtex|cff`.


