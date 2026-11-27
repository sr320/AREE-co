import json
import sqlite3
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import pytest

from aree.export.database import (
    export_jsonld_catalog,
    export_parquet_bundle,
    export_schema_sql,
    export_sqlite,
)
from aree.reporting.evidence_cards import build_evidence_cards
from aree.reporting.forest_plot import render_forest_svg, write_forest_plot_svg

ROOT = Path(__file__).resolve().parents[1]


def test_forest_svg_valid_xml_escaping_and_determinism(tmp_path):
    rows = [
        {"label": "PRJNA1 <heat>", "effect": 1.45, "se": 0.25, "subtext": "FDR=1e-4"},
        {"label": "PRJNA2 & control", "effect": -0.82, "se": 0.30, "subtext": "FDR=0.02"},
        {"label": "Exploratory study", "effect": 0.65, "se": None, "subtext": "nominal"},
    ]
    pooled = {"label": "Pooled REML", "effect": 0.42, "ci_low": -0.1, "ci_high": 0.94, "subtext": "I2=62%"}
    svg1 = render_forest_svg(rows, pooled=pooled, title="Gene <HSP70> & family")
    svg2 = render_forest_svg(rows, pooled=pooled, title="Gene <HSP70> & family")
    assert svg1 == svg2
    elem = ET.fromstring(svg1)
    assert elem.tag.endswith("svg")
    assert "<polygon" in svg1
    assert "<circle" in svg1

    out = write_forest_plot_svg(tmp_path / "forest.svg", rows=rows, pooled=pooled, title="Gene <HSP70> & family")
    assert out.exists() and out.read_text() == svg1


def test_evidence_cards_plots_flag_opt_in(tmp_path):
    no_plots_dir = tmp_path / "cards_plain"
    build_evidence_cards(output_dir=no_plots_dir, include_plots=False)
    assert not (no_plots_dir / "plots").exists()

    with_plots_dir = tmp_path / "cards_plots"
    cards = build_evidence_cards(output_dir=with_plots_dir, include_plots=True)
    assert cards
    svgs = list((with_plots_dir / "plots").glob("*.svg"))
    assert len(svgs) == len(cards)
    assert "![Forest plot" in cards[0].read_text()


def test_database_exports_sqlite_parquet_and_jsonld(tmp_path):
    db_path = tmp_path / "aree_test.sqlite"
    sql_path = tmp_path / "schema.sql"
    export_schema_sql(sql_path)
    assert "CREATE TABLE IF NOT EXISTS harmonized_evidence" in sql_path.read_text()

    info = export_sqlite(db_path, root=ROOT)
    assert info["tables"]["studies"] >= 30
    assert info["tables"]["refseq_references"] == 8
    assert info["tables"]["harmonized_evidence"] > 100000
    assert info["tables"]["synthesis_strata"] == 5

    conn = sqlite3.connect(str(db_path))
    try:
        n_refs = conn.execute("SELECT COUNT(*) FROM refseq_references").fetchone()[0]
        assert n_refs == 8
        top = conn.execute(
            "SELECT feature_id_standardized FROM synthesis_meta WHERE stratum_id = 'acute_heat_response_gill' LIMIT 5"
        ).fetchall()
        assert len(top) == 5
    finally:
        conn.close()

    jsonld_path = tmp_path / "catalog.jsonld"
    jinfo = export_jsonld_catalog(jsonld_path, root=ROOT)
    assert jinfo["n_datasets"] >= 20
    payload = json.loads(jsonld_path.read_text())
    assert payload["@context"] == "https://schema.org"
    assert payload["@type"] == "DataCatalog"

    pyarrow = pytest.importorskip("pyarrow")
    assert pyarrow is not None
    pq_dir = tmp_path / "pq_bundle"
    pinfo = export_parquet_bundle(pq_dir, root=ROOT)
    pkg = json.loads(Path(pinfo["datapackage_path"]).read_text())
    assert pkg["profile"] == "tabular-data-package"
    names = {r["name"] for r in pkg["resources"]}
    assert {"studies", "refseq_references", "harmonized_evidence", "synthesis_strata", "synthesis_meta"}.issubset(names)
    df_refs = pd.read_parquet(pq_dir / "refseq_references.parquet")
    assert len(df_refs) == 8
