from pathlib import Path

import pytest

pytest.importorskip("streamlit.testing.v1", reason="app smoke test needs streamlit >= 1.28 (the [app] extra)")
from streamlit.testing.v1 import AppTest  # noqa: E402


APP = Path(__file__).resolve().parents[1] / "app" / "main.py"


def test_app_renders_real_mode_by_default():
    app = AppTest.from_file(str(APP), default_timeout=60).run()
    assert not app.exception
    assert [tab.label for tab in app.tabs] == ["Studies", "Evidence", "Findings & Cards", "Pipeline Status"]
    assert len(app.dataframe) >= 4
    studies = app.dataframe[0].value
    assert "CGIG_HEAT_RNASEQ_001" not in set(studies["study_id"])
    assert "CGIG_HEAT_RNASEQ_PRJNA516762" in set(studies["study_id"])
    assert "score" not in app.dataframe[1].value.columns
    assert "score" not in app.dataframe[3].value.columns
    assert any("not validated biomarkers" in msg.value for msg in app.info)


def test_real_mode_search_and_card_selection():
    app = AppTest.from_file(str(APP), default_timeout=60).run()
    app.text_input[0].input("105327837").run()
    assert not app.exception
    evidence = app.dataframe[1].value
    assert len(evidence) and evidence["feature_id_standardized"].eq("NCBI:GeneID:105327837").all()
    app.selectbox[1].select("NCBI:GeneID:105327837").run()
    assert not app.exception
    assert any("Real Evidence Card: NCBI:GeneID:105327837" in block.value for block in app.markdown)
    assert any("nominal only" in block.value for block in app.markdown)


def test_synthetic_demo_mode_search_and_card_selection():
    app = AppTest.from_file(str(APP), default_timeout=60).run()
    app.radio[0].set_value("Synthetic demo").run()
    assert not app.exception
    assert any("simulated demo data" in msg.value for msg in app.warning)
    demo_studies = app.dataframe[0].value
    assert "CGIG_HEAT_RNASEQ_001" in set(demo_studies["study_id"])
    assert "CGIG_HEAT_RNASEQ_PRJNA516762" not in set(demo_studies["study_id"])

    app.text_input[0].input("LOC105317001").run()
    assert not app.exception
    evidence = app.dataframe[1].value
    assert len(evidence) and evidence["feature_id_standardized"].eq("NCBI:LOC105317001").all()
    app.selectbox[0].select("NCBI:LOC105317001").run()
    assert not app.exception
    assert any("Evidence Card: NCBI:LOC105317001" in block.value for block in app.markdown)


@pytest.mark.parametrize("query", ["[", "(", ".*"])
def test_search_treats_user_text_as_literal(query):
    app = AppTest.from_file(str(APP), default_timeout=60).run()
    app.text_input[0].input(query).run()
    assert not app.exception
    real_hits = app.dataframe[1].value
    if query == ".*":
        assert real_hits.empty
    else:
        assert not real_hits.empty
        cols = ["feature_id_standardized", "description", "study_id", "sample_comparison", "quality_flags"]
        assert real_hits[cols].astype(str).apply(
            lambda col: col.str.contains(query, case=False, na=False, regex=False)
        ).any(axis=1).all()

    app.radio[0].set_value("Synthetic demo").run()
    app.text_input[0].input(query).run()
    assert not app.exception
    assert app.dataframe[1].value.empty
