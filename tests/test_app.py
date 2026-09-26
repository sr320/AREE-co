from pathlib import Path

import pytest

pytest.importorskip("streamlit.testing.v1", reason="app smoke test needs streamlit >= 1.28 (the [app] extra)")
from streamlit.testing.v1 import AppTest  # noqa: E402


APP = Path(__file__).resolve().parents[1] / "app" / "main.py"


def test_app_renders_every_tab_without_errors():
    app = AppTest.from_file(str(APP), default_timeout=60).run()
    assert not app.exception
    assert [tab.label for tab in app.tabs] == ["Studies", "Evidence", "Candidates", "Pipeline Status"]
    assert len(app.dataframe) >= 3


def test_evidence_search_and_card_selection(tmp_path):
    app = AppTest.from_file(str(APP), default_timeout=60).run()
    app.text_input[0].input("LOC105317001").run()
    assert not app.exception
    evidence = app.dataframe[1].value
    assert len(evidence) and evidence["feature_id_standardized"].eq("NCBI:LOC105317001").all()
    app.selectbox[0].select("NCBI:LOC105317001").run()
    assert not app.exception
    assert any("Evidence Card: NCBI:LOC105317001" in block.value for block in app.markdown)
