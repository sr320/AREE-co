from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def root_path(*parts):
    return REPO_ROOT.joinpath(*parts)

