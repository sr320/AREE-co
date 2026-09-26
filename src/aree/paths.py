import os
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parent
_SOURCE_CHECKOUT = PACKAGE_DIR.parents[1]


def _is_project(path):
    return (path / "registry" / "studies").is_dir()


def project_root():
    """The AREE project directory holding registry/, data/ and reports/.

    Resolved on each call: ``$AREE_ROOT`` if set; else the source checkout when the package is
    installed from one (editable installs); else the nearest directory at or above the current
    working directory that contains ``registry/studies``; else the current working directory.
    A regular (non-editable) install therefore works when run from inside a project.
    """
    configured = os.environ.get("AREE_ROOT")
    if configured:
        return Path(configured).resolve()
    if _is_project(_SOURCE_CHECKOUT):
        return _SOURCE_CHECKOUT
    cwd = Path.cwd().resolve()
    for candidate in (cwd, *cwd.parents):
        if _is_project(candidate):
            return candidate
    return cwd


def root_path(*parts):
    return project_root().joinpath(*parts)


def package_path(*parts):
    """Files shipped inside the package (e.g. JSON schemas), independent of the project directory."""
    return PACKAGE_DIR.joinpath(*parts)
