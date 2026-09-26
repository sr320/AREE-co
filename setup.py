# Kept so `pip install -e .` works with pip < 21.3, which cannot do editable installs from
# pyproject.toml alone (PEP 660). All metadata lives in pyproject.toml.
from setuptools import setup


setup()
