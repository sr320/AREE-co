import json
from pathlib import Path

import jsonschema
import yaml

from aree.paths import root_path


def load_schema(name):
    with root_path("schemas", name).open() as handle:
        return json.load(handle)


def load_yaml(path):
    with Path(path).open() as handle:
        return yaml.safe_load(handle)


def validate_study_file(path):
    data = load_yaml(path)
    jsonschema.Draft7Validator(load_schema("study.schema.json")).validate(data)
    return data


def validate_evidence_records(records):
    validator = jsonschema.Draft7Validator(load_schema("evidence.schema.json"))
    for record in records:
        validator.validate(record)
    return True

