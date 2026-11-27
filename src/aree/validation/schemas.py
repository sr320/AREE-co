import json
from pathlib import Path

import jsonschema
import yaml

from aree.paths import package_path, root_path


CONTROLLED_VOCABULARIES = (
    ("assay_type", "assays.yaml", "assays"),
    ("phenotype", "phenotypes.yaml", "phenotypes"),
    ("stressor_class", "stressors.yaml", "stressors"),
)


def load_schema(name):
    with package_path("schemas", name).open() as handle:
        return json.load(handle)


def load_yaml(path):
    with Path(path).open() as handle:
        return yaml.safe_load(handle)


def load_controlled_vocabularies():
    vocab_dir = root_path("registry", "controlled_vocabularies")
    if not vocab_dir.is_dir():
        return {}
    vocabularies = {}
    for field, filename, key in CONTROLLED_VOCABULARIES:
        path = vocab_dir / filename
        if path.is_file():
            vocabularies[field] = set(load_yaml(path).get(key, []))
    return vocabularies


def validate_study_file(path, check_filename=True):
    path = Path(path)
    data = load_yaml(path)
    jsonschema.Draft7Validator(load_schema("study.schema.json")).validate(data)
    if check_filename and path.stem != data["study_id"]:
        raise ValueError(
            "Study filename {!r} does not match study_id {!r}".format(path.name, data["study_id"])
        )
    for field, allowed in load_controlled_vocabularies().items():
        if data[field] not in allowed:
            raise ValueError(
                "Invalid {} {!r} in {}; allowed values: {}".format(
                    field, data[field], data["study_id"], ", ".join(sorted(allowed))
                )
            )
    return data


def validate_evidence_records(records):
    validator = jsonschema.Draft7Validator(load_schema("evidence.schema.json"))
    for record in records:
        validator.validate(record)
    return True

