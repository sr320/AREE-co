import gzip
import re
from urllib.parse import unquote


def parse_attributes(text):
    """Parse a GFF3 attribute column (``key=value;key=value``) into a dict with URL-decoded values."""
    attributes = {}
    for item in text.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = unquote(value)
    return attributes


def gene_annotations(path, annotation_release):
    """Gene and pseudogene annotations from a gzipped RefSeq GFF, keyed by ``NCBI:GeneID:<id>``."""
    annotations = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "pseudogene"}:
                continue
            attributes = parse_attributes(fields[8])
            match = re.search(r"(?:^|,)GeneID:(\d+)(?:,|$)", attributes.get("Dbxref", ""))
            if not match:
                continue
            standardized = "NCBI:GeneID:" + match.group(1)
            annotations[standardized] = {
                "feature_id_standardized": standardized,
                "gene_symbol": attributes.get("gene", attributes.get("Name", "")),
                "description": attributes.get("description", ""),
                "gene_biotype": attributes.get("gene_biotype", ""),
                "annotation_release": annotation_release,
            }
    return annotations
