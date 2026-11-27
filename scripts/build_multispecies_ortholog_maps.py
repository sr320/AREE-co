"""Build versioned multi-species bivalve ortholog mapping tables from OrthoFinder outputs."""

import argparse
from pathlib import Path

from aree.harmonize.orthologs import DEFAULT_RELEASE, build_all_multispecies_maps


DEFAULT_ORTHOFINDER_DIR = Path(
    "/Volumes/ugreen-2TB/aree/orthogroups/orthofinder3.1.5_2026-11-26/proteomes/OrthoFinder/Results_Nov26"
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthofinder-dir", type=Path, default=DEFAULT_ORTHOFINDER_DIR)
    parser.add_argument("--output-dir", type=Path, default=Path("data/mappings/orthologs"))
    parser.add_argument("--release", default=DEFAULT_RELEASE)
    args = parser.parse_args()

    summaries = build_all_multispecies_maps(
        results_dir=args.orthofinder_dir,
        output_dir=args.output_dir,
        release=args.release,
    )
    for item in summaries:
        print(
            "{species_code} ({species}): {total_mapped_genes} genes "
            "[1:1={one_to_one}, M:1={many_to_one}, 1:M={one_to_many}, M:N={many_to_many_inferred}]".format(**item)
        )


if __name__ == "__main__":
    main()
