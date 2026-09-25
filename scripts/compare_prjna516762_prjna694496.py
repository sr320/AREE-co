"""Create a non-pooled, cross-context comparison of the two real RNA-seq studies."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def compare(heat_path, selection_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    heat = pd.read_csv(heat_path, sep="\t")
    selection = pd.read_csv(selection_path, sep="\t")
    required = {"feature_id_standardized", "log2FoldChange", "lfcSE", "stat", "padj", "gene_symbol"}
    for name, table in (("heat", heat), ("selection", selection)):
        if not required.issubset(table.columns):
            raise ValueError("{} results lack required columns".format(name))
        if table.feature_id_standardized.duplicated().any():
            raise ValueError("{} results contain duplicate identifiers".format(name))
    keep = ["feature_id_standardized", "gene_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]
    common = heat[keep].merge(selection[keep], on="feature_id_standardized", how="inner",
                              validate="one_to_one", suffixes=("_acute_heat", "_selection"))
    common["significant_acute_heat"] = common.padj_acute_heat < 0.05
    common["significant_selection"] = common.padj_selection < 0.05
    common["same_direction"] = common.log2FoldChange_acute_heat * common.log2FoldChange_selection > 0
    common["cross_context_class"] = np.select(
        [
            common.significant_acute_heat & common.significant_selection & common.same_direction,
            common.significant_acute_heat & common.significant_selection & ~common.same_direction,
            common.significant_acute_heat & ~common.significant_selection,
            ~common.significant_acute_heat & common.significant_selection,
        ],
        ["significant_both_concordant", "significant_both_opposite",
         "acute_heat_only", "selection_only"],
        default="not_significant_either",
    )
    common["cross_context_score"] = np.minimum(common.stat_acute_heat.abs(), common.stat_selection.abs())
    common = common.sort_values(
        ["cross_context_class", "cross_context_score"], ascending=[True, False])
    common.to_csv(output_dir / "gene_level_cross_context_comparison.tsv", sep="\t", index=False)
    candidates = common[common.cross_context_class == "significant_both_concordant"].copy()
    candidates.sort_values("cross_context_score", ascending=False).to_csv(
        output_dir / "cross_context_concordant_candidates.tsv", sep="\t", index=False)

    valid = common.dropna(subset=["log2FoldChange_acute_heat", "log2FoldChange_selection"])
    pearson = valid.log2FoldChange_acute_heat.corr(valid.log2FoldChange_selection, method="pearson")
    spearman = valid.log2FoldChange_acute_heat.rank().corr(
        valid.log2FoldChange_selection.rank(), method="pearson")
    counts = common.cross_context_class.value_counts()
    lines = [
        "# Cross-context comparison: PRJNA516762 and PRJNA694496", "",
        "This comparison contrasts acute heat response in juvenile gill with constitutive expression associated with an artificially selected thermotolerant larval population. It ranks concordant signals but does not pool effect sizes because the biological contrasts, tissues, and life stages differ.", "",
        "- Common retained current-reference genes: {:,}.".format(len(common)),
        "- Pearson effect correlation: {:.3f}.".format(pearson),
        "- Spearman effect correlation: {:.3f}.".format(spearman),
        "- Significant in both with the same direction: {:,}.".format(
            counts.get("significant_both_concordant", 0)),
        "- Significant in both with opposite directions: {:,}.".format(
            counts.get("significant_both_opposite", 0)),
        "- Significant only after acute heat shock: {:,}.".format(counts.get("acute_heat_only", 0)),
        "- Significant only in the selected population comparison: {:,}.".format(
            counts.get("selection_only", 0)), "",
        "`cross_context_score` is the smaller absolute DESeq2 Wald statistic across the two studies. It ranks the strength of the weaker signal and is not a pooled effect or validation probability.", "",
    ]
    (output_dir / "cross_context_report.md").write_text("\n".join(lines))
    return common


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--heat", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    compare(args.heat, args.selection, args.output_dir)
