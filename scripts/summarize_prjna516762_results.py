"""Summarize the current-reference PRJNA516762 control/heat reanalysis."""

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def summarize(root, output, publication_evidence):
    output.mkdir(parents=True, exist_ok=True)
    de_dir = root / "deseq2"
    result_path = de_dir / "heat_vs_control_deseq2_all_genes.tsv"
    effects = pd.read_csv(result_path, sep="\t")
    ranking_path = de_dir / "heat_vs_control_deseq2_apeglm_ranking.tsv"
    ranking = pd.read_csv(ranking_path, sep="\t")
    qc = pd.read_csv(root / "salmon/salmon_qc_summary.tsv", sep="\t")
    pca = pd.read_csv(de_dir / "sample_pca.csv")
    pca_variance = pd.read_csv(de_dir / "sample_pca_variance.csv")
    correlations = pd.read_csv(de_dir / "sample_vst_correlations.csv", index_col=0)
    sample_qc = pd.read_csv(de_dir / "sample_deseq2_qc.tsv", sep="\t")
    if len(qc) != 6 or qc["sample"].nunique() != 6:
        raise ValueError("Six unique quantified libraries are required")
    if effects.feature_id_standardized.duplicated().any():
        raise ValueError("Duplicate gene IDs in DESeq2 output")

    tx2gene = pd.read_csv(
        root.parent / "CGIG_THERMOTOL_RNASEQ_PRJNA694496/reference/GCF_963853765.1_tx2gene.tsv",
        sep="\t",
    )
    annotation = tx2gene.groupby("gene_id", as_index=False).agg(
        source_transcript_ids=("transcript_id", lambda x: ";".join(sorted(set(x)))),
        gene_symbol=("gene_symbol", lambda x: ";".join(sorted(set(x.dropna())))),
    )
    annotated = effects.merge(annotation, left_on="feature_id_standardized", right_on="gene_id",
                              how="left", validate="one_to_one")
    annotated = annotated.merge(
        ranking[["feature_id_standardized", "log2FoldChange"]].rename(
            columns={"log2FoldChange": "apeglm_log2FoldChange_ranking_only"}),
        on="feature_id_standardized", validate="one_to_one")
    annotated.to_csv(output / "gene_results_annotated.tsv", sep="\t", index=False)
    significant = annotated[annotated.padj < 0.05].copy()
    significant["absolute_shrunken_effect"] = significant.apeglm_log2FoldChange_ranking_only.abs()
    significant.sort_values("absolute_shrunken_effect", ascending=False).head(25).to_csv(
        output / "top_25_genes.tsv", sep="\t", index=False)

    published = pd.read_csv(publication_evidence, sep="\t")
    published = published[(published.study_id == "CGIG_HEAT_RNASEQ_PRJNA516762") &
                          (published.mapping_confidence != "unresolved")].copy()
    reproduction = published[["feature_id_original", "feature_id_standardized", "effect_size",
                              "adjusted_p_value", "mapping_confidence"]].merge(
        annotated, on="feature_id_standardized", how="inner", validate="one_to_one",
        suffixes=("_published", "_reanalysis"))
    reproduction["same_molecular_direction"] = (
        reproduction.effect_size * reproduction.log2FoldChange > 0)
    reproduction.to_csv(output / "publication_reproduction.tsv", sep="\t", index=False)

    for name in ["analysis_summary.tsv", "sample_pca.csv", "sample_pca_variance.csv",
                 "sample_vst_correlations.csv", "sample_deseq2_qc.tsv", "sample_pca.png",
                 "heat_vs_control_MA.png", "sample_correlations.png", "sessionInfo.txt"]:
        shutil.copy2(de_dir / name, output / name)
    shutil.copy2(root / "salmon/salmon_qc_summary.tsv", output / "salmon_qc_summary.tsv")
    inputs = [result_path, ranking_path, tx2gene.attrs.get("source", root.parent /
              "CGIG_THERMOTOL_RNASEQ_PRJNA694496/reference/GCF_963853765.1_tx2gene.tsv"),
              publication_evidence]
    provenance = {str(Path(path)): hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in inputs}
    (output / "input_sha256.json").write_text(json.dumps(provenance, indent=2) + "\n")

    pc1 = 100 * pca_variance.loc[pca_variance.component == "PC1", "variance_fraction"].iloc[0]
    control_pc1 = pca.loc[pca.condition == "control", "PC1"]
    heat_pc1 = pca.loc[pca.condition == "heat", "PC1"]
    separated = control_pc1.max() < heat_pc1.min() or heat_pc1.max() < control_pc1.min()
    lines = [
        "# PRJNA516762: heat-shock-versus-control reanalysis", "",
        "Generated {}.".format(datetime.now(timezone.utc).isoformat()), "",
        "## Design and methods", "",
        "Three untreated-control and three 35 C, 2 h heat-shock juvenile-gill libraries were quantified against GCF_963853765.1 / RS_2024_06. Each library pools five oysters. Salmon used the version-matched decoy-aware index with sequence and GC bias correction; tximport summarized transcript estimates to genes; DESeq2 fit `~ condition`. Genes require at least 10 counts in at least three libraries. Positive effects indicate higher abundance after heat shock. BH-adjusted p < 0.05 defines significance; unshrunk effects and standard errors are retained for evidence, while apeglm effects are used only for ranking.", "",
        "## Results", "",
        "- Genes retained after count filtering: {:,}.".format(len(effects)),
        "- Genes with nonmissing adjusted p-values: {:,}.".format(effects.padj.notna().sum()),
        "- Significant genes: {:,}; {:,} higher and {:,} lower after heat shock.".format(
            len(significant), (significant.log2FoldChange > 0).sum(),
            (significant.log2FoldChange < 0).sum()),
        "- Mapping rates: {:.2f}%–{:.2f}%.".format(qc.mapping_rate_percent.min(),
                                                   qc.mapping_rate_percent.max()),
        "- PC1 explains {:.1f}% of variance; complete condition separation on PC1: {}.".format(
            pc1, "yes" if separated else "no"),
        "- Genes above the sample-level Cook's-distance threshold range from {} to {}.".format(
            int(sample_qc.genes_above_cooks_cutoff.min()),
            int(sample_qc.genes_above_cooks_cutoff.max())), "",
        "![Sample PCA](sample_pca.png)", "", "![MA plot](heat_vs_control_MA.png)", "",
        "## Publication comparison", "",
        "The publication reported 150 significant genes, of which {} have conservative current-GeneID mappings. The current-reference reanalysis retained {} of those mapped genes; {} are significant at FDR < 0.05 and {} have the same effect direction.".format(
            len(published), len(reproduction), int((reproduction.padj < 0.05).sum()),
            int(reproduction.same_molecular_direction.sum())), "",
        "Differences can reflect the newer reference annotation, transcript-level quantification, gene filtering, and incomplete mapping of legacy CGI identifiers. Publication agreement is a reproducibility check, not a criterion for deleting current-reference results.", "",
        "## Interpretation limits", "",
        "The libraries are biological replicates of pooled animals, so individual-level variability cannot be estimated. This contrast measures acute molecular heat response and does not directly measure survival or inherited thermotolerance. It should be compared with PRJNA694496 as cross-context evidence and should not be pooled as if both studies estimate the same biological contrast.", "",
        "## Files", "",
        "- `gene_results_annotated.tsv`: complete retained-gene results with standard errors and transcript provenance.",
        "- `publication_reproduction.tsv`: mapped comparison with the published significant-gene list.",
        "- `top_25_genes.tsv`: significant genes ranked by absolute apeglm effect.",
        "- QC tables, figures, software session information, and input SHA-256 checksums accompany this report.", "",
    ]
    (output / "analysis_report.md").write_text("\n".join(lines))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--publication-evidence", type=Path,
                        default=Path("data/harmonized/evidence.tsv"))
    args = parser.parse_args()
    summarize(args.analysis_root, args.output, args.publication_evidence)
