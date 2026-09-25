"""Report the six-library reanalysis and descriptive overlap with PRJNA516762."""

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def summarize(root, output, prior_evidence):
    output.mkdir(parents=True, exist_ok=True)
    de_dir = root / "deseq2"
    result_path = de_dir / "selected_vs_control_deseq2_all_genes.tsv"
    effects = pd.read_csv(result_path, sep="\t")
    ranking = pd.read_csv(de_dir / "selected_vs_control_deseq2_apeglm_ranking.tsv", sep="\t")
    qc = pd.read_csv(root / "salmon/salmon_qc_summary.tsv", sep="\t")
    if len(qc) != 6 or qc["sample"].nunique() != 6:
        raise ValueError("Six unique quantified libraries are required")
    if effects.feature_id_standardized.duplicated().any():
        raise ValueError("Duplicate gene IDs in DESeq2 output")
    annotation = pd.read_csv(root / "reference/GCF_963853765.1_tx2gene.tsv", sep="\t")
    annotation = annotation.groupby("gene_id", as_index=False).agg(
        source_transcript_ids=("transcript_id", lambda x: ";".join(sorted(set(x)))),
        gene_symbol=("gene_symbol", lambda x: ";".join(sorted(set(x.dropna())))))
    annotated = effects.merge(annotation, left_on="feature_id_standardized", right_on="gene_id", how="left", validate="one_to_one")
    annotated = annotated.merge(ranking[["feature_id_standardized", "log2FoldChange"]].rename(
        columns={"log2FoldChange": "apeglm_log2FoldChange_ranking_only"}), on="feature_id_standardized", validate="one_to_one")
    annotated.to_csv(output / "gene_results_annotated.tsv", sep="\t", index=False)
    significant = annotated[annotated.padj < 0.05].copy()
    significant["absolute_shrunken_effect"] = significant.apeglm_log2FoldChange_ranking_only.abs()
    top = significant.sort_values("absolute_shrunken_effect", ascending=False).head(25)
    top.to_csv(output / "top_25_genes.tsv", sep="\t", index=False)
    prior = pd.read_csv(prior_evidence, sep="\t")
    prior = prior[(prior.study_id == "CGIG_HEAT_RNASEQ_PRJNA516762") &
                  (prior.mapping_confidence != "unresolved")].copy()
    overlap = prior[["feature_id_original", "feature_id_standardized", "effect_size", "adjusted_p_value", "mapping_confidence"]].merge(
        annotated, on="feature_id_standardized", how="inner", validate="one_to_one")
    overlap["same_molecular_direction"] = overlap.effect_size * overlap.log2FoldChange > 0
    overlap.to_csv(output / "cross_study_descriptive_overlap.tsv", sep="\t", index=False)
    for name in ["analysis_summary.tsv", "sample_pca.csv", "sample_pca_variance.csv",
                 "sample_vst_correlations.csv", "sample_deseq2_qc.tsv", "sample_pca.png",
                 "selected_vs_control_MA.png", "sample_correlations.png", "sessionInfo.txt"]:
        shutil.copy2(de_dir / name, output / name)
    shutil.copy2(root / "salmon/salmon_qc_summary.tsv", output / "salmon_qc_summary.tsv")
    provenance = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in
                  [result_path, de_dir / "selected_vs_control_deseq2_apeglm_ranking.tsv",
                   root / "reference/GCF_963853765.1_tx2gene.tsv", prior_evidence]}
    (output / "input_sha256.json").write_text(json.dumps(provenance, indent=2) + "\n")
    lines = ["# PRJNA694496: selected-versus-control reanalysis", "",
             "Generated " + datetime.now(timezone.utc).isoformat() + ".", "",
             "Exploratory study-level results; QC must be reviewed before evidence harmonization.", "",
             "## Design and methods", "",
             "Three selected and three control whole-larva libraries; positive effects indicate higher constitutive expression in the selected population. Reference: GCF_963853765.1, RS_2024_06. Salmon 1.10.3 with decoys, sequence and GC bias correction; tximport gene-level summarization; DESeq2 ~ condition. Genes require at least 10 counts in at least three libraries. BH-adjusted p < 0.05 defines significance. Unshrunk effects and standard errors are retained; apeglm effects are for ranking only.", "",
             "## Results", "",
             f"- Genes retained after count filtering: {len(effects):,}.",
             f"- Genes with nonmissing raw p-values: {effects.pvalue.notna().sum():,}.",
             f"- Genes with nonmissing adjusted p-values: {effects.padj.notna().sum():,}.",
             f"- Significant genes: {len(significant):,}; {(significant.log2FoldChange > 0).sum():,} higher and {(significant.log2FoldChange < 0).sum():,} lower in selected larvae.",
             f"- Mapping rates: {qc.mapping_rate_percent.min():.2f}%–{qc.mapping_rate_percent.max():.2f}%.", "",
             "![Sample PCA](sample_pca.png)", "", "![MA plot](selected_vs_control_MA.png)", "",
             "## Descriptive comparison with the earlier heat-response study", "",
             f"Of {len(prior)} mapped, published significant heat-response genes in PRJNA516762, {len(overlap)} occur among the genes retained here; {(overlap.padj < 0.05).sum()} are significant here. The overlap table retains both effect directions and the earlier study's inferred mapping confidence.", "",
             "No pooled effect or enrichment test is calculated. PRJNA516762 provides a significant-only gene list without standard errors, and its acute heat response in juvenile gill differs from constitutive expression in selected larvae. Direction agreement is descriptive and does not establish a common protective mechanism.", "",
             "## Interpretation limits", "",
             "Selection is confounded with population background. Survival was assessed in related later-stage cohorts, not the sequenced larval pools. These results are associations, not validated predictors of individual survival. Mapping-rate variation, sample correlation/PCA, and Cook's distances require joint review; no samples were removed automatically. Missing p-values and adjusted p-values remain missing.", "",
             "## Files", "",
             "- `gene_results_annotated.tsv`: all retained genes, effects, uncertainty, transcript provenance, and ranking-only shrunken effects.",
             "- `top_25_genes.tsv`: significant genes ordered by absolute shrunken effect.",
             "- `cross_study_descriptive_overlap.tsv`: mapped overlap, without statistical pooling.",
             "- QC tables, figures, software session information, and input SHA-256 checksums are alongside this report.", ""]
    (output / "analysis_report.md").write_text("\n".join(lines))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--prior-evidence", type=Path, default=Path("data/harmonized/evidence.tsv"))
    args = parser.parse_args()
    summarize(args.analysis_root, args.output, args.prior_evidence)
