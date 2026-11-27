# PRJNA735889: tank-aware juvenile ocean acidification reanalysis

## Design and methods

Fifty whole-animal juvenile Pacific oyster (*Crassostrea gigas*) RNA-seq libraries from 10 tanks (5 control tanks at pHT 7.4–7.8, 5 acidified tanks at pHT 6.5–6.8; 5 oysters per tank; library `13Ebis` excluded as a technical resequencing of `13E`) were quantified with Salmon 2.3.4 against `GCF_963853765.1` / `RS_2024_06`.

Because oysters are nested within tanks (`one_tank_per_ph_level`) and male gametogenesis is unbalanced between bands (active in 11/25 control oysters across all 5 control tanks vs 2/25 acidified oysters in `t7` and `t14`), the primary inference model uses `tximport` (`countsFromAbundance = "lengthScaledTPM"`), DESeq2 median-of-ratios effective library sizes, `limma::voom` precision weights, and `variancePartition::dream` mixed-effects modeling with a random tank intercept and Satterthwaite degrees of freedom:

`~ gametogenesis + condition + (1 | tank)`

Genes passing the count filter (`>= 10` raw counts in `>= 3` libraries; 28,811 genes) were tested with `eBayes(..., robust = TRUE)`.

## Primary results (`dream`)

- Tested genes: 28,811 (non-missing adjusted p-values: 28,811).
- Significant genes at FDR < 0.05: **264** (**15** higher in acidified, **249** lower in acidified).
- Median `log2FoldChange` standard error: **0.2801** (versus 0.2164 under the prior unmodeled oyster-level DESeq2 run).
- Satterthwaite degrees of freedom (`conditionacidified`): median **49.03**, minimum **5.19** (for genes with substantial tank-to-tank variance), maximum **49.44**.
- SE consistency check against `dream` moderated *t*-statistics: max `|logFC / lfcSE - t| = 1.776e-15`.

## Model sensitivity comparison

| Model | Genes tested | FDR < 0.05 | Higher | Lower | Median SE | Median df | Pearson r vs primary | Overlap with primary | Direction agreement |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `dream_gametogenesis_plus_condition_tank_re` | 28,811 | 264 | 15 | 249 | 0.2801 | 49.03 | 1.0000 | 264 | 1.0 |
| `limma_voom_duplicateCorrelation` | 28,811 | 1,271 | 169 | 1,102 | 0.2616 | 52.41 | 0.9993 | 264 | 1.0 |
| `tank_pseudobulk_deseq2_gametogenesis_frac` | 28,657 | 134 | 78 | 56 | 0.4200 | 7.00 | 0.4654 | 28 | 1.0 |
| `tank_pseudobulk_deseq2_unadjusted` | 28,083 | 2,830 | 1,810 | 1,020 | 0.2418 | 8.00 | 0.7801 | 221 | 1.0 |
| `oyster_deseq2_unmodeled_tank_prior` | 28,248 | 4,725 | 2,735 | 1,990 | 0.2164 | 47.00 | 0.7799 | 246 | 1.0 |

Modeling the tank random effect removes 4,479 nominal calls from the earlier unmodeled oyster-level DESeq2 model while preserving 100% directional agreement on overlapping significant genes (`246/264`). Tank-level pseudobulk DESeq2 with `gametogenesis_frac` suffers strong collinearity between tank-level gametogenesis fraction and pH band (`median SE = 0.4200` on 7 residual df), whereas `dream` estimates the gametogenesis effect within tanks across all 50 oysters while shrinking condition degrees of freedom per gene via Satterthwaite approximation when tank variance is present.

## Remaining design limitations

- Five intermediate tanks in the pHT 7.3–6.9 tipping-point window are excluded, so effects contrast two pH bands rather than estimating a continuous pH slope.
- No tank was held at ambient coastal pH (~8.0); the reference band is pHT 7.4–7.8.
- Individual pH levels are unreplicated (`one_tank_per_ph_level`, with `t9` and `t15` both at pHT 6.5).
- `gametogenesis` is a transcriptional proxy derived from unadjusted PC1 (`scripts/score_prjna735889_gametogenesis.py`), not a histological sex/maturation assay.
