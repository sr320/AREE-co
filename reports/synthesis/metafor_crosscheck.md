# R `metafor` Cross-Check Report

Verified Python `aree.meta_analysis.random_effects` (DL, REML, and block-covariance shared-control GLS) against R `metafor` (`rma.uni` and `rma.mv`).

- Cases tested: 152 (12 synthetic edge/covariance cases + 140 real-study gene groups from `data/harmonized/`).
- Max absolute difference in `pooled_effect`: `9.846e-11`
- Max absolute difference in `pooled_standard_error`: `1.551e-10`
- Max absolute difference in `tau2`: `1.163e-10`
- Max absolute difference in `q`: `5.684e-14`
- Max absolute difference in `i2_percent`: `1.030e-07`
- Max relative difference in `p_value`: `2.289e-09`
- Max absolute difference vs `metafor::rma.mv(yi, V, random = ~ 1 | group, method = 'REML')`: `6.678e-07`

All checks passed within tolerance.
