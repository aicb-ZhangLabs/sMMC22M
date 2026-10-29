# Full 23-category results

CONCH + Ridge; HD18 plus five additional **Xenium data only** categories.
Each macro weights categories equally. Acquisitions are not independent donors.
Overlapping Xenium categories are not pooled into the HD rows. No combined CI is claimed.

| Category | n | Data | Gene P | Coord P | KNN P | Gene S | Cell P | F1 | R² | RMSE |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Brain | 6 | Visium HD | 0.505 | 0.076 | 0.072 | 0.492 | 0.471 | 0.617 | 0.238 | 0.575 |
| Breast | 5 | Visium HD | 0.395 | 0.074 | 0.057 | 0.378 | 0.613 | 0.503 | 0.160 | 0.605 |
| Colon | 3 | Visium HD | 0.424 | 0.110 | 0.116 | 0.400 | 0.526 | 0.472 | 0.042 | 0.505 |
| Embryo | 3 | Visium HD | 0.320 | 0.041 | 0.060 | 0.273 | 0.519 | 0.418 | -0.303 | 0.496 |
| Head | 1 | Visium HD | 0.088 | -0.012 | -0.001 | 0.057 | 0.131 | 0.023 | -1.891 | 0.138 |
| Heart | 1 | Visium HD | 0.372 | 0.111 | 0.081 | 0.312 | 0.446 | 0.287 | 0.097 | 0.403 |
| Kidney | 4 | Visium HD | 0.368 | 0.061 | 0.045 | 0.347 | 0.545 | 0.374 | 0.127 | 0.461 |
| Liver | 2 | Visium HD | 0.477 | 0.008 | 0.023 | 0.466 | 0.627 | 0.781 | 0.226 | 0.622 |
| Lung | 6 | Visium HD | 0.432 | 0.009 | 0.058 | 0.370 | 0.497 | 0.479 | 0.140 | 0.494 |
| Lymph Node | 2 | Visium HD | 0.151 | 0.015 | 0.018 | 0.140 | 0.334 | 0.036 | -0.001 | 0.232 |
| Ovary | 4 | Visium HD | 0.393 | 0.021 | 0.036 | 0.385 | 0.550 | 0.517 | 0.130 | 0.524 |
| Pancreas | 3 | Visium HD | 0.403 | 0.101 | 0.092 | 0.368 | 0.724 | 0.528 | 0.117 | 0.544 |
| Plant (Arabidopsis) | 1 | Visium HD | 0.388 | -0.033 | -0.013 | 0.376 | 0.233 | 0.519 | -1.051 | 0.716 |
| Plant (Soybean) | 1 | Visium HD | 0.142 | 0.057 | 0.033 | 0.128 | 0.251 | 0.026 | -0.133 | 0.218 |
| Prostate | 1 | Visium HD | 0.260 | 0.033 | 0.065 | 0.236 | 0.292 | 0.090 | -0.174 | 0.257 |
| Small Intestine | 1 | Visium HD | 0.509 | 0.047 | 0.024 | 0.481 | 0.537 | 0.622 | 0.246 | 0.526 |
| Tonsil | 3 | Visium HD | 0.282 | 0.023 | 0.040 | 0.252 | 0.510 | 0.382 | 0.047 | 0.530 |
| Xenograft | 1 | Visium HD | 0.655 | 0.171 | 0.198 | 0.538 | 0.196 | 0.731 | 0.364 | 0.402 |
| Bone | 4 | Xenium data only | 0.417 | 0.058 | 0.089 | 0.380 | 0.524 | 0.480 | 0.157 | 0.538 |
| Bone Marrow | 1 | Xenium data only | 0.116 | 0.009 | 0.017 | 0.101 | 0.249 | 0.022 | -0.012 | 0.221 |
| Cervix | 1 | Xenium data only | 0.308 | 0.042 | 0.041 | 0.268 | 0.368 | 0.231 | 0.072 | 0.329 |
| Skin | 5 | Xenium data only | 0.402 | 0.060 | 0.102 | 0.347 | 0.476 | 0.464 | -0.460 | 0.515 |
| Whole Organism | 2 | Xenium data only | 0.274 | 0.038 | 0.056 | 0.218 | 0.477 | 0.249 | -0.517 | 0.393 |
| 18-category HD macro mean | 48 | Category-weighted mean | 0.365 | 0.051 | 0.056 | 0.333 | 0.445 | 0.411 | -0.090 | 0.458 |
| 5-category Xenium-only macro mean | 13 | Category-weighted mean | 0.303 | 0.042 | 0.061 | 0.263 | 0.419 | 0.289 | -0.152 | 0.399 |
| 23-category macro mean | 61 | Category-weighted mean | 0.351 | 0.049 | 0.057 | 0.318 | 0.439 | 0.385 | -0.103 | 0.445 |

Generated from the unrounded CSVs by `scripts/summarize_results.py`.
Coord P: quadratic-coordinate Ridge Gene Pearson; KNN P: spatial-KNN Gene Pearson.
See [protocols](../docs/protocols.md) and [all results](README.md) for scope and uncertainty.
