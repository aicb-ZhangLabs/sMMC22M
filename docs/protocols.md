# Evaluation protocols and interpretation

## Modern within-acquisition benchmark

The original sampling retains at most 12,000 gene-positive cells with a centroid
inside the selected H&E image per acquisition, using a 20 × 20 spatial sampler.
Its frozen seed is 20260723 plus the acquisition's frozen index. Exact sampled
cell IDs and row links are retained in the data-side fixtures; the public fold
plans define the original targets and partition references.

Seeds 0–3 select high-x, high-y, low-x, and low-y edge holdouts. Each test region
contains approximately 20% of sampled cells. Training rows are separated by a
buffer of 5% of that coordinate axis's full span. These are contiguous spatial
holdouts, not random single-cell splits. Buffer rows are not training or test.

For each organ/species/fold group, intersect the **measured gene IDs** and rank
training-only log1p population variance with equal acquisition weights. Select
50 HVGs; exact gene-ID ascending order breaks ties. No library-size normalization
or zero-filled union panel is used. Skin's five Xenium acquisitions share only
47 measured IDs; the documented pre-evaluation amendment defines four exact-panel
subgroups, each with 50 training-selected targets. No Skin acquisition is discarded.

Predictors are refit per fold. Standardization, PCA and Ridge see training rows
only. Coordinate Ridge and KNN also use training expression only. Targets,
sampling, crop padding, and evaluator precision remain frozen for reproducibility.

## Cross-acquisition benchmark

HD holds out a complete acquisition within the same organ/species stratum.
There are **34 held-out acquisitions in 10 strata**; 50 HVGs are selected from
external training acquisitions, and image preprocessing/head fitting is confined
to those acquisitions. The primary test set is the full frozen sampled
acquisition. Macro Gene Pearson is **0.2603575231**.

The **matched-edge** comparison uses the same external-training-selected genes
and the same edge test cells for within- and cross-acquisition heads:
**0.353957 → 0.248129** macro Gene Pearson. It has a different denominator from
the full-acquisition primary score. Training-set size differs too, so this is
not a pure causal estimate of domain shift.

Known related sections/material can occur across original sample-held-out
splits. A source-page-out lung sensitivity changes Gene Pearson **0.423849 →
0.227551**. The runtime can construct guarded folds and reselect HVGs after
excluding known related training material. These relations do not demonstrate
independent patients. New guarded runs are not historical table values.

## Legacy Xenium benchmark

The historical 30-acquisition within-sample cohort normalizes to 16 tissue labels;
it differs from the native-Xenium five-category extension. Its candidate axis
uses measured-panel overlap/metadata and markers (up to 600 symbols), then
training-ranked targets. The cross-sample campaign has 18 targets in six
categories and three training seeds. Primary summaries use 50 targets.

All primary top-50 target lists are Gene Expression. For 15 of 18 historical
cross targets, larger **secondary** target-list tails contain negative-control
features. Those full-target historical metrics are explicitly mixed-feature
results; they do not redefine the primary gene-only tables. Exact original
target order and feature types remain in the fold and packet plans.

HLDX010's old image-pair scores are superseded. Historical tables retain an
explicit status and **NA** in current presentation; old values are not presented
as results for the corrected liver image. Not all 71 current Xenium acquisitions
have a modern benchmark score. Legacy and modern protocol results are not pooled.

## Metric definitions

Predictions are on **log1p raw expression**, without clipping negative
predictions. The six primary metrics are:

| Metric | Definition |
| --- | --- |
| Gene Pearson | Pearson across test cells per gene, then finite-gene mean |
| Gene Spearman | Spearman across test cells per gene, then finite-gene mean |
| Cell Pearson | Pearson across selected genes per test cell, then finite-cell mean |
| F1 | Gene-wise detection F1; truth raw count > 0, predicted log1p > log1p(0.5) |
| R² | Gene-wise coefficient of determination, then finite-gene mean |
| RMSE | Root mean squared log1p error over the selected test cells/genes |

Negative R² is retained, including severe failures. Constant-target R² is
undefined. Modern protocols explicitly exclude constant truth/prediction gene
correlations. Legacy numerical conventions are retained for faithful replay.
HD cross scoring uses float64; modern within-acquisition scoring uses original
float32. Legacy metrics are evaluated on the original full axis before selecting
the primary target slice, preserving historical arithmetic/array layout.

Average folds/seeds within acquisition, then acquisitions within category, then
categories equally. The **23-category** descriptive mean combines HD18 with the
five additional Xenium-only categories; overlapping Xenium categories are not
added to HD category rows. It is not a 61-sample-weighted or all-124-sample mean.

## Confidence intervals

HD18 conditional intervals use 10,000 spatial block resamples, primary width
500 µm, seed 20261019, with shared block multiplicities across folds/models.
Multiple-acquisition categories resample acquisitions within species strata.
Spearman midranks are recomputed under multiplicities. Low block support is
flagged; single-acquisition categories below the declared support thresholds
have unavailable limits. The original 250/1000 µm sensitivities are archived.

HD cross intervals use 10,000 fixed-cohort test-tissue spatial reweightings,
primary width 500 µm, seed 20261020, with original acquisition weights, fitted
heads, and HVGs fixed. They are **not new-acquisition or donor-generalization
confidence intervals**. Unsupported metrics have withheld limits; SD across
acquisitions is a different quantity and is not relabeled as a CI.

No new bootstrap was run for this repository publication. Full historical
numerical interval tables and implementation/configuration remain in the data
package. No formal combined 23-category or generalization macro CI is claimed.
