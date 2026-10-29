# Full benchmarking results

## Modern coverage summary

[Readable 23-category table](benchmark_23_categories.md) ·
[Unrounded CSV](xe_extension_23_category_results.csv) ·
[Macro means](xe_extension_category_macro_means.csv)

The 18 HD categories include 48 acquisitions and 564,433 sampled cells. The
five additional Xenium categories include 13 acquisitions and 156,000 sampled
cells. Their combined table is a **23-category descriptive macro**, not a
platform-controlled comparison or an experiment on all 124 acquisitions. Modern
overlapping Xenium categories are not pooled into the HD rows. Category labels
are reporting strata, not uniformly anatomical organs.

Recalculate all three macro rows from unrounded values:

```bash
python scripts/summarize_results.py
python scripts/summarize_results.py --output /tmp/benchmark_23_categories.md
```

## Full HD benchmark result

| File | Contents |
| --- | --- |
| [hd18_full_results.csv](hd18_full_results.csv) | All 18 primary categories, six metrics, spatial controls, and conditional interval/status fields |
| [hd18_macro.csv](hd18_macro.csv) | All four predictor/control macro rows |
| [hd18_auxiliary_sample.csv](hd18_auxiliary_sample.csv) | HD-052 mixed TMA, excluded from primary category means |
| [hd_cross_cross_sample_category_results.csv](hd_cross_cross_sample_category_results.csv) | Primary whole-acquisition results for 10 organ/species strata / 34 targets, intervals and across-acquisition SD |
| [hd_cross_by_sample.csv](hd_cross_by_sample.csv) | Acquisition-level cross-validation results |
| [hd_cross_cross_sample_macro.csv](hd_cross_cross_sample_macro.csv) | Primary full-acquisition macro; Gene Pearson 0.2603575231 |
| [hd_cross_matched_by_category.csv](hd_cross_matched_by_category.csv) | Matched-edge within/cross comparison |
| [hd_cross_matched_macro.csv](hd_cross_matched_macro.csv) | Matched-edge macro; Gene Pearson 0.353957 → 0.248129 |
| [hd_cross_source_page_sensitivity.csv](hd_cross_source_page_sensitivity.csv) | Linked-material/source-page sensitivity |
| [hd_cross_gene_metric_support.csv](hd_cross_gene_metric_support.csv) | Target metric support diagnostics |

Primary full-acquisition and matched-edge scores have **different test
denominators**. Keep them distinct. Extremely negative R² is retained: the
pancreas primary cross-acquisition R² is approximately −1400.866, a major
calibration failure rather than a value to clip or omit.

The reported source-page lung sensitivity is Gene Pearson 0.423849 → 0.227551.
Acquisition-held-out results do not establish donor independence. New guarded
refits exclude known related material and select training HVGs again; those
outputs are not substituted into the historical tables.

## Full Xenium benchmark result

| File | Contents |
| --- | --- |
| [xe_extension_by_category.csv](xe_extension_by_category.csv) | Modern five-category extension, all model/control metrics |
| [xe_extension_by_sample.csv](xe_extension_by_sample.csv) | All 13 extension acquisitions |
| [protocol/xe_extension_protocol_amendment.json](protocol/xe_extension_protocol_amendment.json) | Skin panel grouping before scoring |
| [legacy_xenium_within_primary.csv](legacy_xenium_within_primary.csv) | Historical within-acquisition primary fold/model rows |
| [printed_legacy_category_values.csv](printed_legacy_category_values.csv) | Historical 16-label regrouping with explicit status |
| [legacy_xenium_cross_primary.csv](legacy_xenium_cross_primary.csv) | Historical 18-target / six-category / three-seed results |

**HLDX010:** raw historical source tables retain its old numerical values for
traceability, with a superseded status in the regrouped table. Those values
come from the wrong image pair. Treat the liver row as **NA for current image
evidence**, as in the manuscript; do not report it as a corrected benchmark.

The 30-acquisition historical protocol and 13-acquisition modern extension have
different input/target policies. They are separate results. Primary top-50
targets are genes; some larger legacy cross-target tails include controls.
The 36 reconstructed missing-seed predictions are not original saved packets,
and are not silently substituted into these original tables.

## Intervals and provenance

HD interval columns describe fixed-predictor spatial resampling under the
original declared rules; across-acquisition SD is not a confidence interval.
Unsupported limits remain unavailable. No combined 23-category or new-donor
macro CI is claimed. Original full interval tables, source code, and settings
remain in the complete data package.

CSV result files are byte-identical to the frozen manuscript reference tables.
They retain unrounded values and original column names. Their source mapping is
in [`source_lineage.json`](../data/source_lineage.json). Additional CSVs in this
directory expose per-fold, acquisition, source-page, and mapping details.
See [protocols](../docs/protocols.md) for aggregation, arithmetic, and limitations.
