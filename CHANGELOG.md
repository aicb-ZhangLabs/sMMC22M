# Release notes

## Public code and documentation — 2026-10-29

- Initialize this repository for sMMC-40M, retaining the existing `sMMC22M` URL.
- Add the author-supplied main figure, exact release statistics, 124-acquisition
  catalog, schema, and compact QC/source tables.
- Publish the tested benchmark runtime, source-RNA/CONCH pipeline, metric kernels,
  five protocol descriptions, exact fold/target plans, and numerical results.
- Add offline demo, result-summary validation, release readers, scientific tests,
  and explicit requirements for the separately stored large benchmark assets.
- Keep public code publication separate from public data hosting: no Hugging Face
  dataset upload or data-host identifier is announced in this version.

## Frozen data version — 2026-10-28

- Unify 53 HD and 71 Xenium acquisitions; 38,472,192 gene-positive observations,
  including 2,041,153 HD nuclei, in 23 harmonized reporting categories.
- Integrate source-backed sample context, feature axes, image eligibility,
  completed HD annotation evidence, and corrected HLDX010 image pairing.
- Freeze exact benchmark cohorts/splits/targets, cached inputs, predictions,
  arithmetic profiles, source-material relations, and per-file checksums.
- Retain historical protocols separately. Superseded HLDX010 predictions remain
  archived but are invalid for the current image pair. Missing legacy seed
  predictions are labeled reconstructions, with numerical differences preserved.

The data version is an immutable local freeze. This code publication does not
rebuild the data, change historical scores, rerun all confidence intervals, or
change the paper manuscript.
