# Public-code validation

[`publication_checks.json`](publication_checks.json) records checks of this code
against the complete immutable **2026-10-28** local data-side package:

- Seven scientific-contract tests, including spatial buffers, no test-expression
  leakage into predictors, constant/negative metrics, exact weighted Spearman,
  superseded image rejection, and unrounded category averaging.
- One original prediction packet replayed from each of the five tracks.
- Three real predictor refits: HD-001 within-acquisition, XE-055
  within-acquisition, and HD-001 cross-acquisition. All compared metrics met the
  original 5e-5 bound.
- Actual raw-RNA preparation for four cells in each modern track, two real H&E
  image crops, four image/object-branch readers, and two sample-only readers.
- Rejection of missing assets, invalid HLDX010 image evidence, known related
  material without an explicit protocol choice, and output inside the immutable
  data package.

The scientific runtime functions, pipeline, evaluators, interval kernel, and
published source CSVs were compared with the frozen inputs. The main image was
checked byte-for-byte against the supplied original. Repository links, syntax,
statistics, source lineage, private machine paths/credential patterns, and file
sizes were checked before publication.

The full freeze had already validated 692 eligible original packets and 12,672
metric values. This publication performs selected end-to-end checks; it does not
claim another complete replay, all-image re-encoding, a new biological accuracy
audit, or public data-download verification.

[`source_checksums.json`](source_checksums.json) records the public files before
that checksum file was written; it excludes itself. The exact Git commit pins
the complete published tree.
