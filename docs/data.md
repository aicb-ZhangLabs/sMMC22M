# Data version, structure, and reading

## Availability

The frozen data version is **2026-10-28**, schema **2.0.0**. This Git repository
contains code, documentation, source metadata, and summary results. The full
versioned data folder is approximately **931.6 GB** in logical size, including
images, masks, expression, and reproducibility assets. It has been assembled and
checked locally; **a public data upload has not yet been performed**.

No Hugging Face dataset ID or dataset-download command is provided until the
actual remote release exists. Original acquisition sources are linked in
[`sample_catalog.csv`](../data/sample_catalog.csv); upstream downloads alone do
not contain this resource's processed masks, aggregated RNA, metadata, or frozen
benchmark partitions. The complete processed package is needed for exact replay.

The data-side benchmark folder includes approximately 2.93 GB of sampled input
fixtures, 138 MB of membership Parquet files, original prediction packets, exact
target lists, and frozen encoder code. These are separate from the small code
repository. Model weights are acquired from the official upstream source.

## Counts and denominators

| Quantity | Count | Meaning |
| --- | ---: | --- |
| Acquisitions | 124 | 53 HD + 71 Xenium; not 124 independent donors/studies |
| Indexed candidates | 46,284,626 | Every segmentation/native candidate |
| Matrix observations with nonzero features | 38,486,654 | May include only non-gene controls/protein |
| Gene-positive observations | 38,472,192 | At least one positive Gene Expression count |
| Gene-positive whole cells | 36,431,039 | HD H&E + both Xenium branches |
| Gene-positive nuclei | 2,041,153 | HD IF; not an additional cell count |
| Technically eligible H&E cell pairs | 32,682,585 | Gene-positive, selected-image centroid/mask available |
| Complete 224 px H&E cell crops | 32,682,304 | Additional image-edge criterion |

Missing expression links mean **no matrix row**, not a measured all-zero gene
vector. Native measured-zero Xenium candidates have an explicit status. Keep
non-gene controls separate from Gene Expression. HD UMI and Xenium gene-assigned
transcript counts are distinct measurement types.

## Complete release layout

```text
2026-10-28/
  release.json                 # version and exact counts
  manifest.json, SHA256SUMS     # versioned file inventory
  metadata/
    samples.parquet            # 124-row starting index; root-relative paths
    schema.json
    biological_context.parquet
    field_evidence.parquet
    eligibility/{HD,Xenium}/    # per-candidate technical flags
    source_material_relationships.json
  HD/
    metadata/{cells,nuclei}/    # Parquet observation shards
    samples/HD-001/ ...
  Xenium/
    metadata/cells/
    samples/HLDX010/ ...        # corrected image pair is the default
  benchmark/
    folds.json, tracks.json
    splits/, targets/, fixtures/, packets/
    reference/, provenance/, vendor/, validation/
  qc/
  provenance/
```

The root [`schema.json`](../data/schema.json) explains path bases. The root sample
table uses paths relative to the release folder. Preserved platform-global cell
shards use HD- or Xenium-relative H5AD paths; resolve them with the root sample
record rather than concatenating paths blindly. Sample-local metadata resolves
paths relative to the acquisition folder, not its `metadata/` child.

## Expression and metadata

H5AD stores raw nonnegative integer counts as CSR on the original measured
feature axis. Use backed reading for metadata and selected rows; do not densify
millions of observations. There are 84 feature-axis definitions across the
release. A shared `gene_table_id` indicates an identical measured axis. Missing
panel genes are **unmeasured**, not zero. Avoid a zero-filled global gene union.

The stable join is `sample_id` plus `cell_id` at observation level, and
`sample_id` for sample context. Some preserved HD cell shards say
`platform="Visium HD"`, while the root table uses `HD`; aliases are declared in
[`observation_views.json`](../data/observation_views.json). `h5ad_row` is nullable
for candidates without RNA. Original feature types, gene IDs, and symbols remain
separate. `obsm/spatial` is in selected-image pixels; `obsm/spatial_um` is in µm.

```bash
python -m pip install -r benchmark/requirements-read.txt
python examples/read_release.py --release-root /path/to/2026-10-28 --sample HD-001
python examples/read_sample.py --sample-root /path/to/2026-10-28/HD/samples/HD-001
```

Both readers retrieve one real sparse RNA row and its metadata. The release
reader uses bounded Parquet batches; it does not load the whole cell index into
memory. A sample-only reader does not require the full resource.

## Context and annotations

The catalog retains original tissue labels and harmonized reporting categories.
Mixed HD-052 has multiple source-backed tissue categories and no invented
per-cell organ assignment. Do not pool developmental, whole-organism, plant, or
xenograft labels into a count of anatomical organs.

Clinical fields are sample/tissue context. Missing values have explicit
reported/not-reported/not-applicable/conflict/inherited status where available;
unknown does not imply healthy. Age is available for 34 acquisitions and sex for
41; one acquisition has a verified public paper-scoped patient label. Five known
material groups span 16 acquisitions. These relations cannot establish that all
other acquisitions represent independent donors.

HD annotation evidence includes marker support, high/medium/low prediction
tiers, uncalibrated scores, version, profile, and `Unresolved`. Predicted labels
are not manually annotated ground truth; sample tumor status is not a cell-level
malignancy label. Xenium retains its original annotation status.

## Integrity and versioning

The root release manifest SHA256 is
`f4f20ea1a5acfe660f1988184337921a92e1f5fcb519609324bfb921cfbda460`.
The full release contains nested manifests and per-file checksums. New/changed
files were freshly hashed; unchanged large inputs retained earlier verified
hashes under file-identity guards. This was not a new TB-scale full rehash.
Use a new version for changes and never modify an immutable data package in place.

Public projections and copied source files are described in
[`source_lineage.json`](../data/source_lineage.json). No raw biological arrays,
private machine paths, credentials, or weights are required to inspect the public
catalog and result summaries.
