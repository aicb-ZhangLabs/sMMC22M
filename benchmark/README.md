# Full HD and Xenium benchmark reproduction

This folder publishes the running scientific code for the **2026-10-28** frozen
benchmark. The large data-side assets are stored separately. The public Git
repository is not a substitute for the complete processed data release.

## Installation and offline checks

Use **Python 3.10** in a fresh environment:

```bash
python -m pip install -r benchmark/requirements-cpu.txt
python benchmark/smcc_benchmark.py list
python benchmark/smcc_benchmark.py demo --output /tmp/smmc-demo
python scripts/summarize_results.py
python -m unittest discover -s tests -v
```

The synthetic demo trains predictors and exercises controls/metrics, using
generated inputs. It requires no checkpoint, network, or biological dataset and
does **not** reproduce a reported tissue score. `summarize_results.py` verifies
the published 23-category arithmetic directly from the included numerical tables.

CPU reference pins are NumPy 1.23.5, SciPy 1.8.1, pandas 1.4.3,
scikit-learn 1.3.0, h5py 3.9.0, and pyarrow 23.0.1. The legacy reference profile
is separate in `requirements-legacy.txt`; do not install both profiles into one
environment. Image/RNA examples additionally use `requirements-read.txt`.
Encoding uses `requirements-encoder.txt` and a compatible native OpenSlide/C++
runtime. `environment.yml` specifies the native image-reader runtime.

## Exact asset package

For biological replay or refitting, obtain the **complete frozen data-side
`benchmark/` package**. A public download URL has not yet been announced; see
[availability](../docs/data.md#availability). It must contain:

```text
benchmark/
  folds.json, tracks.json, prediction_index.json
  splits/        # exact row memberships, including buffers
  targets/       # original gene IDs/symbols and training-only group plans
  fixtures/      # sampled counts/expression, cell axes, cached CONCH features
  packets/       # original predictions and matching truth/IDs
  metadata/      # documented source-material relationships
  models/        # official checkpoint declaration, not model weights
  frozen/        # evaluator snapshot
  reference/     # original numerical summaries
  provenance/    # original configs/code/environments/conditional intervals
  vendor/        # encoder source snapshot used in the frozen experiments
  SEALED.json, package_manifest.json, SHA256SUMS
```

[`plans/`](plans/) publishes exact fold and prediction-index JSON plans for
inspection. It does not contain the referenced large Parquet/NPZ/HDF5 assets.
There are 596 fold plans and 696 original prediction packets; 692 packets are
valid for current image evidence and four old HLDX010 packets are superseded.
The 36 missing historical cross-sample seed predictions are indexed separately
as **reconstructions**, not original saved predictions.

```bash
export SMMC_RELEASE_ROOT=/path/to/2026-10-28
# Or set SMMC_BENCHMARK_ROOT=/path/to/benchmark for cache-only operations.
# --benchmark-root and --release-root override environment defaults.

python benchmark/smcc_benchmark.py replay --track hd18 --sample HD-001 \
  --output /tmp/smmc-replay
python benchmark/smcc_benchmark.py refit --fold hd18_HD-001_seed0 \
  --output /tmp/smmc-refit
```

`replay` recomputes six primary metrics from **original** saved predictions and
checks them against frozen references at the unchanged absolute **5e-5** bound.
It does not retrain a model or re-encode images. `--limit 1` limits replay to one
prediction packet, not one cell. By default all eligible original packets are
checked. Superseded image pairs are excluded.

`refit` fits the predictor again on exactly frozen rows and targets from cached
CONCH features. It requires no GPU/checkpoint. Its output records the numerical
environment and new scores; a refit is not assumed bitwise identical to the
original training run. `--fold` identifies one plan in `plans/folds.json`.

## Track definitions

| Track | Scope | Targets and test set |
| --- | --- | --- |
| `hd18` | 49 H&E acquisitions in the cache; 48 in the 18-category primary table | Four directional spatial holdouts, 50 training-only organ/species/fold HVGs; mixed TMA HD-052 is auxiliary |
| `hd_cross` | 34 targets in 10 organ/species strata | Whole-acquisition LOSO, matched-edge comparisons, and source-page-out sensitivity; use each plan's regime |
| `xenium5` | 13 acquisitions in five additional categories | Native Xenium RNA, four spatial holdouts, 50 training-only HVGs; Skin exact-panel amendment |
| `xenium_legacy` | Historical 30-acquisition campaign | Original metadata-selected candidate-symbol axis and training-ranked targets; HLDX010 image scores invalid |
| `xenium_legacy_cross` | Historical 18 targets / six categories / three seeds | Archived cross-acquisition inputs and targets; only seed-0 predictions were originally retained |

These are acquisition counts, not independent-patient counts. Current HD18 and
Xenium5 results do not constitute an all-124-acquisition prediction experiment.
Legacy and modern protocols are preserved separately; no silent score pooling.

## Model and inputs

Modern inputs are **150 × 150 native-pixel RGB H&E crops**, with centres rounded
ties-to-even, half-open windows, and white edge padding. Physical crop size
varies with image MPP. The official 448-pixel CONCH transform gives 512-dimensional
**preprojection, unnormalized** frozen features. Targets are **log1p(raw integer
Gene Expression counts)**, without library-size normalization. HD counts use the
strict majority-area 2 µm aggregation; Xenium uses native cellular counts.

The image head is StandardScaler → randomized PCA (up to 128 components) → Ridge
(`alpha=10`). Spatial controls are quadratic-coordinate Ridge and training-only
8-neighbour expression averaging. The training mean is a separate control.
No gene measurements from a held-out acquisition enter training-HVG ranking.

See [protocols](../docs/protocols.md) for sampling, exact gene-panel handling,
float32/float64 arithmetic, aggregation order, metrics, and confidence intervals.

## Regenerating modern inputs from raw release assets

Install the encoder dependencies. Download the official checkpoint under its
upstream terms; no checkpoint is bundled here. The reference is
`MahmoodLab/CONCH`, revision `f9ca9f877171a28ade80228fb195ac5d79003357`,
file `pytorch_model.bin`, 802,235,437 bytes, SHA256
`40a9644b9ba0e83a74576e0a5e5f7313599fa9c9cdaf3c20f8a3e271b0e9ae7c`.
The full data-side benchmark retains the exact encoder implementation snapshot.

```python
from huggingface_hub import hf_hub_download
checkpoint = hf_hub_download(
    repo_id="MahmoodLab/CONCH",
    revision="f9ca9f877171a28ade80228fb195ac5d79003357",
    filename="pytorch_model.bin",
)
print(checkpoint)
```

```bash
# Set CHECKPOINT to the returned local path. RUN must be outside the release.
python -m pip install -r benchmark/requirements-encoder.txt
export RUN=/tmp/smmc-raw-run
export CHECKPOINT=/path/to/pytorch_model.bin

# Fast real source-RNA/image smoke, 32 original seed-0 test cells.
python benchmark/smcc_benchmark.py prepare --track hd18 --sample HD-001 \
  --limit 32 --output "$RUN"
python benchmark/smcc_benchmark.py encode --track hd18 --sample HD-001 \
  --checkpoint "$CHECKPOINT" --device cuda --output "$RUN"

# Full cohort preparation/encoding for a true predictor refit:
python benchmark/smcc_benchmark.py prepare --track hd18 --sample all --output "$RUN"
python benchmark/smcc_benchmark.py encode --track hd18 --sample all \
  --checkpoint "$CHECKPOINT" --device cuda --output "$RUN"
python benchmark/smcc_benchmark.py select-hvgs --track hd18 \
  --inputs-root "$RUN" --output "$RUN"
python benchmark/smcc_benchmark.py refit --fold hd18_HD-001_seed0 \
  --inputs-root "$RUN" --output /tmp/smmc-fresh-refit
```

For the native-Xenium extension, use `--track xenium5` and `--sample XE-055`, or
`all` for the full cohort. `--device cpu` is supported for encoding; it can be
slow. `--disable-cudnn` uses native CUDA convolution if a system's cuDNN ABI is
incompatible. Use allocated compute resources for full image encoding.

`prepare` verifies exact cell/gene axes and integer RNA against the frozen cache.
`encode` verifies the checkpoint hash and new features against frozen features.
`select-hvgs` recomputes training-only variances for complete acquisition groups
and requires exact target agreement. A 32-cell smoke cannot replace a full
cohort or full HVG/refit run. Fresh preparation/encoding/refit is implemented for
modern spatial tracks; other tracks use their documented frozen input snapshots.

## Known related source material

```bash
python benchmark/smcc_benchmark.py groups --output /tmp/smmc-material-groups
python benchmark/smcc_benchmark.py refit --fold loso_HD-010 \
  --material-group-guard --output /tmp/smmc-guarded
```

Default cross-acquisition refit rejects folds with known shared source material.
`--material-group-guard` excludes related training acquisitions and **reselects
HVGs** from the remaining training acquisitions before fitting. This is a new
guarded protocol, not a replacement historical score. An explicit
`--allow-related-material` reproduces the disclosed historical setup. No remaining
training acquisition gives an unsupported fold rather than an invented result.
Unknown donor relationships remain unknown.

## Validation and boundaries

The data freeze checked 692 current-valid original packets, 2,112 model rows and
12,672 primary metric comparisons. Separately, 36 missing legacy seed
reconstructions were labeled and checked; one image R² differed by
5.49626e-5, above the original replay bound. The original score was not changed
and the strict replay tolerance was not relaxed. Reconstructed predictions are
not included in strict original-prediction replay.

Raw-source checks at the data freeze included 32 HD-001 and 32 XE-055 cells
through RNA/image → CONCH → predictor, and eight HVG selections. They were
selected end-to-end checks, not a full re-encoding of the data library. The
public-code validation checks are documented in
[`validation/publication_checks.json`](../validation/publication_checks.json).

Historical interval tables/configuration and original CI programs remain in
the data package's `provenance/*conditional_intervals/`. The public
[`intervals/cross_sample_interval_kernel.py`](intervals/cross_sample_interval_kernel.py)
is the unchanged numerical spatial-reweighting kernel, not a standalone
reproduction of every historical aggregation. Replay validates point metrics;
it does not rerun bootstrap intervals. There is no new 23-category macro CI.
