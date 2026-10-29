# Contributing

Please open a repository issue with the release version, protocol track,
sample/fold IDs, software environment, and a small reproducible example.
Do not attach patient data, credentials, model weights, whole slides, or the
full expression/feature cache to issues or Git commits.

For code changes, run:

```bash
python -m pip install -r benchmark/requirements-cpu.txt
python -m unittest discover -s tests -v
python benchmark/smcc_benchmark.py demo --output /tmp/smmc-demo
python scripts/summarize_results.py
```

With the complete data-side benchmark package, also replay the affected track
and compare a real predictor refit. Keep published source results immutable.
Changes to cell membership, target selection, image preprocessing, metric
arithmetic, biological annotation, or source-material exclusion define a new
version/protocol and should not silently replace an old score.

Use explicit missing/unsupported states rather than fabricated zeros or cell
labels. Sample context and predicted cell types are not cell-level clinical
ground truth. Preserve upstream source attribution and applicable usage terms.
