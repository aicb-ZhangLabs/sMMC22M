"""Synthetic offline smoke; scores have no biological interpretation."""
from pathlib import Path
import json

import numpy as np


def run(output):
    from smcc_benchmark import METRICS, protocol, score

    output = Path(output)
    rng = np.random.RandomState(20261028)
    n, dimensions, n_genes = 1600, 32, 80
    coords = rng.uniform(0, 2000, (n, 2)).astype(np.float32)
    features = rng.normal(size=(n, dimensions)).astype(np.float32)
    rate = np.exp(np.clip(features[:, :8] @ rng.normal(0, .15, (8, n_genes)), -2, 2))
    raw = rng.poisson(rate).astype(np.int64)
    raw[:, -1] = 0  # Includes an unexpressed feature to exercise target filtering.
    mod = protocol(Path(__file__).resolve().parent)
    train, test, split = mod.spatial_split(coords, 0, .2, .05)
    logged = np.log1p(raw.astype(np.float32))
    # One synthetic acquisition: training-only population log-variance and
    # stable gene-ID ties. This does not claim to simulate the release's cohorts.
    genes = [f'synthetic_gene_{i:03d}' for i in range(n_genes)]
    variance = np.var(logged[train].astype(np.float64), axis=0)
    selected = sorted((i for i in range(n_genes) if variance[i] > 0),
                      key=lambda i: (-variance[i], genes[i]))[:50]
    y = logged[:, selected]
    predictions = mod.make_predictions(features, y, coords, train, test, 16, 10., 8, 0)
    values = {name: score(mod, raw[test][:, selected], y[test], pred,
                          [genes[i] for i in selected], True)
              for name, pred in predictions.items()}

    def clean(value):
        if isinstance(value, dict):
            return {k: clean(v) for k, v in value.items()}
        if isinstance(value, (float, np.floating)):
            return float(value) if np.isfinite(value) else None
        return value

    result = clean(dict(status='passed', kind='synthetic_demo_not_reported_benchmark',
                        observations=n, selected_genes=len(selected),
                        n_train=len(train), n_test=len(test), split=split,
                        primary_metrics=METRICS, metrics=values))
    output.mkdir(parents=True, exist_ok=True)
    (output/'demo.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result, indent=2))
    return result
