"""Checks for spatial leakage, undefined metrics, weighted ranks, and scope."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'benchmark'))
sys.path.insert(0, str(ROOT/'scripts'))
from smcc_benchmark import protocol, refit, score
from summarize_results import summarize


class ScientificContracts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = protocol(ROOT/'benchmark')

    def test_four_spatial_holdouts_have_separate_buffer(self):
        coords = np.stack(np.meshgrid(np.arange(40), np.arange(40)), axis=-1).reshape(-1, 2)
        seen = set()
        for seed in range(4):
            train, test, info = self.mod.spatial_split(coords, seed, .2, .05)
            self.assertFalse(set(train) & set(test))
            axis = 0 if info['axis'] == 'x' else 1
            buffer = .05 * np.ptp(coords[:, axis])
            if info['side'] == 'high':
                self.assertLess(float(coords[train, axis].max()), info['cut']-buffer)
                self.assertGreaterEqual(float(coords[test, axis].min()), info['cut'])
            else:
                self.assertGreater(float(coords[train, axis].min()), info['cut']+buffer)
                self.assertLessEqual(float(coords[test, axis].max()), info['cut'])
            seen.add((info['axis'], info['side']))
        self.assertEqual(len(seen), 4)

    def test_predictors_never_read_test_expression(self):
        rng = np.random.RandomState(7)
        coords = rng.uniform(0, 1000, (600, 2)).astype(np.float32)
        x = rng.normal(size=(600, 12)).astype(np.float32)
        y = rng.normal(size=(600, 8)).astype(np.float32)
        train, test, _ = self.mod.spatial_split(coords, 0, .2, .05)
        original = self.mod.make_predictions(x, y, coords, train, test, 8, 10., 8, 0)
        changed = y.copy()
        changed[test] = rng.normal(1000, 100, changed[test].shape)
        repeated = self.mod.make_predictions(x, changed, coords, train, test, 8, 10., 8, 0)
        for name in original:
            np.testing.assert_array_equal(original[name], repeated[name], err_msg=name)

    def test_perfect_predictions_and_constant_targets(self):
        raw = np.array([[0, 2, 0], [1, 0, 0], [3, 1, 0], [2, 4, 0]], dtype=np.int64)
        logged = np.log1p(raw.astype(np.float32))
        values = score(self.mod, raw, logged, logged.copy(), ['a', 'b', 'constant'], True, 'float64')
        self.assertAlmostEqual(values['pearson_gene_macro'], 1.)
        self.assertAlmostEqual(values['spearman_gene_macro'], 1.)
        self.assertAlmostEqual(values['r2_gene_macro'], 1.)
        self.assertEqual(values['rmse_log1p'], 0.)
        constant = np.repeat(logged.mean(axis=0, keepdims=True), len(raw), axis=0)
        self.assertTrue(np.isnan(score(self.mod, raw, logged, constant, ['a','b','constant'], True)['pearson_gene_macro']))

    def test_negative_r2_is_preserved(self):
        raw = np.arange(1, 13).reshape(6, 2)
        logged = np.log1p(raw.astype(np.float32))
        values = score(self.mod, raw, logged, logged+100, ['a','b'], True, 'float64')
        self.assertLess(values['r2_gene_macro'], -1000)

    def test_weighted_spearman_matches_literal_cell_repetition(self):
        source = ROOT/'benchmark/intervals/cross_sample_interval_kernel.py'
        spec = importlib.util.spec_from_file_location('ci_kernel_test', source)
        kernel = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(kernel)
        truth = np.array([[0,1],[1,0],[1,2],[2,1],[3,4],[0,2]], dtype=float)
        pred = np.array([[1,0],[0,1],[0,1],[3,2],[3,3],[2,2]], dtype=float)
        blocks = np.array([0,0,1,1,2,2])
        weights = np.array([[1,1,1],[2,0,1],[0,1,2]])
        got, eligible = kernel.exact_spearman(truth, pred, blocks, weights)
        for i, w in enumerate(weights):
            take = np.repeat(np.arange(len(truth)), w[blocks])
            expected = np.nanmean([spearmanr(truth[take,j], pred[take,j]).statistic
                                   if hasattr(spearmanr(truth[take,j], pred[take,j]), 'statistic')
                                   else spearmanr(truth[take,j], pred[take,j]).correlation
                                   for j in range(truth.shape[1])])
            self.assertAlmostEqual(got[i], expected, places=12)

    def test_superseded_image_pair_rejected_before_input_access(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path/'folds.json').write_text(json.dumps([dict(fold_id='old_liver', valid_for_current_release=False)]))
            with self.assertRaisesRegex(ValueError, 'invalid/superseded'):
                refit(path, 'old_liver', path/'output')

    def test_reported_macros_are_unrounded_category_means(self):
        rows, means = summarize(ROOT)
        self.assertEqual((len(rows), means[-1]['n']), (23, 61))
        self.assertAlmostEqual(means[-1]['Gene P'], .351293415401868, places=14)
        sample_weighted = sum(float(r['Gene P'])*int(r['n']) for r in rows)/61
        self.assertGreater(abs(sample_weighted-means[-1]['Gene P']), .001)


if __name__ == '__main__':
    unittest.main()
