#!/usr/bin/env python
"""Independently check and render the published category-weighted results."""
import argparse
import csv
import json
import math
from pathlib import Path

METRICS = ['Gene P', 'Coord P', 'KNN P', 'Gene S', 'Cell P', 'F1', 'R²', 'RMSE']
ROOT = Path(__file__).resolve().parents[1]


def read_csv(path):
    with Path(path).open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def summarize(root=ROOT):
    root = Path(root)
    rows = read_csv(root/'results/xe_extension_23_category_results.csv')
    expected = read_csv(root/'results/xe_extension_category_macro_means.csv')
    if len(rows) != 23 or len({r['Category'] for r in rows}) != 23:
        raise ValueError('Expected 23 unique category rows')
    hd = [r for r in rows if r['Data'] == 'Visium HD']
    xe = [r for r in rows if r['Data'] == 'Xenium data only']
    if len(hd) != 18 or len(xe) != 5:
        raise ValueError('Unexpected HD/Xenium cohort definition')
    means = []
    for group, wanted in zip([hd, xe, rows], expected):
        mean = {m: math.fsum(float(r[m]) for r in group)/len(group) for m in METRICS}
        mean.update(Category=wanted['Category'], n=sum(int(r['n']) for r in group),
                    Data='Category-weighted mean')
        if mean['n'] != int(wanted['n']):
            raise ValueError('Acquisition denominator mismatch')
        for m in METRICS:
            if not math.isclose(mean[m], float(wanted[m]), abs_tol=1e-12, rel_tol=0):
                raise ValueError(f'Macro mismatch: {wanted["Category"]} / {m}')
        means.append(mean)
    return rows, means


def markdown(rows, means):
    columns = ['Category', 'n', 'Data']+METRICS
    lines = ['# Full 23-category results', '',
             'CONCH + Ridge; HD18 plus five additional **Xenium data only** categories.',
             'Each macro weights categories equally. Acquisitions are not independent donors.',
             'Overlapping Xenium categories are not pooled into the HD rows. No combined CI is claimed.', '',
             '| '+' | '.join(columns)+' |',
             '| '+' | '.join(['---', '---:', '---']+['---:']*len(METRICS))+' |']
    for row in rows+means:
        cells = [str(row['Category']), str(row['n']), str(row['Data'])]
        cells += [f'{float(row[m]):.3f}' for m in METRICS]
        lines.append('| '+' | '.join(cells)+' |')
    lines += ['', 'Generated from the unrounded CSVs by `scripts/summarize_results.py`.',
              'Coord P: quadratic-coordinate Ridge Gene Pearson; KNN P: spatial-KNN Gene Pearson.',
              'See [protocols](../docs/protocols.md) and [all results](README.md) for scope and uncertainty.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path, help='Optional Markdown output file')
    args = parser.parse_args()
    rows, means = summarize(args.root)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(markdown(rows, means))
    print(json.dumps(dict(status='passed', categories=len(rows), acquisitions=means[-1]['n'],
                          aggregation='equal categories, computed from unrounded estimates', macros=means), indent=2))


if __name__ == '__main__':
    main()
