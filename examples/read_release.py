#!/usr/bin/env python
"""Read one real raw CSR observation without loading a whole acquisition index."""
import argparse
import json
import os
from pathlib import Path


def first_expression_row(index):
    import pyarrow.parquet as pq
    for batch in pq.ParquetFile(index).iter_batches(
            batch_size=8192, columns=['cell_id', 'h5ad_row', 'x_px', 'y_px']):
        for row in batch.to_pylist():
            if row['h5ad_row'] is not None:
                return row
    raise ValueError('No expression-linked observation found')


def read(root, sample_id):
    import anndata as ad
    import numpy as np
    import pyarrow.parquet as pq

    root = Path(root).expanduser().resolve()
    table = pq.read_table(root/'metadata/samples.parquet').to_pylist()
    selected = [r for r in table if r['sample_id'] == sample_id]
    if len(selected) != 1:
        raise ValueError(f'Expected one sample record for {sample_id}')
    sample = selected[0]
    cell = first_expression_row(root/sample['cell_index_file'])
    matrix = ad.read_h5ad(root/sample['expression_file'], backed='r')
    try:
        row = int(cell['h5ad_row'])
        if str(matrix.obs_names[row]) != str(cell['cell_id']):
            raise ValueError('Cell ID / H5AD row link mismatch')
        values = matrix.X[row:row+1]
        if not np.issubdtype(values.dtype, np.integer) or np.any(values.data < 0):
            raise ValueError('Expected raw nonnegative integer CSR counts')
        genes = pq.read_table(root/sample['genes_file']).to_pylist()
        if len(genes) != matrix.n_vars:
            raise ValueError('Measured feature axis does not match expression')
        positive = [dict(gene_id=str(genes[int(i)]['gene_id']),
                         feature_type=genes[int(i)].get('feature_type'), count=int(v))
                    for i, v in zip(values.indices[:10], values.data[:10])]
        result = dict(sample_id=sample_id, stain=sample['stain'],
                      segmentation_kind=sample['segmentation_kind'], shape=list(matrix.shape),
                      cell=cell, nonzero_features=int(values.nnz), raw_feature_count_sum=int(values.sum()),
                      first_nonzero_features=positive, source_url=sample['source_url'],
                      note='Feature counts include the measured feature axis; select Gene Expression explicitly for gene-only analysis.')
    finally:
        matrix.file.close()
    print(json.dumps(result, indent=2, default=str))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release-root', type=Path,
                        default=Path(os.environ['SMMC_RELEASE_ROOT']) if os.environ.get('SMMC_RELEASE_ROOT') else None)
    parser.add_argument('--sample', required=True)
    args = parser.parse_args()
    if args.release_root is None:
        parser.error('Set --release-root or SMMC_RELEASE_ROOT')
    read(args.release_root, args.sample)


if __name__ == '__main__':
    main()
