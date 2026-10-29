#!/usr/bin/env python
"""Read a sample-only package; its asset paths are acquisition-root relative."""
import argparse
import json
from pathlib import Path


def main():
    import anndata as ad
    import pyarrow.parquet as pq
    from read_release import first_expression_row

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample-root', required=True, type=Path)
    args = parser.parse_args()
    root = args.sample_root.expanduser().resolve()
    records = pq.read_table(root/'metadata/samples.parquet').to_pylist()
    if len(records) != 1:
        raise ValueError('Expected one acquisition in a sample-local package')
    sample = records[0]
    # Sample-local cells.parquet, unlike the root sample table, lives in metadata/.
    cell = first_expression_row(root/'metadata/cells.parquet')
    data = ad.read_h5ad(root/sample['expression_file'], backed='r')
    try:
        row = int(cell['h5ad_row'])
        if str(data.obs_names[row]) != str(cell['cell_id']):
            raise ValueError('Cell ID / matrix row mismatch')
        values = data.X[row:row+1]
        print(json.dumps(dict(sample_id=sample['sample_id'], shape=list(data.shape),
                              cell=cell, raw_feature_count_sum=int(values.sum()),
                              nonzero_features=int(values.nnz), source_url=sample['source_url']), indent=2, default=str))
    finally:
        data.file.close()


if __name__ == '__main__':
    main()
