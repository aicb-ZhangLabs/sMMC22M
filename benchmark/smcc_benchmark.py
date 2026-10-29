"""Relocatable frozen-feature benchmark and saved-prediction replay.

Public runtime for the 2026-10-28 frozen data-side benchmark package.
Set SMMC_RELEASE_ROOT or pass --benchmark-root for exact replay/refit.
The offline demo and track list do not require biological data assets.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

METRICS = ['pearson_gene_macro', 'spearman_gene_macro', 'pearson_cell_macro',
           'f1_gene_macro', 'r2_gene_macro', 'rmse_log1p']


def json_read(path):
    return json.loads(Path(path).read_text())


def protocol(root):
    # Use the source version shipped with this code, not an arbitrary external
    # module supplied alongside data. This evaluator is an unchanged snapshot.
    p = Path(__file__).resolve().parent / 'frozen/evaluate_spatial.py'
    spec = importlib.util.spec_from_file_location('smcc_frozen_spatial', p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def score_float32(mod, raw, logged, pred, genes, modern):
    per_gene = mod.metric_rows(raw, logged, pred, genes)
    if modern:
        constant = (np.ptp(pred, axis=0) == 0) | (np.ptp(logged, axis=0) == 0)
        per_gene.loc[constant, ['pearson', 'spearman']] = np.nan
    return mod.summarize(per_gene, logged, pred)


def score(mod, raw, logged, pred, genes, modern, metric_dtype='float32'):
    if metric_dtype=='float32':
        return score_float32(mod,raw,logged,pred,genes,modern)
    if metric_dtype!='float64':raise ValueError('Unknown metric arithmetic dtype '+metric_dtype)
    # The original HD cross-acquisition evaluator explicitly uses float64 for
    # all scoring and excludes R² for exactly constant target columns.
    raw=raw.astype(np.float64);logged=logged.astype(np.float64);pred=pred.astype(np.float64)
    per_gene=mod.metric_rows(raw,logged,pred,genes)
    constant=(np.ptp(pred,axis=0)==0)|(np.ptp(logged,axis=0)==0)
    per_gene.loc[constant,['pearson','spearman']]=np.nan
    per_gene.loc[np.ptp(logged,axis=0)==0,'r2']=np.nan
    return mod.summarize(per_gene,logged,pred)


def legacy_primary_score(mod,raw,logged,pred,genes,primary=50):
    """Exact legacy sequence: full-feature metrics, then primary-column subset."""
    per_gene=mod.metric_rows(raw,logged,pred,genes)
    take=np.arange(min(primary,len(genes)))
    return mod.summarize(per_gene.iloc[take],logged[:,take],pred[:,take])


def replay(root, track=None, sample=None, limit=None, reuse_paths=None):
    mod = protocol(root)
    records = json_read(root / 'prediction_index.json')
    selected = [r for r in records if (track is None or r['track'] == track)
                and (sample is None or r['sample_id'] == sample)
                and r.get('valid_for_current_release', True)]
    if limit:
        selected = selected[:limit]
    if not selected:
        raise ValueError('No eligible prediction packets match the requested track/sample')
    checked = []
    reused=[r for r in selected if reuse_paths is not None and r['path'] in reuse_paths]
    for record in selected:
        if record in reused:continue
        with np.load(root / record['path'], allow_pickle=False) as packet:
            n = record.get('primary_genes', len(packet['gene_ids']))
            raw = packet['y_true_raw'] if not record['modern_metrics'] else packet['y_true_raw'][:, :n]
            logged = packet['y_true_log1p'] if not record['modern_metrics'] else packet['y_true_log1p'][:, :n]
            genes = packet['gene_ids'].astype(str).tolist() if not record['modern_metrics'] else packet['gene_ids'][:n].astype(str).tolist()
            assert len(set(packet['release_cell_ids'].astype(str))) == len(raw)
            for model, expected in record['expected'].items():
                pred = packet['pred_' + model] if not record['modern_metrics'] else packet['pred_' + model][:, :n]
                values = legacy_primary_score(mod,raw,logged,pred,genes,n) if not record['modern_metrics'] else score(mod, raw, logged, pred, genes, True,record.get('metric_dtype','float32'))
                delta = []
                for metric in METRICS:
                    want = expected.get(metric)
                    got = values[metric]
                    if want is None or not np.isfinite(want):
                        assert not np.isfinite(got), (record['path'], model, metric, want, got)
                    else:
                        assert np.isfinite(got), (record['path'], model, metric, want, got)
                        error = abs(float(got) - float(want))
                        assert error <= 5e-5, (record['path'], model, metric, want, got, error)
                        delta.append(error)
                checked.append({'packet': record['path'], 'track': record['track'],
                                'sample_id': record['sample_id'], 'model': model,
                                'maximum_absolute_error': max(delta, default=0), **{m:float(values[m]) for m in METRICS}})
        print(json.dumps({'replayed': record['path'], 'models': len(record['expected'])}), flush=True)
    reused_models=sum(len(r['expected']) for r in reused)
    return {'status': 'passed', 'packets': len(selected), 'model_rows': len(checked)+reused_models,
            'newly_replayed_packets':len(selected)-len(reused),'reused_verified_packets':len(reused),
            'reused_verified_model_rows':reused_models,'reused_verified_error_bound':5e-5 if reused else None,
            'metric_checks': (len(checked)+reused_models)*len(METRICS),
            'max_absolute_error': max(r['maximum_absolute_error'] for r in checked), 'rows': checked}


def cached_sample(root, fixture, sid, inputs_root=None):
    base = (inputs_root if inputs_root is not None else root) / 'fixtures' / fixture / sid
    cells = pd.read_parquet(base / 'cells.parquet')
    with h5py.File(base / 'features.h5', 'r') as handle:
        x = np.asarray(handle['embeddings'][:], dtype=np.float32)
        ids = [v.decode() if isinstance(v, bytes) else str(v) for v in handle['cell_ids'][:]]
        assert ids == cells['benchmark_cell_id'].astype(str).tolist()
        if 'expressions' in handle:
            y = np.asarray(handle['expressions'][:], dtype=np.float32)
            genes = [v.decode() if isinstance(v, bytes) else str(v) for v in handle['gene_names'][:]]
        else:
            y = sparse.load_npz(base / 'counts.npz')
            genes = pd.read_parquet(base / 'genes.parquet')['gene_id'].astype(str).tolist()
    return cells, x, y, genes


def target_columns(raw, all_genes, chosen):
    lookup = {gid:i for i,gid in enumerate(all_genes)}
    assert len(lookup) == len(all_genes)
    indices = [lookup[x] for x in chosen]
    out = raw[:, indices]
    return out.toarray() if sparse.issparse(out) else out


def feature_types(root,fixture,sid,chosen):
    if fixture!='xenium_legacy':return ['Gene Expression']*len(chosen)
    mapping=pd.read_parquet(root/'fixtures/xenium_legacy'/sid/'gene_mapping.parquet')
    lookup=dict(zip(mapping.legacy_gene_symbol,mapping.feature_type))
    return [lookup.get(g,'Unknown') or 'Unknown' for g in chosen]


def group_map(root):
    out={}
    for group in json_read(root/'metadata/source_material_relationships.json'):
        for sid in group['sample_ids']:
            if sid in out:raise ValueError('A sample belongs to multiple material groups')
            out[sid]=group['source_material_group_id']
    return out


def group_audit(root, output):
    groups=group_map(root);report=[];guarded=[]
    for f in json_read(root/'folds.json'):
        if f['regime']=='spatial':continue
        target_group=groups.get(f['sample_id'])
        removed=[s for s in f['train_sample_ids'] if target_group is not None and groups.get(s)==target_group]
        remaining=[s for s in f['train_sample_ids'] if s not in removed]
        report.append({'fold_id':f['fold_id'],'track':f['track'],'test_sample_id':f['sample_id'],
                       'source_material_group_id':target_group,'known_material_group_overlap':bool(removed),
                       'overlapping_train_samples':removed,'guarded_train_sample_ids':remaining,
                       'guarded_split_supported':bool(remaining),'independent_donors_verified':False})
        output.mkdir(parents=True,exist_ok=True)
        original_membership=pd.read_parquet(root/f['membership_file'])
        changed=original_membership.copy()
        excluded=(changed['partition']=='train') & changed['sample_id'].isin(removed)
        changed.loc[excluded,'partition']='excluded_related_material'
        membership_path=output/'memberships'/(f['fold_id']+'.parquet')
        membership_path.parent.mkdir(parents=True,exist_ok=True)
        changed.to_parquet(membership_path,index=False,compression='zstd')
        stored_membership=str(membership_path.relative_to(root)) if membership_path.is_relative_to(root) else str(membership_path.resolve())
        guarded.append(dict(f,train_sample_ids=remaining,removed_related_samples=removed,
                            gene_ids=None, original_gene_ids=f['gene_ids'],
                            membership_file=stored_membership,
                            status='requires_training_only_HVG_reselection_and_refit' if remaining else 'unsupported_no_unrelated_training_acquisition',
                            new_protocol=True, independent_donors_verified=False))
    output.mkdir(parents=True,exist_ok=True)
    (output/'material_guarded_fold_plans.json').write_text(json.dumps(guarded,indent=2)+'\n')
    result={'status':'passed','cross_folds_checked':len(report),
            'historical_folds_with_known_related_material':sum(r['known_material_group_overlap'] for r in report),
            'guarded_supported_folds':sum(r['guarded_split_supported'] for r in report),
            'scope':'Known source-material groups only; unknown donor IDs do not establish donor independence. Original published splits remain immutable. New guarded refits reselect targets using remaining training samples.',
            'folds':report}
    return result


def refit(root, fold_id, output, material_group_guard=False, allow_related_material=False, inputs_root=None):
    """Refit an exactly frozen fold from shared cached CONCH features on CPU."""
    folds = json_read(root / 'folds.json')
    fold = dict(next(x for x in folds if x['fold_id'] == fold_id))
    if not fold['valid_for_current_release']:
        raise ValueError('This is an invalid/superseded historical image pair')
    mod = protocol(root)
    if inputs_root is not None and fold['track'] not in ['hd18','xenium5']:
        raise ValueError('Fresh raw-input refits currently cover modern hd18/xenium5 spatial tracks; other tracks use their fully frozen input snapshots')
    cells, x, y, genes = cached_sample(root, fold['fixture'], fold['sample_id'],inputs_root)
    if inputs_root is not None:
        reference=pd.read_parquet(root/'fixtures'/fold['fixture']/fold['sample_id']/'cells.parquet')
        if cells.cell_id.tolist()!=reference.cell_id.tolist():raise ValueError('Refit requires the full frozen cohort; --limit outputs are preparation/encoding smokes only')
    chosen = fold['gene_ids']
    if inputs_root is not None and not material_group_guard:
        selected=json_read(inputs_root/'hvgs'/fold['fixture']/fold['sample_id']/(f"fold_{fold['seed']}.json"))
        assert selected['gene_ids']==chosen,'Fresh training-only HVGs changed'
    if fold['regime']!='spatial':
        gg=group_map(root);target_group=gg.get(fold['sample_id'])
        related=[s for s in fold['train_sample_ids'] if target_group is not None and gg.get(s)==target_group]
        if related and not material_group_guard and not allow_related_material:
            raise ValueError('Historical sample-held-out fold shares known material: '+','.join(related)+'. Use --material-group-guard for a new guarded refit, or explicitly --allow-related-material to reproduce this historical score.')
    if material_group_guard:
        if fold['regime']=='spatial':raise ValueError('A material-group guard applies to cross-acquisition folds')
        groups=group_map(root);g=groups.get(fold['sample_id'])
        original=fold['train_sample_ids']
        fold['train_sample_ids']=[s for s in original if g is None or groups.get(s)!=g]
        if not fold['train_sample_ids']:raise ValueError('No unrelated training acquisition remains after the known-material guard')
        sources=[cached_sample(root,fold['fixture'],s,inputs_root) for s in fold['train_sample_ids']]
        common=set(genes)
        for _,_,_,gg in sources:common.intersection_update(gg)
        common=sorted(common) if fold['modern_metrics'] else [g for g in sources[0][3] if g in common]
        if fold['modern_metrics']:
            vv=np.zeros(len(common),dtype=np.float64)
            for _,_,rr,gg in sources:
                logged=np.log1p(target_columns(rr,gg,common).astype(np.float64))
                vv+=np.var(logged,axis=0)/len(sources)
            ii=sorted(range(len(common)),key=lambda i:(-vv[i],common[i]))[:50]
        else:
            raw_sources=np.concatenate([target_columns(rr,gg,common) for _,_,rr,gg in sources])
            ii=mod.top_variable_genes(raw_sources,200).tolist()
        chosen=[common[i] for i in ii]
        if len(chosen)<10:raise ValueError('Insufficient variable shared genes after material guard')
    raw = target_columns(y, genes, chosen)
    logged = np.log1p(raw.astype(np.float32)).astype(np.float32)
    membership = pd.read_parquet(root / fold['membership_file'])
    local = membership[membership['sample_id'] == fold['sample_id']]
    test = local.loc[local['partition'] == 'test', 'benchmark_row'].to_numpy(dtype=np.int64)
    if fold['regime'] == 'spatial':
        train = local.loc[local['partition'] == 'train', 'benchmark_row'].to_numpy(dtype=np.int64)
        coords = cells[['coord_x','coord_y']].to_numpy(dtype=np.float32)
        expected_train, expected_test, info = mod.spatial_split(coords, fold['seed'], .2, .05)
        assert np.array_equal(train, expected_train) and np.array_equal(test, expected_test)
        predictions = mod.make_predictions(x, logged, coords, train, test, 128, 10., 8, fold['seed'])
    else:
        from sklearn.decomposition import PCA
        from sklearn.linear_model import Ridge
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
        xx, yy = [], []
        legacy_inputs=[]
        for sid in fold['train_sample_ids']:
            scells, sx, sy, sg = cached_sample(root, fold['fixture'], sid,inputs_root)
            sub = membership[(membership['sample_id'] == sid) & (membership['partition'] == 'train')]
            train = sub['benchmark_row'].to_numpy(dtype=np.int64)
            xx.append(sx[train]); yy.append(np.log1p(target_columns(sy, sg, chosen)[train].astype(np.float32)).astype(np.float32))
            if fold['fixture']=='xenium_legacy':legacy_inputs.append((sy,sg,train))
        xx = np.concatenate(xx); yy = np.concatenate(yy)
        if fold['fixture']=='xenium_legacy':
            # Preserve the original cross evaluator's raw-column concatenation,
            # F-order log1p array, and final global train-row selection sequence.
            common_genes=common if material_group_guard else json_read(root/'targets/xenium_legacy_cross'/(fold['sample_id']+'.json'))['common_gene_symbols']
            joined=np.concatenate([target_columns(sy,sg,common_genes) for sy,sg,_ in legacy_inputs],axis=0)
            lookup={g:i for i,g in enumerate(common_genes)}
            joined=joined[:,[lookup[g] for g in chosen]]
            logged_joined=np.log1p(np.clip(joined,0,None)).astype(np.float32)
            offset=0;global_rows=[]
            for sy,sg,train in legacy_inputs:
                global_rows.append(train+offset);offset+=len(sy)
            yy=logged_joined[np.concatenate(global_rows)]
        pipeline = Pipeline([('scale',StandardScaler()),('pca',PCA(n_components=min(128,xx.shape[1],len(xx)-1),svd_solver='randomized',random_state=fold['seed'])),('ridge',Ridge(alpha=10.))])
        pipeline.fit(xx, yy)
        pred = pipeline.predict(x[test]).astype(np.float32)
        predictions = {'image': pred, 'train_mean': np.repeat(np.mean(yy,axis=0,keepdims=True),len(test),axis=0)}
    output.mkdir(parents=True, exist_ok=True)
    full_legacy_target=fold['fixture']=='xenium_legacy' and fold['regime']!='spatial'
    score_raw=np.asfortranarray(raw) if full_legacy_target else raw[test]
    score_logged=np.asfortranarray(logged) if full_legacy_target else logged[test]
    payload = dict(cell_ids=cells.iloc[test]['benchmark_cell_id'].to_numpy(dtype=str),
                   release_cell_ids=cells.iloc[test]['cell_id'].to_numpy(dtype=str),
                   source_h5ad_rows=cells.iloc[test]['source_h5ad_row'].to_numpy(dtype=np.int64),
                   gene_ids=np.asarray(chosen,dtype=str), y_true_raw=score_raw,y_true_log1p=score_logged,
                   **{'pred_'+k:v for k,v in predictions.items()})
    types=feature_types(root,fold['fixture'],fold['sample_id'],chosen)
    payload['target_feature_type']=np.asarray(types,dtype=str)
    payload['target_identifier_kind']=np.asarray('historical_gene_symbol_or_control_feature_name' if fold['fixture']=='xenium_legacy' else 'source_gene_id')
    payload['n_gene_expression_targets']=np.asarray(sum(t=='Gene Expression' for t in types))
    payload['n_control_targets']=np.asarray(sum(t!='Gene Expression' and any(w in t.lower() for w in ['control','codeword','blank']) for t in types))
    payload['primary_n_gene_expression_targets']=np.asarray(sum(t=='Gene Expression' for t in types[:50]))
    payload['primary_n_control_targets']=np.asarray(sum(t!='Gene Expression' and any(w in t.lower() for w in ['control','codeword','blank']) for t in types[:50]))
    payload['n_other_or_unknown_targets']=np.asarray(len(types)-int(payload['n_gene_expression_targets'])-int(payload['n_control_targets']))
    payload['primary_n_other_or_unknown_targets']=np.asarray(len(types[:50])-int(payload['primary_n_gene_expression_targets'])-int(payload['primary_n_control_targets']))
    np.savez_compressed(output / (fold_id + '.npz'), **payload)
    primary = min(fold.get('primary_genes',len(chosen)),len(chosen))
    values = {m:(legacy_primary_score(mod,score_raw,score_logged,p,chosen,primary) if not fold['modern_metrics'] else score(mod,raw[test,:primary],logged[test,:primary],p[:,:primary],chosen[:primary],True,fold.get('metric_dtype','float64' if fold['track']=='hd_cross' else 'float32'))) for m,p in predictions.items()}
    result = {'fold_id':fold_id,'metrics':values,'seed':fold['seed'],'trained':True,
              'material_group_guard':material_group_guard,'new_protocol':material_group_guard,
              'train_sample_ids':fold.get('train_sample_ids',[fold['sample_id']]),'gene_ids':chosen,
              'target_feature_type':types,
              **{k:int(payload[k]) for k in ['n_gene_expression_targets','n_control_targets','n_other_or_unknown_targets','primary_n_gene_expression_targets','primary_n_control_targets','primary_n_other_or_unknown_targets']},
              'input_type':'explicit prepared CONCH/RNA inputs; see raw_preparation/encoding_validation' if inputs_root is not None else 'frozen cached CONCH embeddings; image encoder not rerun'}
    import sklearn, scipy
    result['numeric_environment']={'numpy':np.__version__,'pandas':pd.__version__,'scikit_learn':sklearn.__version__,'scipy':scipy.__version__}
    (output / (fold_id + '.json')).write_text(json.dumps(result,indent=2,default=lambda x:float(x))+'\n')
    return result


def release_example(root, release_root, sample_id, output):
    """Read real release cell/raw RNA and the same native-pixel RGB crop rule."""
    sample_index = pd.read_parquet(release_root / 'metadata/samples.parquet')
    sample = sample_index.loc[sample_index.sample_id == sample_id].iloc[0]
    if str(sample.get('stain','')).lower() not in ['h&e','he','h\u0026e']:
        raise ValueError('The RGB benchmark example requires an H&E sample; IF/DAPI eligibility is a separate branch')
    cell_path = release_root / sample['cell_index_file']
    cells = pd.read_parquet(cell_path, columns=['cell_id','h5ad_row','x_px','y_px'])
    cells = cells[cells.h5ad_row.notna()]
    chosen = cells.iloc[len(cells)//2]
    from PIL import Image
    from image_reader import RGBReader
    reader = RGBReader(release_root / sample['image_file'], (sample['image_height_px'],sample['image_width_px']))
    try:
        image = reader.crop(chosen.x_px, chosen.y_px, 150)
        backend=reader.backend
    finally:
        reader.close()
    output.mkdir(parents=True,exist_ok=True)
    image.save(output / (sample_id+'_150px.png'))
    with h5py.File(release_root/sample['expression_file'],'r') as handle:
        row=int(chosen.h5ad_row);indptr=handle['X/indptr'][row:row+2];start,end=map(int,indptr)
        counts=handle['X/data'][start:end];cols=handle['X/indices'][start:end]
        index_name=handle['obs'].attrs.get('_index','_index')
        index_name=index_name.decode() if isinstance(index_name,bytes) else str(index_name)
        obs=handle['obs'][index_name][row];obs=obs.decode() if isinstance(obs,bytes) else str(obs)
        assert obs==chosen.cell_id
    result={'status':'passed','sample_id':sample_id,'cell_id':chosen.cell_id,'h5ad_row':row,'positive_raw_counts':int(np.sum(counts)),'n_nonzero_features':len(cols),'crop_native_px':150,'crop_shape':list(np.asarray(image).shape),'image_read_backend':backend,'crop_rule':'round centre (ties-to-even), left=round(x)-75, half-open150px; white padding'}
    (output/(sample_id+'_example.json')).write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    if not __debug__:
        raise SystemExit('Run without python -O; scientific integrity assertions must be enabled.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['list','demo','replay','refit','example','groups','prepare','encode','select-hvgs'])
    parser.add_argument('--benchmark-root',type=Path,help='Data-side benchmark asset directory (not this code checkout)')
    parser.add_argument('--release-root',type=Path)
    parser.add_argument('--track');parser.add_argument('--sample');parser.add_argument('--fold')
    parser.add_argument('--inputs-root',type=Path)
    parser.add_argument('--checkpoint',type=Path)
    parser.add_argument('--device',choices=['cpu','cuda'],default='cuda')
    parser.add_argument('--batch-size',type=int,default=8)
    parser.add_argument('--disable-cudnn',action='store_true')
    parser.add_argument('--material-group-guard',action='store_true')
    parser.add_argument('--allow-related-material',action='store_true')
    parser.add_argument('--limit',type=int);parser.add_argument('--output',type=Path,default=Path('benchmark-output'))
    args=parser.parse_args()
    if args.limit is not None and args.limit <= 0:
        parser.error('--limit must be positive')
    if args.batch_size <= 0:
        parser.error('--batch-size must be positive')
    if args.command=='list':
        print(json.dumps(json_read(Path(__file__).resolve().parent/'tracks.json'),indent=2));return
    if args.command=='demo':
        from demo import run
        run(args.output);return
    if args.release_root is None and os.environ.get('SMMC_RELEASE_ROOT'):
        args.release_root=Path(os.environ['SMMC_RELEASE_ROOT'])
    assets=args.benchmark_root
    if assets is None and os.environ.get('SMMC_BENCHMARK_ROOT'):
        assets=Path(os.environ['SMMC_BENCHMARK_ROOT'])
    if assets is None and args.release_root is not None:
        assets=args.release_root/'benchmark'
    if assets is None:
        parser.error('Exact reproduction requires the separately stored frozen asset package. Set SMMC_RELEASE_ROOT, SMMC_BENCHMARK_ROOT, or --benchmark-root. See benchmark/README.md. Use demo for an offline smoke.')
    root=assets.expanduser().resolve()
    if not (root/'folds.json').is_file() or not (root/'fixtures').is_dir():
        parser.error('Incomplete benchmark assets: expected folds.json and fixtures/ in '+str(root)+'. This Git checkout contains code/plans, not the large data-side assets.')
    if args.release_root is not None:
        args.release_root=args.release_root.expanduser().resolve()
    destination=args.output.expanduser().resolve()
    protected=[root]+([args.release_root] if args.release_root is not None else [])
    if any(destination==p or p in destination.parents for p in protected):
        parser.error('Outputs must be outside the immutable data/benchmark package')
    args.output=destination
    if args.command=='refit' and args.fold is None:
        parser.error('refit requires --fold (see benchmark/plans/folds.json)')
    if args.command=='replay':result=replay(root,args.track,args.sample,args.limit)
    elif args.command=='refit':result=refit(root,args.fold,args.output,args.material_group_guard,args.allow_related_material,args.inputs_root)
    elif args.command=='groups':result=group_audit(root,args.output)
    elif args.command in ['prepare','encode','select-hvgs']:
        import pipeline
        inputs=args.inputs_root or args.output
        if args.command=='select-hvgs':result=pipeline.select_hvgs(root,args.track,inputs,args.sample,args.output)
        else:
            if args.release_root is None or args.track is None or args.sample is None:parser.error('prepare/encode need --release-root, --track and --sample')
            chosen=pipeline.samples(root,args.track) if args.sample=='all' else [args.sample]
            if args.command=='prepare':rr=[pipeline.prepare(root,args.release_root,args.track,s,inputs,args.limit) for s in chosen]
            else:
                if args.checkpoint is None:parser.error('encode requires an external --checkpoint with the declared SHA256')
                rr=[pipeline.encode(root,args.release_root,args.track,s,inputs,args.checkpoint,args.device,args.batch_size,args.disable_cudnn) for s in chosen]
            result={'status':'passed','samples':len(rr),'results':rr}
    else:
        if args.release_root is None or args.sample is None:parser.error('example needs --release-root and --sample (H&E)')
        result=release_example(root,args.release_root.resolve(),args.sample,args.output)
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/(args.command+'_validation.json')).write_text(json.dumps(result,indent=2,allow_nan=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},default=lambda x:float(x)),flush=True)


if __name__=='__main__':main()
