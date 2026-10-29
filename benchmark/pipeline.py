"""Modern raw-release -> counts -> CONCH -> training-HVG pipeline.

The exported module uses only --release-root, external checkpoint, and package-
relative assets. Frozen sampling identities are reused to reproduce the reported
experiment; no new outcome-dependent sample selection occurs.
"""
from __future__ import annotations
import contextlib
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def js(p):return json.loads(Path(p).read_text())
def dump(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
def decode(x):return [v.decode() if isinstance(v,bytes) else str(v) for v in x]
def fixture_path(inputs_root,track,sid):return inputs_root/'fixtures'/track/sid
def samples(root,track):
    if track not in ['hd18','xenium5']:raise ValueError('Raw prepare/encode/HVG is defined for modern hd18 and xenium5; legacy caches preserve a separate historical input protocol')
    return sorted({f['sample_id'] for f in js(root/'folds.json') if f['track']==track})
def source_sample(release_root,sid):
    frame=pd.read_parquet(release_root/'metadata/samples.parquet');selected=frame[frame.sample_id==sid]
    if len(selected)!=1:raise ValueError('Release sample index is not unique for '+sid)
    return selected.iloc[0]


def prepare(root,release_root,track,sid,inputs_root,limit=None):
    if sid not in samples(root,track):raise ValueError('Unknown modern acquisition '+sid)
    reference=pd.read_parquet(root/'fixtures'/track/sid/'cells.parquet')
    if limit:
        # Spread the smoke across the original seed0 test region, allowing direct
        # prediction parity against the archived edge-test predictions.
        f=next(f for f in js(root/'folds.json') if f['fold_id']==f'{track}_{sid}_seed0')
        members=pd.read_parquet(root/f['membership_file'])
        candidates=members[members.partition=='test'].benchmark_row.to_numpy(dtype=np.int64)
        selected=candidates[np.linspace(0,len(candidates)-1,min(limit,len(candidates)),dtype=int)]
    else:selected=np.arange(len(reference),dtype=np.int64)
    cells=reference.iloc[selected].copy().reset_index(drop=True)
    cells['frozen_benchmark_row']=selected
    cells['benchmark_row']=np.arange(len(cells),dtype=np.int64)
    gene_axis=pd.read_parquet(root/'fixtures'/track/sid/'genes.parquet')
    source=source_sample(release_root,sid)
    current_genes=pd.read_parquet(release_root/source['genes_file'])
    source_columns=gene_axis.source_gene_index.to_numpy(dtype=np.int64)
    assert gene_axis.gene_id.astype(str).tolist()==current_genes.iloc[source_columns].gene_id.astype(str).tolist()
    assert set(current_genes.iloc[source_columns].feature_type)=={'Gene Expression'}
    data=[];indices=[];ptr=[0];pixel=[];physical=[]
    with h5py.File(release_root/source['expression_file'],'r') as handle:
        obs_key=handle['obs'].attrs.get('_index','_index');obs_key=obs_key.decode() if isinstance(obs_key,bytes) else str(obs_key)
        shape=tuple(handle['X'].attrs['shape'])
        assert max(source_columns)<shape[1]
        assert np.issubdtype(handle['X/data'].dtype,np.integer)
        for record in cells.itertuples(index=False):
            row=int(record.source_h5ad_row);cell=handle['obs'][obs_key][row];cell=cell.decode() if isinstance(cell,bytes) else str(cell)
            assert cell==record.cell_id
            a,b=map(int,handle['X/indptr'][row:row+2]);vv=handle['X/data'][a:b];ii=handle['X/indices'][a:b]
            assert np.all(vv>=0)
            data.append(vv.astype(np.int64));indices.append(ii.astype(np.int64));ptr.append(ptr[-1]+len(ii))
            pixel.append(handle['obsm/spatial'][row]);physical.append(handle['obsm/spatial_um'][row])
        raw=sparse.csr_matrix((np.concatenate(data),np.concatenate(indices),np.asarray(ptr,dtype=np.int64)),shape=(len(cells),shape[1]))[:,source_columns].tocsr()
    raw.sum_duplicates();raw.sort_indices()
    cache=sparse.load_npz(root/'fixtures'/track/sid/'counts.npz')[selected]
    assert (raw!=cache).nnz==0,'Fresh release CSR differs from the archived sampled integer RNA'
    pixel=np.asarray(pixel);physical=np.asarray(physical)
    assert np.allclose(pixel,cells[['x_px','y_px']].to_numpy(),rtol=0,atol=1e-7)
    assert np.allclose(physical.astype(np.float32),cells[['coord_x','coord_y']].to_numpy(dtype=np.float32),rtol=0,atol=1e-4)
    out=fixture_path(inputs_root,track,sid);out.mkdir(parents=True,exist_ok=True)
    sparse.save_npz(out/'counts.npz',raw)
    cells.to_parquet(out/'cells.parquet',index=False,compression='zstd');gene_axis.to_parquet(out/'genes.parquet',index=False,compression='zstd')
    report={'status':'passed','sample_id':sid,'track':track,'observations':len(cells),'genes':raw.shape[1],'raw_count_sum':int(raw.sum()),
            'RNA_source':'versioned release raw integer CSR; exact frozen cell/gene ID axes',
            'matches_archived_counts_exactly':True,'max_pixel_coordinate_error':float(np.max(np.abs(pixel-cells[['x_px','y_px']].to_numpy()))),
            'full_frozen_sample':not bool(limit),'frozen_benchmark_rows':selected.tolist(),
            'source_expression_file':source['expression_file'],'source_genes_file':source['genes_file'],
            'target_feature_type':'Gene Expression only; no cross-panel zero padding or normalization'}
    dump(out/'raw_preparation.json',report);return report


def allocated_gpu(torch):
    if not os.environ.get('SLURM_JOB_ID'):return None
    assert torch.cuda.device_count()==1,'The smoke must use one allocated CUDA-visible device'
    import subprocess
    p=torch.cuda.get_device_properties(0);uuid=str(p.uuid)
    if not uuid.startswith('GPU-'):uuid='GPU-'+uuid
    q=subprocess.run(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,process_name,used_gpu_memory','--format=csv,noheader,nounits'],capture_output=True,text=True,check=True)
    active=[]
    for line in q.stdout.splitlines():
        v=[x.strip() for x in line.split(',')]
        if len(v)>=2 and v[0]==uuid and int(v[1])!=os.getpid():active.append(v)
    if active:raise RuntimeError('Assigned GPU already has another compute process: '+str(active))
    return {'uuid':uuid,'name':p.name,'slurm_job_id':os.environ['SLURM_JOB_ID'],'node':os.uname().nodename}


def encoder(root,checkpoint,device,disable_cudnn=False):
    import torch
    spec=js(root/'models/conch_v1_checkpoint.json')
    h=hashlib.sha256()
    with checkpoint.open('rb') as stream:
        for part in iter(lambda:stream.read(16*1024*1024),b''):h.update(part)
    assert checkpoint.stat().st_size==spec['size_bytes'] and h.hexdigest()==spec['sha256']
    sys.path.insert(0,str(root/'vendor'))
    from conch.open_clip_custom import factory
    from conch.open_clip_custom.coca_model import _build_vision_tower
    config=js(root/'vendor/conch/open_clip_custom/model_configs/conch_ViT-B-16.json');config.pop('custom_text',None)
    state=torch.load(checkpoint,map_location='cpu',weights_only=True,mmap=True)
    if 'state_dict' in state:state=state['state_dict']
    if next(iter(state)).startswith('module.'):state={k[7:]:v for k,v in state.items()}
    model=_build_vision_tower(config['embed_dim'],config['vision_cfg'],embed_dim_caption=config['embed_dim_caption'])
    model.load_state_dict({k[len('visual.'):]:v for k,v in state.items() if k.startswith('visual.')},strict=True)
    model.eval()
    for p in model.parameters():p.requires_grad_(False)
    if disable_cudnn:torch.backends.cudnn.enabled=False
    model.to(device)
    transform=factory.image_transform(model.image_size,mean=factory.OPENAI_DATASET_MEAN,std=factory.OPENAI_DATASET_STD)
    return model,transform,spec['sha256']


def encode(root,release_root,track,sid,inputs_root,checkpoint,device='cuda',batch_size=8,disable_cudnn=False):
    os.environ.setdefault('HF_HOME',str(inputs_root/'cache/huggingface'))
    os.environ.setdefault('TORCH_HOME',str(inputs_root/'cache/torch'))
    import torch
    from image_reader import RGBReader
    if device=='cuda' and not torch.cuda.is_available():raise ValueError('No CUDA device; choose --device cpu or allocate one GPU')
    gpu=allocated_gpu(torch) if device=='cuda' else None
    torch.manual_seed(0);torch.set_num_threads(max(1,int(os.environ.get('OMP_NUM_THREADS','4'))))
    directory=fixture_path(inputs_root,track,sid)
    cells=pd.read_parquet(directory/'cells.parquet');source=source_sample(release_root,sid)
    if str(source['stain']).lower()!='h&e':raise ValueError('Modern CONCH requires the H&E branch')
    model,transform,sha=encoder(root,checkpoint,device,disable_cudnn)
    reader=RGBReader(release_root/source['image_file'],(source['image_height_px'],source['image_width_px']))
    vectors=[];started=time.monotonic()
    try:
        with torch.inference_mode():
            for first in range(0,len(cells),batch_size):
                part=cells.iloc[first:first+batch_size]
                imgs=torch.stack([transform(reader.crop(row.x_px,row.y_px,150)) for row in part.itertuples(index=False)])
                with torch.autocast('cuda',dtype=torch.float16) if device=='cuda' else contextlib.nullcontext():
                    xx=model.forward_no_head(imgs.to(device),normalize=False)
                vectors.append(xx.float().cpu().numpy())
    finally:reader.close()
    x=np.concatenate(vectors).astype(np.float32);assert x.shape==(len(cells),512) and np.isfinite(x).all()
    with h5py.File(root/'fixtures'/track/sid/'features.h5','r') as f:
        all_reference=np.asarray(f['embeddings'][:],dtype=np.float32)
    selected=cells.frozen_benchmark_row.to_numpy(dtype=np.int64);reference=all_reference[selected]
    cosine=np.sum(x.astype(float)*reference,axis=1)/(np.linalg.norm(x,axis=1)*np.linalg.norm(reference,axis=1))
    rel=np.linalg.norm(x.astype(float)-reference,axis=1)/np.maximum(np.linalg.norm(reference,axis=1),1e-15)
    assert np.min(cosine)>=.999 and np.max(rel)<=.01,'New raw-image CONCH features do not match the frozen feature protocol'
    with h5py.File(directory/'features.h5','w') as f:
        f.create_dataset('embeddings',data=x,compression='lzf')
        f.create_dataset('cell_ids',data=np.asarray(cells.benchmark_cell_id.astype(str),dtype=object),dtype=h5py.string_dtype('utf-8'))
        f.attrs.update(checkpoint_sha256=sha,crop_native_px=150,preprojection=True,feature_normalized=False,feature_dimensions=512,device=device)
    # Compare a freshly fitted seed0 image predictor on fresh versus archived
    # features, and compare their matching test rows to the saved predictions.
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.linear_model import Ridge
    fold=next(f for f in js(root/'folds.json') if f['fold_id']==f'{track}_{sid}_seed0')
    membership=pd.read_parquet(root/fold['membership_file']);train=membership[membership.partition=='train'].benchmark_row.to_numpy(dtype=np.int64)
    cache_genes=pd.read_parquet(root/'fixtures'/track/sid/'genes.parquet').gene_id.astype(str).tolist();lookup={g:i for i,g in enumerate(cache_genes)}
    counts=sparse.load_npz(root/'fixtures'/track/sid/'counts.npz')[:,[lookup[g] for g in fold['gene_ids']]].toarray()
    y=np.log1p(counts.astype(np.float32)).astype(np.float32)
    pipe=Pipeline([('scale',StandardScaler()),('pca',PCA(n_components=min(128,len(train)-1,all_reference.shape[1]),svd_solver='randomized',random_state=0)),('ridge',Ridge(alpha=10.))])
    pipe.fit(all_reference[train],y[train]);new_pred=pipe.predict(x).astype(np.float32);old_pred=pipe.predict(reference).astype(np.float32)
    pred_diff=float(np.sqrt(np.mean((new_pred-old_pred)**2)))
    assert pred_diff<=.02,'Raw image to gene prediction differs excessively'
    record=next(r for r in js(root/'prediction_index.json') if r['track']==track and r['sample_id']==sid and r['seed']==0)
    with np.load(root/record['path'],allow_pickle=False) as packet:
        positions={int(row):i for i,row in enumerate(packet['test_indices'])}
        overlap=[i for i,row in enumerate(selected) if int(row) in positions]
        archived=packet['pred_image'][[positions[int(selected[i])] for i in overlap]]
        archive_error=float(np.max(np.abs(old_pred[overlap]-archived)))
        fresh_archive_rmse=float(np.sqrt(np.mean((new_pred[overlap]-archived)**2)))
    assert archive_error<=1e-3 and fresh_archive_rmse<=.02
    report={'status':'passed','track':track,'sample_id':sid,'observations':len(cells),'checkpoint_sha256':sha,'gpu':gpu,'device':device,'disable_cudnn':disable_cudnn,
            'feature_dimensions':512,'minimum_feature_cosine':float(np.min(cosine)),'maximum_relative_feature_L2':float(np.max(rel)),
            'feature_max_absolute_error':float(np.max(np.abs(x-reference))), 'predictor_old_vs_fresh_feature_RMSE':pred_diff,
            'refitted_original_features_vs_archived_prediction_max_abs_error':archive_error,'fresh_image_features_vs_archived_prediction_RMSE':fresh_archive_rmse,
            'archived_test_rows_compared':len(overlap),'elapsed_encoding_seconds':time.monotonic()-started,
            'scope':'Raw release image -> native150 crop -> official448 transform -> strict visual checkpoint ->512 features -> frozen seed0 ridge target; approximate cross-precision feature parity, not a full-library rerun'}
    dump(directory/'encoding_validation.json',report);return report


def select_hvgs(root,track,inputs_root,sample_id=None,output_root=None):
    # Frozen group membership is defined by measured-panel metadata, not outcomes.
    samples(root,track)
    folds=[f for f in js(root/'folds.json') if f['track']==track and (sample_id in [None,'all'] or f['sample_id']==sample_id)]
    if not folds:raise ValueError('No modern track')
    protocol_spec=importlib.util.spec_from_file_location('hvg_spatial_protocol',root/'frozen/evaluate_spatial.py')
    module=importlib.util.module_from_spec(protocol_spec);protocol_spec.loader.exec_module(module)
    planned={}
    for f in folds:
        target=js(root/f['target_file']);members=target['source_samples'];key=(tuple(members),f['seed'])
        planned.setdefault(key,[]).append(f)
    results=[]
    for (members,seed),ff in sorted(planned.items()):
        axes={};common=None
        for sid in members:
            p=fixture_path(inputs_root,track,sid);cells=pd.read_parquet(p/'cells.parquet');genes=pd.read_parquet(p/'genes.parquet')
            if len(cells)!=len(pd.read_parquet(root/'fixtures'/track/sid/'cells.parquet')):raise ValueError('HVG selection requires every frozen sampled cell; --limit is only an encoding/preparation smoke')
            axes[sid]=(cells,genes);available=set(genes.gene_id.astype(str));common=available if common is None else common&available
        universe=sorted(common)
        if len(universe)<50:raise ValueError('Group has fewer than50 shared measured genes; retain the frozen Skin exact-panel group amendment')
        variance=np.zeros(len(universe),dtype=np.float64)
        for sid,(cells,genes) in axes.items():
            lookup={g:i for i,g in enumerate(genes.gene_id.astype(str))};raw=sparse.load_npz(fixture_path(inputs_root,track,sid)/'counts.npz')[:,[lookup[g] for g in universe]]
            train,test,info=module.spatial_split(cells[['coord_x','coord_y']].to_numpy(dtype=np.float32),seed,.2,.05)
            log=raw[train].astype(np.float64);log.data=np.log1p(log.data);mu=np.asarray(log.mean(axis=0)).ravel();second=log.copy();second.data**=2;second=np.asarray(second.mean(axis=0)).ravel();vv=np.maximum(second-mu**2,0);vv[vv<=64*np.finfo(float).eps*np.maximum(second,1)]=0
            variance+=vv/len(members)
        valid=[i for i in range(len(universe)) if variance[i]>0];ordered=sorted(valid,key=lambda i:(-variance[i],universe[i]))[:50]
        chosen=[universe[i] for i in ordered]
        if len(chosen)!=50:raise ValueError('Fewer than50 variable measured genes')
        for f in ff:
            assert chosen==f['gene_ids'],'Fresh training-only HVGs differ from the frozen target list'
            record={'status':'passed','track':track,'sample_id':f['sample_id'],'seed':seed,'gene_ids':chosen,'source_samples':list(members),'selection':'equal-sample mean population log1p variance in training regions only; gene-ID ascending tie break','matches_frozen_targets':True}
            dump((output_root or inputs_root)/'hvgs'/track/f['sample_id']/(f'fold_{seed}.json'),record);results.append(record)
    report={'status':'passed','track':track,'groups_times_folds':len(planned),'sample_folds':len(results),'source_expression':'explicit input root (fresh prepare outputs or frozen sampled counts)','targets_match_frozen_reference':True}
    dump((output_root or inputs_root)/'hvg_selection_validation.json',report);return report
