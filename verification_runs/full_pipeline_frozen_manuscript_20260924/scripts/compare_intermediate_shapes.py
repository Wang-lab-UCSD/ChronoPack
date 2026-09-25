from pathlib import Path
import json, gzip, hashlib
import anndata as ad
import pandas as pd

RUN=Path('/new-stg/home/alex/tmp/ChronoPack/verification_runs/manuscript_full_pipeline_20260919_165539')
REF=Path('/new-stg/home/alex/projects/chrono_seq/analyses')

h5ads = [
('k562_pre_full','k562-tnfa.hg38.rep/output/adata_full.1e4-variant.chk0.h5ad'),
('k562_pre_hvg','k562-tnfa.hg38.rep/output/adata_hvg.1e4-variant.chk0.h5ad'),
('k562_traj_full','k562-tnfa.hg38.rep/output/adata_full.1e4-variant.chk1.h5ad'),
('k562_traj_hvg','k562-tnfa.hg38.rep/output/adata_hvg.1e4-variant.chk1.h5ad'),
('k562_prefilter','k562-tnfa.hg38.rep/output/adata_prefiltered.1e4-variant.fixed.h5ad'),
('puro_traj_full','k562-tnfa-puromycin.hg38/output/adata_full.h5ad'),
('puro_traj_hvg','k562-tnfa-puromycin.hg38/output/adata_hvg.h5ad'),
('lps_pre_full','thp1-k562-lps.hg38.rep/output/adata_full.chk0.h5ad'),
('lps_pre_hvg','thp1-k562-lps.hg38.rep/output/adata_hvg.chk0.h5ad'),
('lps_traj_full','thp1-k562-lps.hg38.rep/output/adata_full.chk1.h5ad'),
('lps_traj_hvg','thp1-k562-lps.hg38.rep/output/adata_hvg.chk1.h5ad'),
('lps_prefilter_k562','thp1-k562-lps.hg38.rep/output/adata_k562_prefiltered.fixed.h5ad'),
('lps_prefilter_thp1','thp1-k562-lps.hg38.rep/output/adata_thp1_prefiltered.fixed.h5ad'),
]
# puromycin chk0 refs are absent from /projects output, compare against previously verified manuscript run if present
alt_ref = {
 'puro_pre_full': Path('/new-stg/home/alex/tmp/ChronoPack/verification_runs/manuscript_reproduction_20260919/preprocess_puromycin_numba_align/output/adata_full.chk0.h5ad'),
 'puro_pre_hvg': Path('/new-stg/home/alex/tmp/ChronoPack/verification_runs/manuscript_reproduction_20260919/preprocess_puromycin_numba_align/output/adata_hvg.chk0.h5ad'),
 'puro_prefilter': None,
}
h5ads += [('puro_pre_full','k562-tnfa-puromycin.hg38/output/adata_full.chk0.h5ad'),('puro_pre_hvg','k562-tnfa-puromycin.hg38/output/adata_hvg.chk0.h5ad'),('puro_prefilter','k562-tnfa-puromycin.hg38/output/adata_prefiltered.1e4-variant.fixed.h5ad')]

def fp_index(vals):
    h=hashlib.sha256()
    for v in vals:
        h.update(str(v).encode()); h.update(b'\0')
    return h.hexdigest()

def h5_summary(path):
    out={'exists': path is not None and path.exists(), 'path': str(path) if path is not None else None}
    if not out['exists']: return out
    a=ad.read_h5ad(path, backed='r')
    out.update(shape=list(a.shape), obs_hash=fp_index(a.obs_names), var_hash=fp_index(a.var_names))
    for col in ['time','leiden','leiden_unintegrated','traj_clusters','predicted_doublet','cell_type','sample']:
        if col in a.obs:
            out[f'obs_{col}_hash']=fp_index(a.obs[col].astype(str).tolist())
            out[f'obs_{col}_counts']=a.obs[col].astype(str).value_counts(dropna=False).sort_index().to_dict()
    for key in ['ot_graph','t_pruned','W_terminal_cosine_topk','terminal_cell_dist']:
        if key in a.obsp:
            out[f'obsp_{key}_shape']=list(a.obsp[key].shape)
            out[f'obsp_{key}_nnz']=int(getattr(a.obsp[key], 'nnz', -1))
    a.file.close()
    return out

res={'h5ad':{}, 'tradeSeq_inputs':{}}
for name, rel in h5ads:
    rp=RUN/rel
    fp=alt_ref.get(name) if name in alt_ref else REF/rel
    rs=h5_summary(rp); fs=h5_summary(fp)
    keys=set(rs)|set(fs)
    diffs={k:[rs.get(k),fs.get(k)] for k in sorted(keys) if rs.get(k)!=fs.get(k)}
    res['h5ad'][name]={'run':rs,'ref':fs,'diff_keys':list(diffs.keys()),'diffs':diffs}

csvs=[
('k562_genes','k562-tnfa.hg38.rep/intmd/sc-K562-1e4_to_tradeSeq.fixed/subset_genes.csv'),
('k562_genes_full','k562-tnfa.hg38.rep/intmd/sc-K562-1e4_to_tradeSeq.fixed/subset_genes.full.csv'),
('puro_genes','k562-tnfa-puromycin.hg38/intmd/sc-K562-1e6-to-tradeSeq.fixed/subset_genes.csv'),
('puro_genes_full','k562-tnfa-puromycin.hg38/intmd/sc-K562-1e6-to-tradeSeq.fixed/subset_genes.full.csv'),
('lps_k562_genes','thp1-k562-lps.hg38.rep/intmd/sc-K562_to_tradeSeq.fixed/subset_genes.csv'),
('lps_k562_genes_full','thp1-k562-lps.hg38.rep/intmd/sc-K562_to_tradeSeq.fixed/subset_genes.full.csv'),
('lps_thp1_genes','thp1-k562-lps.hg38.rep/intmd/sc-THP1_to_tradeSeq.fixed/subset_genes.csv'),
('lps_thp1_genes_full','thp1-k562-lps.hg38.rep/intmd/sc-THP1_to_tradeSeq.fixed/subset_genes.full.csv'),
]
for name, rel in csvs:
    rp, fp=RUN/rel, REF/rel
    item={'run_exists':rp.exists(),'ref_exists':fp.exists(),'run_path':str(rp),'ref_path':str(fp)}
    if rp.exists() and fp.exists():
        r=pd.read_csv(rp, header=None).iloc[:,0].astype(str).tolist()
        f=pd.read_csv(fp, header=None).iloc[:,0].astype(str).tolist()
        item.update(n_run=len(r), n_ref=len(f), ordered_equal=(r==f), n_common=len(set(r)&set(f)), only_run_head=sorted(set(r)-set(f))[:10], only_ref_head=sorted(set(f)-set(r))[:10])
    res['tradeSeq_inputs'][name]=item

out=RUN/'intermediate_reference_shape_index_comparison.json'
out.write_text(json.dumps(res, indent=2, sort_keys=True))
print(out)
print('\nH5AD concise')
for k,v in res['h5ad'].items():
    r=v['run']; f=v['ref']
    print(k, 'run', r.get('shape'), 'ref', f.get('shape'), 'obs', r.get('obs_hash')==f.get('obs_hash'), 'var', r.get('var_hash')==f.get('var_hash'), 'diffkeys', v['diff_keys'][:8])
print('\nGenes concise')
for k,v in res['tradeSeq_inputs'].items():
    print(k, 'run/ref', v.get('n_run'), v.get('n_ref'), 'ordered', v.get('ordered_equal'), 'common', v.get('n_common'), 'exists', v['run_exists'], v['ref_exists'])
