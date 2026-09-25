from pathlib import Path
import hashlib, json, h5py
RUN=Path('/new-stg/home/alex/tmp/ChronoPack/verification_runs/manuscript_full_pipeline_20260919_165539')
REF=Path('/new-stg/home/alex/projects/chrono_seq/analyses')
ALT={
'puro_pre_full': Path('/new-stg/home/alex/tmp/ChronoPack/verification_runs/manuscript_reproduction_20260919/preprocess_puromycin_numba_align/output/adata_full.chk0.h5ad'),
'puro_pre_hvg': Path('/new-stg/home/alex/tmp/ChronoPack/verification_runs/manuscript_reproduction_20260919/preprocess_puromycin_numba_align/output/adata_hvg.chk0.h5ad'),
}
ITEMS=[
('k562_pre_full','k562-tnfa.hg38.rep/output/adata_full.1e4-variant.chk0.h5ad'),
('k562_pre_hvg','k562-tnfa.hg38.rep/output/adata_hvg.1e4-variant.chk0.h5ad'),
('k562_traj_full','k562-tnfa.hg38.rep/output/adata_full.1e4-variant.chk1.h5ad'),
('k562_traj_hvg','k562-tnfa.hg38.rep/output/adata_hvg.1e4-variant.chk1.h5ad'),
('puro_pre_full','k562-tnfa-puromycin.hg38/output/adata_full.chk0.h5ad'),
('puro_pre_hvg','k562-tnfa-puromycin.hg38/output/adata_hvg.chk0.h5ad'),
('puro_traj_full','k562-tnfa-puromycin.hg38/output/adata_full.h5ad'),
('puro_traj_hvg','k562-tnfa-puromycin.hg38/output/adata_hvg.h5ad'),
('lps_pre_full','thp1-k562-lps.hg38.rep/output/adata_full.chk0.h5ad'),
('lps_pre_hvg','thp1-k562-lps.hg38.rep/output/adata_hvg.chk0.h5ad'),
('lps_traj_full','thp1-k562-lps.hg38.rep/output/adata_full.chk1.h5ad'),
('lps_traj_hvg','thp1-k562-lps.hg38.rep/output/adata_hvg.chk1.h5ad'),
]
ROOTS=['X','layers','obsm','obsp']

def hash_ds(ds):
    h=hashlib.sha256(); h.update(str(ds.shape).encode()); h.update(str(ds.dtype).encode())
    if ds.shape == ():
        v=ds[()]
        h.update(v.tobytes() if hasattr(v,'tobytes') else str(v).encode())
        return h.hexdigest()
    n=ds.shape[0] if ds.shape else 1
    step=min(max(n,1), 256)
    for i in range(0,n,step):
        arr=ds[i:min(i+step,n)]
        h.update(arr.tobytes() if hasattr(arr,'tobytes') else str(arr).encode())
    return h.hexdigest()

def collect(path):
    out={}
    with h5py.File(path,'r') as f:
        def visit(prefix):
            if prefix not in f: return
            obj=f[prefix]
            if isinstance(obj,h5py.Dataset):
                out[prefix]=hash_ds(obj); return
            def cb(name,obj):
                if isinstance(obj,h5py.Dataset):
                    out[prefix+'/'+name]=hash_ds(obj)
            obj.visititems(cb)
        for r in ROOTS: visit(r)
    return out
res={}
for name,rel in ITEMS:
    rp=RUN/rel; fp=ALT.get(name, REF/rel)
    item={'run_path':str(rp),'ref_path':str(fp),'run_exists':rp.exists(),'ref_exists':fp.exists()}
    if rp.exists() and fp.exists():
        rh=collect(rp); fh=collect(fp)
        common=sorted(set(rh)&set(fh))
        diff=[k for k in common if rh[k]!=fh[k]]
        item.update(n_common=len(common), n_diff=len(diff), diff_sample=diff[:30], only_run=sorted(set(rh)-set(fh))[:20], only_ref=sorted(set(fh)-set(rh))[:20])
    res[name]=item
out=RUN/'intermediate_h5ad_value_hash_comparison.json'
out.write_text(json.dumps(res,indent=2,sort_keys=True))
print(out)
for k,v in res.items():
    print(k, 'diff', v.get('n_diff'), 'common', v.get('n_common'), 'sample', v.get('diff_sample'))
