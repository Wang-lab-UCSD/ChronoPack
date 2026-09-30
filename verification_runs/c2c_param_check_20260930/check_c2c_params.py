"""Run cp.tl.c2c.infer with the notebook's arguments and with reduced variants; compare outputs."""
import json, os
import pandas as pd
import chronopack as cp

os.chdir("/new-stg/home/alex/tmp/ChronoPack/c2c")
nb = json.load(open("CellCommunication.ipynb"))
ns = {"__name__": "__main__"}
for i in (1, 2):  # code cells: imports + configuration (paths and parameters)
    exec("".join([c for c in nb["cells"] if c["cell_type"] == "code"][i - 1]["source"]), ns)
import scanpy as sc
adata = sc.read_h5ad(ns["ADATA_PATH"])
print("adata", adata.shape, "| response genes present:", [g for g in ns["RESPONSE_GENES"] if g in adata.var_names])

base_args = dict(
    adata=adata, lr_resource=ns["CELLCHAT_PATH"], association=ns["ASSO_BOTH_PATH"],
    sender=ns["SENDER_CELL_TYPE"], receiver=ns["RECEIVER_CELL_TYPE"], response_genes=ns["RESPONSE_GENES"],
    cell_type_key=ns["CELL_TYPE_KEY"], counts_layer=ns["COUNTS_LAYER"], expression_layer=ns["EXPR_LAYER"],
    min_detect=ns["MIN_DETECT"], causal_window_min=ns["CAUSAL_WINDOW_MIN"], association_p_col="pvalue_2",
    prefer_fits=ns["PREFER_FITS_FOR_SCORING"], tradeseq_export_dir=ns["TS_EXPORT_GRID_DIR"],
)
base = cp.tl.c2c.infer(**base_args)
print("baseline: pairs", base.interactions.shape, "| ligands", base.ligand_scores.shape, "| profiles", base.temporal_profiles.shape)

minimal = {k: base_args[k] for k in ["adata", "lr_resource", "association", "sender", "receiver", "response_genes", "min_detect"]}
variants = {
    'response_genes=["NFKBIA"]': {**base_args, "response_genes": ["NFKBIA"]},
    "no tradeseq_export_dir": {k: v for k, v in base_args.items() if k != "tradeseq_export_dir"},
    "default-valued args dropped": {k: v for k, v in base_args.items() if k not in
        ["cell_type_key", "counts_layer", "expression_layer", "causal_window_min", "association_p_col", "prefer_fits"]},
    "min_detect default (0.0)": {k: v for k, v in base_args.items() if k != "min_detect"},
    "minimal call": minimal,
}
for name, args in variants.items():
    res = cp.tl.c2c.infer(**args)
    same = {t: getattr(res, t).equals(getattr(base, t)) for t in ["interactions", "ligand_scores", "candidate_interactions"]}
    same["temporal_profiles"] = res.temporal_profiles.drop(columns="receiver_response").equals(base.temporal_profiles.drop(columns="receiver_response"))
    same["receiver_response"] = res.temporal_profiles["receiver_response"].equals(base.temporal_profiles["receiver_response"]) if len(res.temporal_profiles) == len(base.temporal_profiles) else False
    print(f"{name:30s} pairs {len(res.interactions):4d} | identical: {same}")
print("defaults in signature:", {k: v.default for k, v in __import__("inspect").signature(cp.tl.c2c.infer).parameters.items() if v.default is not __import__("inspect").Parameter.empty})
print("C2C_CHECK_OK")
