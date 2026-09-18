import json
import sys
from pathlib import Path

nb_path = Path(sys.argv[1])
out_dir = Path(sys.argv[2])
out_dir.mkdir(parents=True, exist_ok=True)
nb = json.loads(nb_path.read_text())
ns = {"__name__": "__main__"}
for i, cell in enumerate(nb["cells"]):
    if cell.get("cell_type") != "code":
        continue
    src = "".join(cell.get("source", []))
    if not src.strip():
        continue
    print(f"--- executing cell {i} ---", flush=True)
    exec(compile(src, f"{nb_path}:cell{i}", "exec"), ns)
    if i == 10 and "adata" in ns:
        cols = [c for c in ["sample", "time", "time_group", "total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_rb", "doublet_score", "leiden"] if c in ns["adata"].obs]
        ns["adata"].obs[cols].to_csv(out_dir / "pre_first_manual_exclusion_obs.csv")
    if i == 16 and "adata_c" in ns:
        cols = [c for c in ["sample", "time", "time_group", "total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_rb", "doublet_score", "leiden"] if c in ns["adata_c"].obs]
        ns["adata_c"].obs[cols].to_csv(out_dir / "pre_second_manual_exclusion_obs.csv")
    if i >= 16:
        break
