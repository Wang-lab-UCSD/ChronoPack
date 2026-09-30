import hashlib
from pathlib import Path
import pandas as pd

NEW = Path(__file__).resolve().parent / "outputs"
REF = Path("/new-stg/home/alex/tmp/ChronoPack/c2c/outputs")
for name in ["ligand_scores", "interactions", "candidate_interactions", "temporal_profiles"]:
    a = pd.read_csv(NEW / "chronopack_c2c" / f"chronopack_c2c.{name}.csv")
    b = pd.read_csv(REF / "chronopack_c2c" / f"chronopack_c2c.{name}.csv")
    dropped = [c for c in b.columns if c not in a.columns]
    b = b[[c for c in b.columns if c in a.columns]]
    print(f"{name:24s} shape {a.shape} vs {b.shape} | identical: {a.equals(b)} | columns no longer written: {dropped}")
lig_new = pd.read_csv(NEW / "chronopack_c2c/chronopack_c2c.ligand_scores.csv")
lig_ref = pd.read_csv(REF / "chronopack_c2c/chronopack_c2c.ligand_scores.csv")
print("ligand ranking identical:", lig_new["ligand"].tolist() == lig_ref["ligand"].tolist(), "| top 10:", lig_new["ligand"].head(10).tolist())
for png in sorted(p.name for p in NEW.glob("*.png")):
    h = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()[:12]
    print(f"{png:44s} identical bytes to previous: {h(NEW / png) == h(REF / png)}")
print("C2C_COMPARE_OK")
