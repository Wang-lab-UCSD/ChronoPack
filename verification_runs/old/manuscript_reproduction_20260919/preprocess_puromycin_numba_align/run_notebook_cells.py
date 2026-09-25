import json
import sys
from pathlib import Path

nb_path = Path(sys.argv[1])
nb = json.loads(nb_path.read_text())
ns = {"__name__": "__main__"}
for i, cell in enumerate(nb["cells"]):
    if cell.get("cell_type") != "code":
        continue
    src = "".join(cell.get("source", []))
    if not src.strip():
        continue
    print(f"\n--- executing cell {i} ---", flush=True)
    exec(compile(src, f"{nb_path}:cell{i}", "exec"), ns)
