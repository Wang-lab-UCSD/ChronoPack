from pathlib import Path
import sys
import nbformat
from nbclient import NotebookClient

notebook_dir = Path(sys.argv[1]).resolve()
for name in ["PreprocessDatasets.ipynb", "GetTrajectories.ipynb"]:
    nb_path = notebook_dir / name
    print(f"EXECUTE {nb_path}", flush=True)
    nb = nbformat.read(nb_path, as_version=4)
    client = NotebookClient(
        nb,
        timeout=14400,
        kernel_name="python3",
        resources={"metadata": {"path": str(notebook_dir)}},
    )
    client.execute()
    out = nb_path.with_name(nb_path.stem + ".executed.ipynb")
    nbformat.write(nb, out)
    print(f"DONE {out}", flush=True)
