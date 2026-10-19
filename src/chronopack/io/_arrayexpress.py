"""Download ChronoSeq input files from ArrayExpress (EMBL-EBI BioStudies)."""

from __future__ import annotations

import json
import os
import shutil
import urllib.request
from pathlib import Path

BIOSTUDIES_INFO_URL = "https://www.ebi.ac.uk/biostudies/api/v1/studies/{accession}/info"


def _files_url(accession: str, timeout: float) -> str:
    """HTTP folder holding a study's data files, as reported by the BioStudies API."""
    with urllib.request.urlopen(BIOSTUDIES_INFO_URL.format(accession=accession), timeout=timeout) as response:
        info = json.load(response)
    if not info.get("httpLink"):
        raise RuntimeError(f"BioStudies returned no download link for {accession}: {info}")
    return info["httpLink"].rstrip("/") + "/Files"


def download_arrayexpress_files(
    accession: str,
    paths,
    *,
    overwrite: bool = False,
    timeout: float = 300.0,
) -> list[Path]:
    """Download the files in ``paths`` that are missing locally from an ArrayExpress study.

    Each file is fetched by its file name (``Path(path).name``) from the study's data files,
    e.g. ``download_arrayexpress_files("E-MTAB-15894", ["../inputs/k562-tnfa.hg38.rep/N701_...csv.gz"])``.
    Existing files are left as they are unless ``overwrite=True``, so the call costs nothing
    once the inputs are in place. Downloads go to ``<path>.<pid>.part`` first and are renamed only
    when complete (and, when the server reports it, of the expected size). Returns the paths.
    """
    paths = [Path(p) for p in ([paths] if isinstance(paths, (str, Path)) else paths)]
    missing = [p for p in paths if overwrite or not p.exists()]
    if not missing:
        return paths
    base = _files_url(accession, timeout)
    for path in missing:
        url = f"{base}/{urllib.request.quote(path.name)}"
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(f"{path.name}.{os.getpid()}.part")  # per process, so parallel runs do not collide
        print(f"Downloading {path.name} from {accession} ...", flush=True)
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response, open(part, "wb") as out:
                shutil.copyfileobj(response, out)
                expected = response.headers.get("Content-Length")
            if expected is not None and part.stat().st_size != int(expected):
                raise OSError(f"incomplete download of {url}: {part.stat().st_size} of {expected} bytes")
        except BaseException:
            part.unlink(missing_ok=True)
            raise
        part.replace(path)
    return paths
