import io
import json

import pytest

import chronopack as cp
from chronopack.io import _arrayexpress


class _Response(io.BytesIO):
    def __init__(self, data, length=None):
        super().__init__(data)
        self.headers = {"Content-Length": str(len(data) if length is None else length)}

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def _fake_server(files, requested, truncate=False):
    def urlopen(url, timeout=None):
        requested.append(url)
        if url.endswith("/info"):
            return _Response(json.dumps({"httpLink": "https://ftp.example/E-MTAB-1"}).encode())
        name = url.rsplit("/", 1)[1]
        data = files[name]
        return _Response(data[:-1] if truncate else data, length=len(data))

    return urlopen


def test_download_arrayexpress_files_fetches_only_missing_files(tmp_path, monkeypatch):
    requested = []
    monkeypatch.setattr(_arrayexpress.urllib.request, "urlopen", _fake_server({"b.csv.gz": b"remote b"}, requested))
    present = tmp_path / "in" / "a.csv.gz"
    present.parent.mkdir()
    present.write_bytes(b"local a")

    out = cp.io.download_arrayexpress_files("E-MTAB-1", [present, tmp_path / "in" / "b.csv.gz"])

    assert [p.read_bytes() for p in out] == [b"local a", b"remote b"]
    assert requested[-1] == "https://ftp.example/E-MTAB-1/Files/b.csv.gz"
    assert not list(tmp_path.rglob("*.part"))

    requested.clear()
    cp.io.download_arrayexpress_files("E-MTAB-1", out)  # everything present: no request at all
    assert requested == []


def test_download_arrayexpress_files_discards_incomplete_download(tmp_path, monkeypatch):
    monkeypatch.setattr(
        _arrayexpress.urllib.request, "urlopen", _fake_server({"b.csv.gz": b"remote b"}, [], truncate=True)
    )
    with pytest.raises(OSError, match="incomplete download"):
        cp.io.download_arrayexpress_files("E-MTAB-1", [tmp_path / "b.csv.gz"])
    assert list(tmp_path.iterdir()) == []
