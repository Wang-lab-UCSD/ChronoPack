import sys

from chronopack._backends.r import _r_home


def _fake_prefix(tmp_path, monkeypatch, *, with_r=True):
    if with_r:
        (tmp_path / "lib" / "R" / "bin").mkdir(parents=True)
    monkeypatch.setattr(sys, "prefix", str(tmp_path))


def test_existing_r_home_is_kept(tmp_path, monkeypatch):
    _fake_prefix(tmp_path, monkeypatch)
    monkeypatch.setenv("R_HOME", "/opt/custom/R")
    assert _r_home.ensure_r_home() == "/opt/custom/R"


def test_r_on_path_is_left_to_rpy2(tmp_path, monkeypatch):
    _fake_prefix(tmp_path, monkeypatch)
    monkeypatch.delenv("R_HOME", raising=False)
    monkeypatch.setattr(_r_home.shutil, "which", lambda name: "/usr/bin/R")
    assert _r_home.ensure_r_home() is None
    assert "R_HOME" not in _r_home.os.environ


def test_unactivated_env_uses_its_own_r(tmp_path, monkeypatch):
    _fake_prefix(tmp_path, monkeypatch)
    monkeypatch.delenv("R_HOME", raising=False)
    monkeypatch.setattr(_r_home.shutil, "which", lambda name: None)
    expected = str(tmp_path / "lib" / "R")
    assert _r_home.ensure_r_home() == expected
    assert _r_home.os.environ["R_HOME"] == expected


def test_no_r_anywhere_leaves_r_home_unset(tmp_path, monkeypatch):
    _fake_prefix(tmp_path, monkeypatch, with_r=False)
    monkeypatch.delenv("R_HOME", raising=False)
    monkeypatch.setattr(_r_home.shutil, "which", lambda name: None)
    assert _r_home.ensure_r_home() is None
    assert "R_HOME" not in _r_home.os.environ
