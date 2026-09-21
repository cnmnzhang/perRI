"""A1C mu-bound / lambda update (v0.4.0) and mu-bound override guard. Synthetic data only."""

import json
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import perri as pri
from perri import defaults

SEXES = ("ALL", "F", "M")
SNAPSHOT_V030 = json.loads((Path(__file__).parent / "fixtures" / "default_params_v0.3.0.json").read_text())

A1C_EXPECTED = {
    "ALL": {"min_mu": 4.5, "max_mu": 12.8, "log_lambda_": 1.4827167639263996},
    "F": {"min_mu": 4.5, "max_mu": 12.8, "log_lambda_": 1.5791345035395428},
    "M": {"min_mu": 4.5, "max_mu": 12.9, "log_lambda_": 1.250166972862566},
}


@pytest.mark.parametrize("sex", SEXES)
def test_a1c_new_bounds_and_lambda(sex):
    p = pri.get_default_params("A1C", sex)
    exp = A1C_EXPECTED[sex]
    assert p["min_mu"] == exp["min_mu"]
    assert p["max_mu"] == exp["max_mu"]
    assert p["log_lambda_"] == exp["log_lambda_"]


@pytest.mark.parametrize("sex", SEXES)
def test_a1c_sigma_bounds_unchanged(sex):
    new, old = pri.get_default_params("A1C", sex), SNAPSHOT_V030[f"A1C|{sex}"]
    assert new["min_sigma"] == old["min_sigma"] == 0.001
    assert new["max_sigma"] == old["max_sigma"]
    assert math.isclose(new["max_sigma"], 0.408, abs_tol=5e-4)
    assert not pri.is_log_transform("A1C")


def test_snapshot_covers_every_marker_and_sex():
    assert set(SNAPSHOT_V030) == {f"{tc}|{s}" for tc in pri.list_supported_markers() for s in SEXES}


@pytest.mark.parametrize("key", sorted(k for k in SNAPSHOT_V030 if not k.startswith("A1C|")))
def test_non_a1c_params_identical_to_v030(key):
    tc, sex = key.split("|")
    assert pri.get_default_params(tc, sex) == SNAPSHOT_V030[key]


def _synthetic_high_a1c_patient(seed=0, n=10, center=11.5):
    rng = np.random.default_rng(seed)
    values = center + rng.normal(0, 0.25, size=n)
    timestamps = pd.date_range("2015-01-01", periods=n, freq="180D")  # all isolated (>90d gaps)
    return values, timestamps


@pytest.mark.parametrize("sex", SEXES)
def test_high_a1c_patient_not_clipped_at_8(sex):
    values, timestamps = _synthetic_high_a1c_patient()
    fit = pri.fit_patient(values, timestamps, test_code="A1C", sex=sex)
    assert fit is not None
    assert fit.mu_history[-1] > 8.5
    assert abs(fit.mu_history[-1] - values[-3:].mean()) < 0.5


def test_old_bounds_would_have_clipped():
    values, timestamps = _synthetic_high_a1c_patient()
    old = dict(SNAPSHOT_V030["A1C|ALL"])
    fit = pri.fit_patient(values, timestamps, params=old)
    assert fit.mu_history[-1] <= old["max_mu"]


# ── override guard ──────────────────────────────────────────────────────────

def _override_df(**row):
    base = {"test_code": "A1C", "sex": "ALL", "min_mu": "4.5", "max_mu": "12.8"}
    base.update(row)
    return pd.DataFrame([base])


@pytest.mark.parametrize(
    "row, match",
    [
        ({"min_mu": ""}, "missing/non-numeric"),
        ({"max_mu": ""}, "missing/non-numeric"),
        ({"max_mu": "abc"}, "missing/non-numeric"),
        ({"min_mu": "nan"}, "non-finite"),
        ({"max_mu": "inf"}, "non-finite"),
        ({"min_mu": "13", "max_mu": "4"}, "min_mu < max_mu"),
        ({"test_code": "NOPE"}, "no matching row"),
    ],
)
def test_override_guard_rejects_bad_bounds(row, match):
    with pytest.raises(ValueError, match=match):
        defaults._validate_mu_bound_overrides(_override_df(**row), defaults._load())


def test_override_guard_rejects_duplicates():
    df = pd.concat([_override_df(), _override_df()], ignore_index=True)
    with pytest.raises(ValueError, match="Duplicate"):
        defaults._validate_mu_bound_overrides(df, defaults._load())


def test_bundled_overrides_are_valid():
    defaults._validate_mu_bound_overrides(pd.read_csv(defaults._DATA_DIR / "mu_bound_overrides.csv", dtype=str, keep_default_na=False), defaults._load())


def test_bad_override_file_fails_loudly_on_param_lookup(tmp_path, monkeypatch):
    shutil.copy(defaults._DATA_DIR / "bayesian_hyperparameters.csv", tmp_path)
    (tmp_path / "mu_bound_overrides.csv").write_text("test_code,sex,min_mu,max_mu\nHB,ALL,,27.75\n")
    monkeypatch.setattr(defaults, "_DATA_DIR", tmp_path)
    monkeypatch.setattr(defaults, "_params_df", None)
    monkeypatch.setattr(defaults, "_overrides_df", None)
    with pytest.raises(ValueError, match=r"\('HB', 'ALL'\)"):
        pri.get_default_params("GLU")


def test_override_applies_to_other_markers(tmp_path, monkeypatch):
    shutil.copy(defaults._DATA_DIR / "bayesian_hyperparameters.csv", tmp_path)
    (tmp_path / "mu_bound_overrides.csv").write_text("test_code,sex,min_mu,max_mu\nHB,F,6.0,20.0\n")
    monkeypatch.setattr(defaults, "_DATA_DIR", tmp_path)
    monkeypatch.setattr(defaults, "_params_df", None)
    monkeypatch.setattr(defaults, "_overrides_df", None)
    p = pri.get_default_params("HB", "F")
    assert (p["min_mu"], p["max_mu"]) == (6.0, 20.0)
    assert p["log_lambda_"] == SNAPSHOT_V030["HB|F"]["log_lambda_"]
    assert pri.get_default_params("HB", "M") == SNAPSHOT_V030["HB|M"]
