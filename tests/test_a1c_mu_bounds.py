"""A1C mu-bound / lambda defaults. Synthetic data only."""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import perri as pri

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


@pytest.mark.parametrize("key", sorted(k for k in SNAPSHOT_V030 if k.split("|")[0] not in {"A1C", "RBC", "RDWCV"}))
def test_other_params_identical_to_v030(key):
    tc, sex = key.split("|")
    assert pri.get_default_params(tc, sex) == SNAPSHOT_V030[key]


@pytest.mark.parametrize("test_code, sex, expected_lambda", [
    ("RBC", "ALL", 0.7841514937236451),
    ("RBC", "F", 0.43041668161665625),
    ("RBC", "M", 1.138549412002305),
    ("RDWCV", "ALL", 1.0022115712102222),
    ("RDWCV", "F", 0.9132405836660368),
    ("RDWCV", "M", 0.9716854868370356),
])
def test_rbc_rdwcv_defaults_match_updated_reference_ranges(test_code, sex, expected_lambda):
    from perri.marker_config import MARKER_CONFIG

    ri = MARKER_CONFIG[test_code]["pop_ri"]
    lower = min(bounds[0] for bounds in ri.values()) if sex == "ALL" else ri[sex][0]
    upper = max(bounds[1] for bounds in ri.values()) if sex == "ALL" else ri[sex][1]
    width = upper - lower
    midpoint = (lower + upper) / 2
    params = pri.get_default_params(test_code, sex)
    assert params["log_lambda_"] == pytest.approx(expected_lambda)
    assert params["min_mu"] == pytest.approx(midpoint - 2 * width)
    assert params["max_mu"] == pytest.approx(midpoint + 2 * width)
    assert params["max_sigma"] == pytest.approx(width / (2 * 1.96))


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
