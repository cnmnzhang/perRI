"""
Default hyperparameters for the Bayesian setpoint model.

Parameters were optimized on a clinical EHR cohort (N ~ 10,000 patients per
marker) by minimizing a composite loss of RMSE + KS + KS-quintile-range on a
held-out test set. See the research repo (models/bayesian.py, scripts/lambda_grid_search.py)
for the optimization pipeline.

The bundled CSV (data/bayesian_hyperparameters.csv) covers 43 lab
markers. Each marker has at minimum an "ALL" (sex-pooled) row and may also
have sex-stratified "M" / "F" rows when the within-person distribution differs
meaningfully by sex (e.g. hemoglobin).

Webapp note
-----------
Expose the five parameters as editable sliders so users can experiment:
  - log_lambda_  : controls recency weighting (slider range: -3 to +2)
  - min_mu / max_mu : grid bounds for the setpoint mean
  - min_sigma / max_sigma : grid bounds for within-person SD
A "Reset to defaults" button should call get_default_params() to restore the
optimized values. Changing any parameter should trigger a re-run of fit_patient().

Mu-bound overrides
------------------
By default min_mu/max_mu in the bundled CSV are derived from the population
reference interval. A marker can opt out by adding rows to
data/mu_bound_overrides.csv: (test_code, sex, min_mu, max_mu) there are used
directly in place of the base CSV's bounds, with the remaining columns recording
provenance (method, percentiles, n_measurements, data_version, source), since
perri cannot recompute them without the source data. Bounds must be in the space
the marker is fit in (log-space for log-transformed markers). Only min_mu/max_mu
are overridden; sigma bounds and log_lambda_ still come from the base CSV.
Overrides are validated on load and raise ValueError if bounds are missing,
non-finite, or not strictly increasing.
"""

import math
from pathlib import Path

from typing import Optional

import pandas as pd

_DATA_DIR = Path(__file__).parent / "data"
_PARAM_COLS = ["log_lambda_", "min_mu", "max_mu", "min_sigma", "max_sigma"]

_params_df: Optional[pd.DataFrame] = None
_overrides_df: Optional[pd.DataFrame] = None
_intra_std_df: Optional[pd.DataFrame] = None


def _load() -> pd.DataFrame:
    global _params_df
    if _params_df is None:
        path = _DATA_DIR / "bayesian_hyperparameters.csv"
        df = pd.read_csv(path, keep_default_na=False)
        df["log_transformed"] = df["log_transformed"].astype(str).str.strip().str.lower().isin({"true", "1", "yes"})
        _params_df = df
    return _params_df


def _validate_mu_bound_overrides(overrides: pd.DataFrame, params: pd.DataFrame) -> None:
    """Raise ValueError unless every override row has finite, increasing bounds for a known (test_code, sex)."""
    known = set(zip(params["test_code"], params["sex"]))
    dupes = overrides[overrides.duplicated(["test_code", "sex"], keep=False)]
    if not dupes.empty:
        raise ValueError(f"Duplicate mu-bound overrides for: {sorted(set(zip(dupes['test_code'], dupes['sex'])))}")
    for row in overrides.itertuples(index=False):
        key = (row.test_code, row.sex)
        if key not in known:
            raise ValueError(f"mu-bound override for {key} has no matching row in bayesian_hyperparameters.csv")
        try:
            lo, hi = float(row.min_mu), float(row.max_mu)
        except (TypeError, ValueError):
            raise ValueError(f"mu-bound override for {key} has missing/non-numeric bounds: min_mu={row.min_mu!r}, max_mu={row.max_mu!r}") from None
        if not (math.isfinite(lo) and math.isfinite(hi)):
            raise ValueError(f"mu-bound override for {key} has non-finite bounds: ({lo}, {hi})")
        if lo >= hi:
            raise ValueError(f"mu-bound override for {key} requires min_mu < max_mu, got ({lo}, {hi})")


def _load_mu_bound_overrides() -> pd.DataFrame:
    global _overrides_df
    if _overrides_df is None:
        path = _DATA_DIR / "mu_bound_overrides.csv"
        df = pd.read_csv(path, keep_default_na=False, dtype=str)
        _validate_mu_bound_overrides(df, _load())
        df["min_mu"] = df["min_mu"].astype(float)
        df["max_mu"] = df["max_mu"].astype(float)
        _overrides_df = df
    return _overrides_df


def _load_intra_std() -> pd.DataFrame:
    global _intra_std_df
    if _intra_std_df is None:
        path = _DATA_DIR / "marker_intra_patient_std.csv"
        _intra_std_df = pd.read_csv(path, keep_default_na=False)
    return _intra_std_df


def get_default_params(test_code: str, sex: str = "ALL") -> dict:
    """
    Return default hyperparameters for a marker.

    Parameters
    ----------
    test_code : str
        Lab marker code, e.g. "HB", "PLT", "GLU". Case-sensitive.
    sex : str
        "ALL", "M", or "F". Falls back to "ALL" if the requested sex is not
        available for this marker.

    Returns
    -------
    dict with keys: log_lambda_, min_mu, max_mu, min_sigma, max_sigma.
    min_mu/max_mu come from data/mu_bound_overrides.csv when that file has a
    row for the resolved (test_code, sex).

    Raises
    ------
    ValueError
        If test_code is not found in the bundled CSV at all.
    """
    df = _load()
    subset = df[df["test_code"] == test_code]
    if subset.empty:
        available = sorted(df["test_code"].unique().tolist())
        raise ValueError(f"test_code '{test_code}' not found. Available markers: {available}")

    row = subset[subset["sex"] == sex]
    if row.empty:
        row = subset[subset["sex"] == "ALL"]
    if row.empty:
        row = subset.iloc[[0]]

    params = row.iloc[0][_PARAM_COLS].to_dict()

    overrides = _load_mu_bound_overrides()
    override = overrides[(overrides["test_code"] == test_code) & (overrides["sex"] == row.iloc[0]["sex"])]
    if not override.empty:
        params["min_mu"] = float(override.iloc[0]["min_mu"])
        params["max_mu"] = float(override.iloc[0]["max_mu"])

    return params


def list_supported_markers() -> list:
    """Return sorted list of test_code values in the bundled parameter file."""
    df = _load()
    return sorted(df["test_code"].unique().tolist())


def is_log_transform(test_code: str) -> bool:
    """
    Return whether `test_code` is fit in log-space by default.

    Reflects the bundled CSV's `log_transformed` column -- the research repo's
    decision that log-space fitting genuinely beat raw-space fitting for that
    marker (min_mu/max_mu are bundled in log-space accordingly for such rows).
    Constant across sex rows for a given test_code. Returns False for markers
    not found in the bundled CSV.
    """
    df = _load()
    subset = df[df["test_code"] == test_code]
    if subset.empty:
        return False
    return bool(subset.iloc[0]["log_transformed"])


def get_intra_patient_std(test_code: str, sex: str = "ALL") -> float:
    """
    Return the median intra-patient standard deviation for a marker.

    Derived from the empirical distribution of per-patient within-person SDs
    on the EHR cohort used for hyperparameter optimization.

    Parameters
    ----------
    test_code : str
        Lab marker code, e.g. "HB". Case-sensitive.
    sex : str
        "ALL", "M", or "F". Falls back to "ALL" if sex-specific value is not
        available.

    Returns
    -------
    float
        Median intra-patient SD in the marker's native units.

    Raises
    ------
    ValueError
        If test_code is not found in the bundled CSV.
    """
    df = _load_intra_std()
    subset = df[df["test_code"] == test_code]
    if subset.empty:
        raise ValueError(f"test_code '{test_code}' not found in intra-patient std data.")
    row = subset[subset["sex"] == sex]
    if row.empty:
        row = subset[subset["sex"] == "ALL"]
    if row.empty:
        row = subset.iloc[[0]]
    return float(row.iloc[0]["intra_patient_std_median"])
