"""
Default hyperparameters for the Bayesian setpoint model.

Parameters were optimized on a clinical EHR cohort (N ~ 10,000 patients per
marker) by minimizing a composite loss of RMSE + KS + KS-quintile-range on a
held-out test set. See the research repo (models/bayesian.py, scripts/lambda_grid_search.py)
for the optimization pipeline.

The bundled CSV (data/bayesian_hyperparameters_v_obj2.csv) covers 43 lab
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
"""

from pathlib import Path

from typing import Optional

import pandas as pd

_DATA_DIR = Path(__file__).parent / "data"
_PARAM_COLS = ["log_lambda_", "min_mu", "max_mu", "min_sigma", "max_sigma"]

_params_df: Optional[pd.DataFrame] = None
_intra_std_df: Optional[pd.DataFrame] = None


def _load() -> pd.DataFrame:
    global _params_df
    if _params_df is None:
        path = _DATA_DIR / "bayesian_hyperparameters_v_obj2.csv"
        _params_df = pd.read_csv(path, keep_default_na=False)
    return _params_df


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
    dict with keys: log_lambda_, min_mu, max_mu, min_sigma, max_sigma

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

    return row.iloc[0][_PARAM_COLS].to_dict()


def list_supported_markers() -> list:
    """Return sorted list of test_code values in the bundled parameter file."""
    df = _load()
    return sorted(df["test_code"].unique().tolist())


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
