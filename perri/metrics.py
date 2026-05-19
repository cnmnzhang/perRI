"""
Evaluation metrics for Bayesian setpoint model predictions.

All three metrics use the convention that mu[t] predicts y[t+1]:
for a patient with measurements [y0, y1, ..., yN] and model outputs
[mu0, ..., muN], the evaluated pairs are (mu0→y1), ..., (mu_{N-1}→yN).

Metrics
-------
RMSE
    Measures prediction accuracy: how close the setpoint estimate is to the
    next observed value on average.

KS statistic
    Measures calibration: if the model's uncertainty (sigma) is well-calibrated,
    the Probability Integral Transform (PIT) values should be Uniform(0, 1).
    KS ≈ 0 means well-calibrated.

KS quintile range
    Measures calibration uniformity across the measurement range. Bins y_next
    into quintiles and computes KS per bin. A high range means the model is
    well-calibrated in the middle of the range but not at extremes (or vice
    versa). KS quintile range ≈ 0 is ideal.

Webapp note
-----------
Display these three numbers in a cohort-level metrics card alongside the
patient-level plots. For a single patient, only RMSE and KS are meaningful
(quintile binning requires many patients). Color-code: green if RMSE < 1 SD
of the population distribution, yellow/red otherwise.
"""

import numpy as np
import pandas as pd
from scipy.stats import kstest, norm


def compute_rmse(mus_list: list, measurements_list: list) -> float:
    """
    Root mean squared error across all patients.

    For each patient, computes (y[t+1] - mu[t])^2 and averages across all
    patients and time steps.

    Parameters
    ----------
    mus_list : list of array-like
        Posterior means per patient.
    measurements_list : list of array-like
        Observed values per patient.

    Returns
    -------
    float
    """
    assert len(mus_list) == len(measurements_list), "mus_list and measurements_list must be the same length"
    se_values = []
    for mus, meas in zip(mus_list, measurements_list):
        mus = np.asarray(mus, dtype=float)
        meas = np.asarray(meas, dtype=float)
        se = (meas[1:] - mus[:-1]) ** 2
        se_values.extend(se)
    return float(np.sqrt(np.mean(se_values)))


def compute_ks(mus_list: list, sigs_list: list, measurements_list: list) -> float:
    """
    KS statistic of pooled PIT values vs Uniform(0, 1).

    PIT value at step t: Phi((y[t+1] - mu[t]) / sigma[t])
    where Phi is the standard normal CDF.

    Returns
    -------
    float
        KS statistic. Lower is better (0 = perfectly calibrated).
    """
    pit_values = _collect_pit_values(mus_list, sigs_list, measurements_list)
    if len(pit_values) == 0:
        return float("nan")
    ks_stat, _ = kstest(pit_values, "uniform")
    return float(ks_stat)


def compute_ks_quintile_range(
    mus_list: list,
    sigs_list: list,
    measurements_list: list,
    n_bins: int = 5,
) -> float:
    """
    Max minus min KS across quintile bins of y_next.

    Measures whether calibration is uniform across the measurement range.
    A high value means the model is well-calibrated at some values but not
    others (e.g. good near the setpoint but poor at extremes).

    Returns
    -------
    float
        KS quintile range. Lower is better. np.nan if insufficient data.
    """
    y_next_all = []
    pit_all = []

    for mus, sigs, meas in zip(mus_list, sigs_list, measurements_list):
        mus = np.asarray(mus, dtype=float)
        sigs = np.asarray(sigs, dtype=float)
        meas = np.asarray(meas, dtype=float)
        if len(meas) < 2:
            continue
        y_next = meas[1:]
        mu_t = mus[:-1]
        sig_t = np.clip(sigs[:-1], 1e-6, None)
        ok = np.isfinite(y_next) & np.isfinite(mu_t) & np.isfinite(sig_t)
        if not np.any(ok):
            continue
        p = norm.cdf(y_next[ok], loc=mu_t[ok], scale=sig_t[ok])
        ok2 = np.isfinite(p)
        if not np.any(ok2):
            continue
        y_next_all.append(y_next[ok][ok2])
        pit_all.append(p[ok2])

    if not y_next_all:
        return float("nan")

    y_next_all = np.concatenate(y_next_all)
    pit_all = np.concatenate(pit_all)

    try:
        bin_ids = pd.qcut(y_next_all, q=n_bins, labels=False, duplicates="drop")
    except Exception:
        return float("nan")

    ks_vals = []
    for b in np.unique(bin_ids):
        mask = bin_ids == b
        if not np.any(mask):
            continue
        ks_vals.append(kstest(pit_all[mask], "uniform").statistic)

    if len(ks_vals) < 2:
        return float("nan")

    return float(np.nanmax(ks_vals) - np.nanmin(ks_vals))


def evaluate_batch(batch_df: pd.DataFrame) -> dict:
    """
    Compute all three metrics from a fit_batch() output DataFrame.

    Parameters
    ----------
    batch_df : pd.DataFrame
        Must have columns: patient_id, value, mu, sigma, measurement_index.
        This is the direct output of fit_batch().

    Returns
    -------
    dict with keys:
        rmse, ks, ks_quintile_range, n_patients, n_predictions
    """
    mus_list = []
    sigs_list = []
    meas_list = []

    for _, group in batch_df.sort_values("measurement_index").groupby("patient_id"):
        mus_list.append(group["mu"].values)
        sigs_list.append(group["sigma"].values)
        meas_list.append(group["value"].values)

    n_predictions = sum(max(0, len(m) - 1) for m in meas_list)

    return {
        "rmse": compute_rmse(mus_list, meas_list),
        "ks": compute_ks(mus_list, sigs_list, meas_list),
        "ks_quintile_range": compute_ks_quintile_range(mus_list, sigs_list, meas_list),
        "n_patients": len(mus_list),
        "n_predictions": n_predictions,
    }


def _collect_pit_values(mus_list, sigs_list, measurements_list) -> np.ndarray:
    """Helper: pool PIT values across all patients."""
    pit = []
    for mus, sigs, meas in zip(mus_list, sigs_list, measurements_list):
        mus = np.asarray(mus, dtype=float)
        sigs = np.asarray(sigs, dtype=float)
        meas = np.asarray(meas, dtype=float)
        if len(meas) < 2:
            continue
        sig_t = np.clip(sigs[:-1], 1e-6, None)
        p = norm.cdf(meas[1:], loc=mus[:-1], scale=sig_t)
        pit.extend(p[np.isfinite(p)].tolist())
    return np.asarray(pit)
