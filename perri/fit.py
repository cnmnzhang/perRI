"""
Entry points for fitting the Bayesian setpoint model.

This is the main module to import. Two functions cover all use cases:

  fit_patient(values, timestamps, test_code)
      Single patient. Returns a SetpointFit with mu_history, sigma_history,
      filtered values, and timestamps. Use this for interactive or real-time
      per-patient display.

      Example::

          result = fit_patient(values, timestamps, test_code="HB")
          # current personalized interval:
          mu, sigma = result.mu_history[-1], result.sigma_history[-1]

  fit_batch(df, value_col, timestamp_col, patient_id_col, test_code)
      Cohort of patients from a DataFrame. Returns a long-format DataFrame
      with one row per (patient, time step). Pass directly to evaluate_batch()
      for RMSE / KS metrics.

      Example::

          batch_df = fit_batch(df, "result_value", "result_date", "patient_id", test_code="HB")
          scores   = evaluate_batch(batch_df)   # {rmse, ks, ks_quintile_range, ...}

Both functions:
  - Apply the isolation filter by default (removes measurements within 90 days
    of any neighbor — prevents acute illness clusters from biasing the setpoint)
  - Load optimized hyperparameters automatically from the bundled CSV when
    test_code is given; pass params= to override
  - Return None / exclude patients with fewer than min_measurements after filtering

"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd

from .bayesian_model import bayesian
from .defaults import get_default_params
from .isolation import filter_isolated


@dataclass
class SetpointFit:
    """
    Output of fit_patient().

    Attributes
    ----------
    mu_history : np.ndarray, shape (n,)
        Posterior mean of the setpoint after each observed measurement.
        The final value mu_history[-1] is the current best estimate.
    sigma_history : np.ndarray, shape (n,)
        Posterior SD (within-person variability) after each measurement.
    values : np.ndarray, shape (n,)
        Measurement values used in fitting (after isolation filtering).
    timestamps : list of pd.Timestamp, length n
        Corresponding measurement dates.

    """

    mu_history: np.ndarray
    sigma_history: np.ndarray
    values: np.ndarray
    timestamps: list


def fit_patient(
    values,
    timestamps,
    test_code: str = None,
    sex: str = "ALL",
    params: dict = None,
    filter_isolated_measurements: bool = True,
    min_gap_days: int = 90,
    min_measurements: int = 5,
) -> SetpointFit | None:
    """
    Fit the Bayesian setpoint model for a single patient.

    Parameters
    ----------
    values : array-like of float
        Measurement values in chronological order.
    timestamps : array-like
        Corresponding dates. Accepts strings parseable by pd.to_datetime.
    test_code : str, optional
        Lab marker code (e.g. "HB", "PLT"). Required when params=None.
    sex : str
        "ALL", "M", or "F". Used to look up sex-stratified defaults.
    params : dict, optional
        Override default hyperparameters. Must include: log_lambda_, min_mu,
        max_mu, min_sigma, max_sigma. When provided, test_code is not needed
        for parameter loading (but may still be used for display/logging).
    filter_isolated_measurements : bool
        If True (default), apply the isolation filter before fitting.
    min_gap_days : int
        Gap threshold for the isolation filter. Ignored when
        filter_isolated_measurements=False.
    min_measurements : int
        Minimum number of measurements required after filtering. Returns None
        if fewer measurements remain.

    Returns
    -------
    SetpointFit or None
        None if fewer than min_measurements remain after filtering.
    """
    if params is None:
        if test_code is None:
            raise ValueError("Either test_code or params must be provided.")
        params = get_default_params(test_code, sex=sex)

    if filter_isolated_measurements:
        values, timestamps = filter_isolated(values, timestamps, min_gap_days=min_gap_days)
    else:
        values = np.asarray(values, dtype=float)
        timestamps = list(pd.to_datetime(timestamps))

    if len(values) < min_measurements:
        return None

    mu_history, sigma_history = bayesian(values, **params)
    return SetpointFit(
        mu_history=mu_history,
        sigma_history=sigma_history,
        values=values,
        timestamps=timestamps,
    )


def fit_batch(
    df: pd.DataFrame,
    value_col: str,
    timestamp_col: str,
    patient_id_col: str,
    test_code: str = None,
    sex: str = "ALL",
    params: dict = None,
    filter_isolated_measurements: bool = True,
    min_gap_days: int = 90,
    min_measurements: int = 5,
) -> pd.DataFrame:
    """
    Fit the Bayesian setpoint model for a cohort of patients.

    Parameters
    ----------
    df : pd.DataFrame
        One row per measurement. Must contain value_col, timestamp_col,
        patient_id_col columns. Additional columns (e.g. test_code, sex) are
        ignored.
    value_col : str
        Column name for the measured value.
    timestamp_col : str
        Column name for the measurement date.
    patient_id_col : str
        Column name for the patient identifier.
    test_code : str, optional
        Lab marker code. Required when params=None.
    sex : str
        "ALL", "M", or "F". Applied uniformly to all patients. For sex-stratified
        fitting, split df by sex and call fit_batch() once per sex.
    params : dict, optional
        Override hyperparameters (same keys as fit_patient()).
    filter_isolated_measurements : bool
        Apply isolation filter before fitting. Default True.
    min_gap_days : int
        Gap threshold for isolation filter.
    min_measurements : int
        Patients with fewer measurements after filtering are excluded.

    Returns
    -------
    pd.DataFrame
        Long-format output with one row per (patient, time step):
          - patient_id
          - timestamp
          - value      (observed measurement)
          - mu         (posterior mean at this step)
          - sigma      (posterior SD at this step)
          - measurement_index  (0-based index within this patient's sequence)

        Patients that do not have enough isolated measurements are silently
        excluded. Check len(result[patient_id_col].unique()) vs the input.

    """
    records = []

    for patient_id, group in df.groupby(patient_id_col):
        result = fit_patient(
            values=group[value_col].values,
            timestamps=group[timestamp_col].values,
            test_code=test_code,
            sex=sex,
            params=params,
            filter_isolated_measurements=filter_isolated_measurements,
            min_gap_days=min_gap_days,
            min_measurements=min_measurements,
        )
        if result is None:
            continue
        for i, (ts, val, mu, sigma) in enumerate(zip(result.timestamps, result.values, result.mu_history, result.sigma_history)):
            records.append(
                {
                    "patient_id": patient_id,
                    "timestamp": ts,
                    "value": val,
                    "mu": mu,
                    "sigma": sigma,
                    "measurement_index": i,
                }
            )

    if not records:
        return pd.DataFrame(columns=["patient_id", "timestamp", "value", "mu", "sigma", "measurement_index"])

    return pd.DataFrame(records)
