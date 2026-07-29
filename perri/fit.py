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


def _fit_one_patient_records(
    patient_id,
    group: pd.DataFrame,
    value_col: str,
    timestamp_col: str,
    test_code: str,
    sex: str,
    params: dict,
    filter_isolated_measurements: bool,
    min_gap_days: int,
    min_measurements: int,
) -> list:
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
        return []
    return [
        {
            "patient_id": patient_id,
            "timestamp": ts,
            "value": val,
            "mu": mu,
            "sigma": sigma,
            "measurement_index": i,
        }
        for i, (ts, val, mu, sigma) in enumerate(zip(result.timestamps, result.values, result.mu_history, result.sigma_history))
    ]


def _fit_chunk_records(
    chunk: list,
    value_col: str,
    timestamp_col: str,
    test_code: str,
    sex: str,
    params: dict,
    filter_isolated_measurements: bool,
    min_gap_days: int,
    min_measurements: int,
) -> list:
    """Fit every (patient_id, group) pair in `chunk` serially within one worker.

    Returns a list of per-patient record-lists (not flattened) -- matches
    fit_batch's serial branch shape so both can be flattened the same way.
    """
    return [
        _fit_one_patient_records(patient_id, group, value_col, timestamp_col, test_code, sex, params, filter_isolated_measurements, min_gap_days, min_measurements)
        for patient_id, group in chunk
    ]


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
    n_jobs: int = 1,
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
    n_jobs : int
        Number of parallel workers (via joblib). Default 1 (serial). Patients
        are split into n_jobs chunks, each fit serially within one worker --
        not one joblib task per patient, which under-amortizes both joblib's
        per-task dispatch overhead and numba's per-process JIT warmup for the
        (small, sub-millisecond) per-patient fit. Chunking gave ~6x wall-clock
        improvement in practice vs. ~2x for one-task-per-patient at the same
        worker count. Pass -1 to use all available cores.

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
    groups = list(df.groupby(patient_id_col))

    if n_jobs == 1 or len(groups) == 0:
        per_chunk = [
            [
                _fit_one_patient_records(
                    patient_id, group, value_col, timestamp_col, test_code, sex, params, filter_isolated_measurements, min_gap_days, min_measurements
                )
                for patient_id, group in groups
            ]
        ]
    else:
        import os

        from joblib import Parallel, delayed

        resolved_n_jobs = os.cpu_count() if n_jobs == -1 else n_jobs
        n_chunks = max(1, min(resolved_n_jobs, len(groups)))
        chunk_size, remainder = divmod(len(groups), n_chunks)
        chunks = []
        start = 0
        for i in range(n_chunks):
            size = chunk_size + (1 if i < remainder else 0)
            chunks.append(groups[start : start + size])
            start += size

        per_chunk = Parallel(n_jobs=n_chunks)(
            delayed(_fit_chunk_records)(
                chunk, value_col, timestamp_col, test_code, sex, params, filter_isolated_measurements, min_gap_days, min_measurements
            )
            for chunk in chunks
        )

    records = [record for chunk_records in per_chunk for patient_records in chunk_records for record in patient_records]

    if not records:
        return pd.DataFrame(columns=["patient_id", "timestamp", "value", "mu", "sigma", "measurement_index"])

    return pd.DataFrame(records)
