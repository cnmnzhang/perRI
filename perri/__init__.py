"""
perri
=====
Bayesian inference of individual-level laboratory reference intervals.

Quick start
-----------
    import perri as pri

    # Fit a single patient
    result = pri.fit_patient(values, timestamps, test_code="HB")
    print(result.mu_history[-1], result.sigma_history[-1])

    # Plot the fit
    fig, ax = pri.plot_fit(result.values, result.timestamps, result.mu_history, result.sigma_history)

    # Fit a cohort and evaluate
    batch_df = pri.fit_batch(df, value_col="result_value", timestamp_col="result_date", patient_id_col="patient_id", test_code="HB")
    scores = pri.evaluate_batch(batch_df)
"""

from .defaults import get_default_params, get_intra_patient_std, list_supported_markers
from .fit import SetpointFit, fit_batch, fit_patient
from .marker_config import BATTERY2TESTCODE, MARKER_CONFIG, MARKER_FULL_NAMES, MARKER_UNITS, get_population_ri
from .metrics import compute_ks, compute_ks_quintile_range, compute_rmse, evaluate_batch
from .plot import plot_fit

__all__ = [
    "fit_patient",
    "fit_batch",
    "SetpointFit",
    "plot_fit",
    "compute_rmse",
    "compute_ks",
    "compute_ks_quintile_range",
    "evaluate_batch",
    "get_default_params",
    "get_intra_patient_std",
    "list_supported_markers",
    "get_population_ri",
    "MARKER_CONFIG",
    "MARKER_UNITS",
    "MARKER_FULL_NAMES",
    "BATTERY2TESTCODE",
]
