"""Resolve model grid bounds from a marker's clinical RI and site calibration data."""

import math

import numpy as np


def grid_from_interval(low: float, high: float, *, log_space: bool = False) -> dict:
    """Expand a finite reference interval into Bayesian mu and sigma bounds."""
    low, high = float(low), float(high)
    if not (math.isfinite(low) and math.isfinite(high) and low < high):
        raise ValueError(f"Expected a finite increasing interval, got ({low}, {high})")
    width = high - low
    center = (low + high) / 2
    min_mu = center - 2 * width
    return {
        "min_mu": min_mu if log_space else max(0.001, min_mu),
        "max_mu": center + 2 * width,
        "min_sigma": 0.001,
        "max_sigma": width / (2 * 1.96),
    }


def resolve_grid_bounds(clinical_ri: tuple, policy: dict, isolated_values=None) -> dict:
    """Apply clinical_ri, empirical_mu, or empirical_open_side policy in raw units.

    `isolated_values` is pooled across the site's complete marker population for
    the requested sex. The caller handles patient-level isolation and persistence.
    """
    method = policy.get("method", "clinical_ri")
    if method == "clinical_ri":
        return grid_from_interval(*clinical_ri)
    if method not in {"empirical_mu", "empirical_open_side"}:
        raise ValueError(f"Unknown mu-bound policy: {method}")
    percentiles = policy.get("percentiles")
    if percentiles is None or len(percentiles) != 2 or not (0 <= percentiles[0] < percentiles[1] <= 100):
        raise ValueError(f"Invalid percentiles for {method}: {percentiles}")
    if isolated_values is None:
        raise ValueError(f"{method} requires site-wide isolated measurements")
    values = np.asarray(isolated_values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2:
        raise ValueError(f"{method} requires at least two finite isolated measurements")
    empirical_low, empirical_high = map(float, np.percentile(values, percentiles))
    if not empirical_low < empirical_high:
        raise ValueError(f"{method} percentiles are not increasing: ({empirical_low}, {empirical_high})")

    if method == "empirical_mu":
        bounds = grid_from_interval(*clinical_ri)  # sigma still uses the clinical RI
        bounds.update(min_mu=empirical_low, max_mu=empirical_high)
        return bounds

    low, high = map(float, clinical_ri)
    if not math.isfinite(low):
        low = empirical_low
    if not math.isfinite(high):
        high = empirical_high
    bounds = grid_from_interval(low, high)
    # BSI estimates sigma from the expanded span when the clinical RI is open.
    bounds["max_sigma"] = (bounds["max_mu"] - bounds["min_mu"]) / 4 / (2 * 1.96)
    return bounds
