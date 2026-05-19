"""
Isolation filter for lab measurements.

A measurement is "isolated" if it is at least min_gap_days before AND after
any other measurement for the same patient. This removes clusters of repeat
draws (e.g. serial phlebotomy, same-day panels) that would bias the setpoint
estimate by over-representing a single time point.

Default gap: 90 days (matches the research repo's utils/setpoints.py).

Webapp note
-----------
For a webapp, expose min_gap_days as a UI slider so users can explore how
the isolation threshold affects the number of retained measurements and the
final setpoint estimate. A value of 0 disables filtering entirely (useful
for markers drawn very infrequently).
"""

import numpy as np
import pandas as pd


def filter_isolated(values, timestamps, min_gap_days: int = 90):
    """
    Keep only measurements that are >= min_gap_days from all neighbors.

    Parameters
    ----------
    values : array-like of float
        Measurement values.
    timestamps : array-like
        Corresponding dates/datetimes. Accepts strings parseable by pd.to_datetime.
    min_gap_days : int
        Minimum gap (in days) required before AND after a measurement to be
        considered isolated. Default 90.

    Returns
    -------
    filtered_values : np.ndarray
        Measurement values after filtering, sorted by date.
    filtered_timestamps : list
        Corresponding timestamps (as pandas Timestamps), sorted.
    """
    values = np.asarray(values, dtype=float)
    timestamps = pd.to_datetime(timestamps)

    order = np.argsort(timestamps)
    timestamps = timestamps[order]
    values = values[order]

    mask = _isolated_mask(timestamps, min_gap_days)
    return values[mask], list(timestamps[mask])


def _isolated_mask(timestamps, min_gap_days: int) -> np.ndarray:
    """Return boolean mask: True where a measurement is isolated."""
    x_days = (timestamps - timestamps[0]).days.to_numpy(dtype=float)
    front_gaps = np.diff(x_days, prepend=-np.inf)
    back_gaps = np.diff(x_days, append=np.inf)
    return (front_gaps > min_gap_days) & (back_gaps > min_gap_days)
